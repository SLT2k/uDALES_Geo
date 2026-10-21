"""Automatic alignment to a reference building raster (align.py) on synthetic
layouts: references are cut with the production convention
domain = translate(rotate(world, A, origin=c), -w0) at known A and w0, in
both turning directions. Offline; 4 m cells keep it fast."""
import json
import math
import os
import random

import numpy as np
import pytest
from rasterio.features import rasterize
from rasterio.transform import from_origin
from shapely.affinity import rotate, scale, translate
from shapely.geometry import box, mapping

from udales_geo import align
from udales_geo.common import ConfigError

X0, Y0 = 500000.0, 5700000.0
CRS = "EPSG:32630"
RES, W, H, R = 4.0, 240.0, 200.0, 260.0        # cell, reference size, search half-width
STEP, REFINE = 6.0, 1.0


def layout(seed=7, n=14, pitch=50.0):
    """An L-shaped district of irregular blocks (no symmetry) round (X0, Y0)."""
    rng = random.Random(seed)
    out = []
    for i in range(n):
        for j in range(n):
            if i > 0.7 * n and j > 0.5 * n:
                continue                           # the missing quadrant makes an L
            if rng.random() < 0.25:
                continue
            w, h = rng.uniform(8, 30), rng.uniform(8, 22)
            x = -n * pitch / 2 + i * pitch + rng.uniform(0, 12)
            y = -n * pitch / 2 + j * pitch + rng.uniform(0, 20)
            p = box(X0 + x, Y0 + y, X0 + x + w, Y0 + y + h)
            if rng.random() < 0.3:
                p = p.difference(box(X0 + x, Y0 + y, X0 + x + w / 2, Y0 + y + h / 2))
            out.append(p)
    return out


def cut_reference(geoms, A, c, w0, res=RES, w=W, h=H):
    """The reference raster of the window at w0 after turning the map by A
    about c (s4 convention), row 0 = top."""
    turned = [translate(rotate(p, A, origin=c), -w0[0], -w0[1]) for p in geoms]
    return rasterize(((p, 1) for p in turned), out_shape=(int(h / res), int(w / res)),
                     transform=from_origin(0, h, res, res), fill=0, dtype="uint8") > 0


GEOMS = layout()
C = (X0 + 10.0, Y0 - 5.0)
W0 = (C[0] - 118.7, C[1] - 97.3)               # off the raster lattice on purpose


def corners(frame):
    xs, ys = frame.xy_to_world([0.0, W, W, 0.0], [0.0, 0.0, H, H])
    return np.column_stack([xs, ys])


def corner_error(res, A=23.0, w0=W0):
    """Largest corner distance [m] between the frame found and the true one."""
    found = align.frame_from_rotation(res["angle"], C, res["window_origin_rotated"])
    true = align.frame_from_rotation(A, C, w0)
    return float(np.hypot(*(corners(found) - corners(true)).T).max())


def acfg_for(mask, **kw):
    cfg = dict(align.DEFAULTS)
    cfg.update(min_score=0.5, mirror_check=False, resolution_m=RES, first_row="top", mask=mask)
    cfg.update(kw)
    return cfg


# --------------------------------------------------------------------------
# the search
# --------------------------------------------------------------------------
@pytest.mark.parametrize("A", [0.0, 23.0, 337.0, 135.0])
def test_angle_and_window_are_recovered_in_both_turning_directions(A):
    ref = cut_reference(GEOMS, A, C, W0)
    res, curves = align.register(GEOMS, ref, RES, C, R, STEP, REFINE, mirror_check=False)
    found = align.frame_from_rotation(res["angle"], C, res["window_origin_rotated"])
    # the sign is pinned by the corners: frame_from_rotation defines x_axis_deg = -A
    assert found.x_axis_deg == pytest.approx((-A) % 360.0, abs=REFINE / 2 + 1e-9)
    assert corner_error(res, A) <= RES                          # every corner within a cell
    assert res["score"] > 0.8
    acc, ratios, reasons, warnings = align.judge(res, acfg_for(ref))
    assert acc == "score" and not reasons and not warnings


def test_prominence_runner_up_and_second_position():
    ref = cut_reference(GEOMS, 23.0, C, W0)
    res, curves = align.register(GEOMS, ref, RES, C, R, STEP, REFINE, mirror_check=False)
    s = np.asarray(curves["coarse"])
    assert res["prominence"] == pytest.approx(s.max() / np.median(s))
    ru = res["runner_up"]
    assert align._circular_gap(ru["angle"], res["coarse"]["angle"]) > align.RUNNER_UP_GAP_DEG
    assert ru["score"] < s.max()
    i = curves["angles"].index(ru["angle"])                    # a peak of its own
    assert s[i] >= s[i - 1] and s[i] >= s[(i + 1) % len(s)]
    assert res["second_position"]["score"] < res["score"]
    k_max = math.ceil(0.5 * STEP / REFINE) + 1                 # half a coarse step + one refine step
    assert res["refine_evaluations"] == (2 * k_max + 1) * 4 - 1


def test_a_weak_upside_down_match_is_refused_and_names_first_row():
    ref = np.flipud(cut_reference(GEOMS, 23.0, C, W0))         # a wrong first_row
    res, _ = align.register(GEOMS, ref, RES, C, R, STEP, REFINE, mirror_check=True)
    acc, ratios, reasons, warnings = align.judge(res, acfg_for(ref))
    assert acc is None                                         # refused by score and prominence
    assert res["mirror"]["score"] > 0.8 > res["score"]
    assert any('"first_row": "bottom"' in r for r in reasons)


def test_reference_of_another_place_is_refused():
    other = layout(seed=99)
    ref = cut_reference(other, 40.0, C, W0)
    res, _ = align.register(GEOMS, ref, RES, C, R, STEP, REFINE, mirror_check=False)
    acc, _, reasons, _ = align.judge(res, acfg_for(ref))
    assert acc is None and "best score" in reasons[0]


def test_a_regular_grid_is_ambiguous():
    grid = [box(X0 + x, Y0 + y, X0 + x + 20, Y0 + y + 20)
            for x in range(-400, 400, 40) for y in range(-400, 400, 40)]
    ref = cut_reference(grid, 0.0, C, (C[0] - 120.0, C[1] - 100.0))
    res, _ = align.register(grid, ref, RES, C, R, STEP, 0.0, mirror_check=False)
    acc, ratios, reasons, warnings = align.judge(res, acfg_for(ref))
    assert ratios["position"] > 0.85
    assert acc == "score" and not reasons                     # accepted by its score: a warning
    assert any("second window position" in w for w in warnings)
    as_peak = _result(0.45, 3.0, position=ratios["position"])   # the same, passing by prominence
    acc, _, reasons, _ = align.judge(as_peak, acfg_for(ref))
    assert acc is None and any("peak prominence only" in r for r in reasons)
    acc, _, reasons, warnings = align.judge(as_peak, acfg_for(ref, max_ambiguity=1.0))
    assert acc == "peak" and not reasons and not warnings      # 1: never refused


def _result(score, prominence, angle=None, position=None, mirror=None, scored=1.0):
    """A register() result with the given score and ambiguity ratios."""
    return {"angle": 337.0, "score": score, "prominence": prominence,
            "coarse": {"angle": 337.0, "score": score}, "scored_fraction": scored,
            "runner_up": None if angle is None else {"angle": 157.0, "score": angle * score},
            "second_position": {"score": None if position is None else position * score,
                                "min_distance_cells": 10},
            "mirror": None if mirror is None else {"angle": 157.0, "score": mirror * score}}


@pytest.mark.parametrize("res, kw, accepted, n_reasons, warn_text", [
    # au_can_v2 (left-right symmetric): right match, upside down 0.93 x -> kept, warned
    (_result(0.91, 3.5, angle=0.81, position=0.58, mirror=0.93), {}, "score", 0,
     '"first_row": "bottom"'),
    # ca_van_v3 read with the wrong first_row: 0.371, prominence 3.41, the right way up
    # 2.33 x better -> accepted by prominence WITH the warning (the documented trade-off)
    (_result(0.371, 3.41, mirror=2.33), {}, "peak", 0, "probably upside down"),
    (_result(0.371, 3.41, mirror=2.33), {"first_row": "bottom"}, "peak", 0, '"first_row": "top"'),
    (_result(0.9, 4.0, angle=0.9), {}, "score", 0, "second angle"),
    (_result(0.9, 4.0, position=0.9), {}, "score", 0, "second window position"),
    (_result(0.35, 4.0, angle=0.9), {}, None, 1, None),          # prominence only: refused
    (_result(0.35, 4.0, position=0.9), {}, None, 1, None),
    (_result(0.35, 4.0, angle=0.95, position=0.97), {"max_ambiguity": 1.0}, "peak", 0, None),
    (_result(0.62, 4.0, mirror=1.17), {"max_ambiguity": 1.0}, "score", 0, "UPSIDE DOWN"),
    (_result(0.9, 4.0, mirror=0.95), {"max_ambiguity": 2.0}, "score", 0, None),
    (_result(0.9, 4.0, scored=0.5), {}, None, 1, None),
    (_result(0.31, 1.55, mirror=2.7), {}, None, 2, None),       # weak, and upside down better
])
def test_judge_warns_or_refuses(res, kw, accepted, n_reasons, warn_text):
    acc, ratios, reasons, warnings = align.judge(res, acfg_for(None, **kw))
    assert acc == accepted and len(reasons) == n_reasons
    if warn_text is None:
        assert not warnings
    else:
        assert len(warnings) == 1 and warn_text in warnings[0]
    if n_reasons == 2:
        assert '"first_row": "bottom"' in reasons[1]


@pytest.mark.parametrize("subset", ["60 % at random", "4 % at random", "a far ring only",
                                    "one quarter only"])
def test_reregistration_on_fewer_footprints_finds_the_same_frame_or_refuses(subset):
    """The 785 lesson: production registered a second time with only the
    buildings that had a height, and that reduced set gave -176 deg at score
    0.17. Fewer footprints than the reference shows must give the same frame or
    be refused, never a wrong frame that passes."""
    ref = cut_reference(GEOMS, 23.0, C, W0)
    rng = random.Random(5)
    pick = {"60 % at random": lambda p: rng.random() < 0.6,
            "4 % at random": lambda p: rng.random() < 0.04,
            "a far ring only": lambda p: max(abs(p.centroid.x - C[0]),
                                             abs(p.centroid.y - C[1])) > 150,
            "one quarter only": lambda p: p.centroid.x < C[0] - 30 and p.centroid.y < C[1] - 30}
    few = [p for p in GEOMS if pick[subset](p)]
    res, _ = align.register(few, ref, RES, C, R, STEP, REFINE, mirror_check=False)
    acc, _, reasons, _ = align.judge(res, acfg_for(ref))
    err = corner_error(res)
    assert acc is None or err <= RES                           # accepted => the right frame
    if subset == "60 % at random":
        assert acc == "score" and err <= RES
    else:                                                       # wrong angle or window: refused
        assert err > 10 * RES and acc is None and reasons


def test_workers_give_the_same_result():
    small = [p for p in GEOMS if abs(p.centroid.x - C[0]) < 200 and abs(p.centroid.y - C[1]) < 200]
    ref = cut_reference(small, 23.0, C, W0)
    one, c1 = align.register(small, ref, RES, C, R, 30.0, 0.0, mirror_check=False, workers=1)
    two, c2 = align.register(small, ref, RES, C, R, 30.0, 0.0, mirror_check=False, workers=2)
    assert one == two and c1 == c2


def test_workers_in_a_script_without_main_guard_stop_with_an_error(tmp_path):
    """Each worker of the spawn pool re-runs the calling script; without an
    `if __name__ == "__main__":` guard it fails there. The pool must then raise
    with the advice, not start new workers forever (multiprocessing.Pool did)."""
    import subprocess
    import sys
    import textwrap
    from pathlib import Path
    repo = str(Path(__file__).resolve().parents[1])
    script = tmp_path / "no_guard.py"
    script.write_text(textwrap.dedent(f"""
        import sys
        sys.path.insert(0, {repo!r})
        import numpy as np
        from shapely.geometry import box
        from udales_geo import align
        ref = np.zeros((20, 20), dtype=bool)
        ref[5:12, 4:9] = True
        try:
            align._run([(0.0, (0.0, 0.0), False)] * 3, [box(0, 0, 8, 10)], ref, None, 1.0,
                       (0.0, 0.0), 30.0, 2)
        except RuntimeError as e:
            print("RAISED:", e)
            sys.exit(3)
        print("NO ERROR")
        """), encoding="utf-8")
    env = {**os.environ, "PYTHONDONTWRITEBYTECODE": "1"}
    proc = subprocess.run([sys.executable, str(script)], capture_output=True, text=True,
                          timeout=300, cwd=tmp_path, env=env)
    assert proc.returncode == 3, proc.stdout + proc.stderr
    assert "worker processes of the search stopped" in proc.stdout
    assert "__main__" in proc.stdout


# --------------------------------------------------------------------------
# the reference file and the config block
# --------------------------------------------------------------------------
MASK = np.zeros((60, 70), dtype=int)
MASK[5:20, 10:30] = 12
MASK[40:55, 50:60] = 7


def test_reference_formats_and_row_order(tmp_path):
    np.savetxt(tmp_path / "a_topo", MASK, fmt="%d")
    np.savetxt(tmp_path / "a.csv", MASK, fmt="%d", delimiter=",")
    np.save(tmp_path / "a.npy", MASK)
    for name in ("a_topo", "a.csv", "a.npy"):
        m = align.read_reference(tmp_path / name)
        assert m.dtype == bool and (m == (MASK > 0)).all()
    assert (align.read_reference(tmp_path / "a.npy", "bottom") == np.flipud(MASK > 0)).all()


@pytest.mark.parametrize("arr, match", [
    (np.where(MASK > 0, np.nan, 0.0), "NaN"),
    (np.zeros((60, 70)), "no building"),
    (np.ones((60, 70)), "no ground"),
    (np.zeros(70), "2-D")])
def test_bad_reference_is_refused(tmp_path, arr, match):
    np.save(tmp_path / "r.npy", arr)
    with pytest.raises(ConfigError, match=match):
        align.read_reference(tmp_path / "r.npy")


def test_text_that_is_no_table_is_refused(tmp_path):
    (tmp_path / "t_topo").write_text("1 2 3\n4 x 6\n")
    with pytest.raises(ConfigError, match="2-D table of numbers"):
        align.read_reference(tmp_path / "t_topo")


NEAR = [-0.1721, 51.5016]


@pytest.mark.parametrize("p3d, match", [
    (" nx = 69,ny = 59,\n dx = 1, dy = 1,\n", None),                       # matches
    (" nx = 69,ny = 59,\n dx = 0.5, dy = 0.5,\n", 'dx = 0.5: set "resolution_m": 0.5'),
    (" nx = 69,ny = 59,\n dy = 0.5,\n", 'dy = 0.5: set "resolution_m": 0.5'),
    (" nx = 383,ny = 59,\n dx = 1, dy = 1,\n", "70 columns but X_p3d says nx = 383"),
    (" nx = 69,ny = 431,\n dx = 1, dy = 1,\n", "60 rows but X_p3d says ny = 431")])
def test_palm_p3d_next_to_the_topo_is_checked(tmp_path, p3d, match):
    np.savetxt(tmp_path / "X_topo", MASK, fmt="%d")
    (tmp_path / "X_p3d").write_text("&initialization_parameters\n" + p3d + "/\n")
    if match is None:
        cfg = align.parse_align({"reference": "X_topo", "near": NEAR}, tmp_path)
        assert (cfg["W"], cfg["H"]) == (70.0, 60.0)
    else:
        with pytest.raises(ConfigError, match=match):
            align.parse_align({"reference": "X_topo", "near": NEAR}, tmp_path)
    (tmp_path / "X_p3d").unlink()                               # no _p3d: no check
    assert align.parse_align({"reference": "X_topo", "near": NEAR}, tmp_path)["W"] == 70.0


@pytest.mark.parametrize("extra, match", [
    ({"centre": NEAR}, 'unknown key "centre"'),
    ({"near": None}, 'missing: "near"'),
    ({"near": [151.2, -233.9]}, "lat-first"),
    ({"reference": 5}, "reference must be a file path"),
    ({"first_row": "up"}, "first_row"),
    ({"resolution_m": 0}, "resolution_m must be > 0"),
    ({"angle_step_deg": 120}, "angle_step_deg"),
    ({"refine_step_deg": 2.0}, "refine_step_deg"),
    ({"min_score": 1.5}, "min_score"),
    ({"min_prominence": 0.5}, "min_prominence"),
    ({"max_ambiguity": 0}, "max_ambiguity"),
    ({"mirror_check": "yes"}, "mirror_check"),
    ({"workers": 0}, "workers"),
    ({"workers": 1.5}, "workers"),
    ({"search_radius_m": 40}, "smaller than half the reference diagonal"),
    ({"max_memory_gb": 0.5}, "memory"),
    ({"frame": {"x_axis_deg": 30, "origin": [X0, Y0], "crs": CRS}}, '"near". Remove them')])
def test_bad_align_block_is_refused(tmp_path, extra, match):
    np.save(tmp_path / "r.npy", MASK)
    with pytest.raises(ConfigError, match=match):
        align.parse_align({"reference": "r.npy", "near": NEAR, **extra}, tmp_path)


FRAME = {"x_axis_deg": 30.0, "origin": [X0, Y0], "crs": CRS}


@pytest.mark.parametrize("extra, match", [
    ({"frame": {"x_axis_deg": 30}}, '"align" frame must be'),
    ({"frame": {**FRAME, "angle": 3}}, '"align" frame must be'),
    ({"frame": {**FRAME, "origin": [X0, "y"]}}, "must be numbers"),
    ({"frame": {**FRAME, "crs": 32630}}, "crs must be a CRS code"),
    ({"frame": {**FRAME, "window_score": 1.5}}, "window_score must be"),
    ({"frame": {**FRAME, "fetch_box": [1, 2, 0, 3]}}, "fetch_box must be"),
    ({"frame": {**FRAME, "footprints_md5": "abc"}}, "footprints_md5 must be"),
    ({"frame": FRAME, "search_radius_m": -5, "workers": 2},
     '"search_radius_m", "workers". Remove them'),
    ({"frame": FRAME, "max_ambiguity": 0.9}, '"max_ambiguity". Remove them')])
def test_bad_frame_is_refused(tmp_path, extra, match):
    np.save(tmp_path / "r.npy", MASK)
    with pytest.raises(ConfigError, match=match):
        align.parse_align({"reference": "r.npy", **extra}, tmp_path)


def test_small_or_missing_reference_is_refused(tmp_path):
    np.save(tmp_path / "s.npy", MASK[:40, :30])
    with pytest.raises(ConfigError, match="only 30 x 40 m"):
        align.parse_align({"reference": "s.npy", "near": NEAR}, tmp_path)
    with pytest.raises(ConfigError, match="reference not found"):
        align.parse_align({"reference": "nope.npy", "near": NEAR}, tmp_path)
    with pytest.raises(ConfigError, match='needs "reference"'):
        align.parse_align({"near": NEAR}, tmp_path)


def test_defaults_and_crs(tmp_path):
    np.save(tmp_path / "r.npy", MASK)
    cfg = align.parse_align({"reference": "r.npy", "near": NEAR, "refine_step_deg": None},
                            tmp_path)
    assert cfg["refine_step_deg"] == 0.1                       # null = the default
    assert cfg["search_radius_m"] == pytest.approx(math.hypot(70, 60) + 100.0)
    assert align.resolve_crs(cfg, "auto") == CRS
    assert align.resolve_crs(cfg, "EPSG:27700") == "EPSG:27700"
    assert align.parse_align({"reference": "r.npy", "near": NEAR, "max_ambiguity": 1.5},
                             tmp_path)["max_ambiguity"] == 1.5
    cfg = align.parse_align({"reference": "r.npy", "frame": FRAME}, tmp_path)   # no "near"
    assert align.resolve_crs(cfg, "auto") == CRS
    assert cfg["frame_record"] == {"window_score": None, "fetch_box": None, "footprints_md5": None}
    with pytest.raises(ConfigError, match="frame is in EPSG:32630"):
        align.resolve_crs(cfg, "EPSG:27700")


def test_angles_are_normalised_below_360():
    from udales_geo.region import make_rotated_region, norm_deg, parse_domain
    assert norm_deg(-1e-20) == 0.0 and (-1e-20) % 360.0 == 360.0   # the float trap
    assert norm_deg(-30.0) == 330.0 and norm_deg(720.5) == 0.5
    centre, size, a = parse_domain({"centre": NEAR, "size_m": [100, 80], "x_axis_deg": -1e-20})
    assert a == 0.0 and make_rotated_region(centre, size, a, crs=CRS).frame is None
    assert align.frame_from_rotation(1e-20, C, W0).x_axis_deg == 0.0
    assert align.frame_from_rotation(-360.0, C, W0).x_axis_deg == 0.0


# --------------------------------------------------------------------------
# stage 3b: align_site (cache, acceptance, warnings, window check, freeze)
# --------------------------------------------------------------------------
def _site(tmp_path, A=23.0, geoms=GEOMS, shift=(0.0, 0.0), w0=None, **kw):
    """A reference file cut at angle A round the projected 'near' (window
    centre about there, moved by `shift` in the turned map), the footprints
    as a GeoJSON file, and the parsed align block."""
    from pyproj import Transformer
    near = Transformer.from_crs(CRS, "EPSG:4326", always_xy=True).transform(*C)
    c = Transformer.from_crs("EPSG:4326", CRS, always_xy=True).transform(*near)
    if w0 is None:
        w0 = (c[0] - 118.7 + shift[0], c[1] - 97.3 + shift[1])
    np.save(tmp_path / "ref.npy", cut_reference(geoms, A, c, w0))
    feats = [{"type": "Feature", "properties": {"id": i}, "geometry": mapping(p)}
             for i, p in enumerate(geoms)]
    gj = tmp_path / "buildings_osm.geojson"
    gj.write_text(json.dumps({"type": "FeatureCollection", "features": feats}), encoding="utf-8")
    raw = {"reference": "ref.npy", "near": list(near), "resolution_m": RES,
           "search_radius_m": R, "angle_step_deg": STEP, "refine_step_deg": REFINE,
           "mirror_check": False, **kw}
    return align.parse_align(raw, tmp_path), feats, gj, raw


def test_align_site_caches_and_rejudges(tmp_path, capsys):
    acfg, feats, gj, raw = _site(tmp_path)
    region, rep = align.align_site(acfg, feats, gj, CRS, tmp_path, 1.0, 4.0)
    assert (tmp_path / "alignment.png").exists() and (tmp_path / "alignment.json").exists()
    assert rep["x_axis_deg"] == pytest.approx(337.0, abs=0.5) and rep["accepted_by"] == "score"
    assert rep["warnings"] == [] and rep["window_score"] == pytest.approx(rep["score"], abs=1e-3)
    assert region.bbox == (0.0, 0.0, 240.0, 200.0) and region.frame is not None
    assert rep["window_step_m"] == RES / 2
    assert "angles x" in capsys.readouterr().out
    cache = json.loads((tmp_path / "alignment.json").read_text(encoding="utf-8"))
    assert cache["judgement"]["accepted_by"] == "score" and cache["judgement"]["warnings"] == []
    region2, rep2 = align.align_site(acfg, feats, gj, CRS, tmp_path, 1.0, 4.0)   # reused
    out = capsys.readouterr().out
    assert "reusing alignment.json" in out and "angles x" not in out
    assert region2.frame == region.frame and rep2 == rep
    strict = align.parse_align({**raw, "min_score": 1.0, "min_prominence": 99.0}, tmp_path)
    with pytest.raises(ConfigError, match="do not match the reference"):        # re-judged only
        align.align_site(strict, feats, gj, CRS, tmp_path, 1.0, 4.0)
    assert "angles x" not in capsys.readouterr().out
    cache = json.loads((tmp_path / "alignment.json").read_text(encoding="utf-8"))
    assert cache["judgement"]["accepted_by"] is None and cache["judgement"]["refused"]
    gj.write_text(gj.read_text() + " ", encoding="utf-8")       # a hand edit: new footprints md5
    align.align_site(acfg, feats, gj, CRS, tmp_path, 1.0, 4.0)
    assert "angles x" in capsys.readouterr().out
    align.align_site(acfg, feats, gj, CRS, tmp_path, 5.0, 4.0)  # min_area is part of the key
    assert "angles x" in capsys.readouterr().out


@pytest.mark.parametrize("change", ["near", "search_radius_m", "angle_step_deg",
                                    "refine_step_deg", "mirror_check", "reference"])
def test_each_search_input_is_part_of_the_cache_key(tmp_path, capsys, change):
    # a cheap search (12 angles); acceptance does not matter here
    acfg, feats, gj, raw = _site(tmp_path, A=330.0, angle_step_deg=30.0, refine_step_deg=0.0)

    def run(raw_block):
        try:
            align.align_site(align.parse_align(raw_block, tmp_path), feats, gj, CRS, tmp_path,
                             1.0, 4.0)
        except ConfigError:
            pass
        return "angles x" in capsys.readouterr().out
    assert run(raw) and not run(raw)                            # searched, then reused
    if change == "reference":
        m = np.load(tmp_path / "ref.npy")
        m[0, 0] = not m[0, 0]
        np.save(tmp_path / "ref.npy", m)
        changed = raw
    else:
        changed = {**raw, change: {"near": [raw["near"][0] + 1e-5, raw["near"][1]],
                                   "search_radius_m": R + 4.0, "angle_step_deg": 45.0,
                                   "refine_step_deg": 5.0, "mirror_check": True}[change]}
    assert run(changed)


def test_a_symmetric_site_passes_with_an_upside_down_warning(tmp_path, capsys):
    """The au_can_v2 lesson: a site that is left-right symmetric scores as
    well with the reference upside down. That must not refuse the right
    match; it warns everywhere a user looks."""
    left = [p for p in GEOMS if p.bounds[2] < C[0] - 2.0]
    sym = left + [scale(p, xfact=-1.0, yfact=1.0, origin=(C[0], 0.0)) for p in left]
    acfg, feats, gj, raw = _site(tmp_path, A=0.0, geoms=sym, w0=(C[0] - W / 2, C[1] - 100.0),
                                 mirror_check=True)
    region, rep = align.align_site(acfg, feats, gj, CRS, tmp_path, 1.0, 4.0)
    assert rep["accepted_by"] == "score" and rep["mirror"]["ratio"] > 0.85
    assert len(rep["warnings"]) == 1 and "upside down" in rep["warnings"][0].lower()
    assert '"first_row": "bottom"' in rep["warnings"][0]
    cache = json.loads((tmp_path / "alignment.json").read_text(encoding="utf-8"))
    assert cache["judgement"]["warnings"] == rep["warnings"]
    out = capsys.readouterr().out
    assert "[align] warning: the reference upside" in out.lower() and "!!!!!!!!" in out
    assert (tmp_path / "alignment.png").stat().st_size > 0
    a = region.frame.x_axis_deg
    assert min(a, 360.0 - a) <= 0.6                            # the right way up, not 180 deg


def test_freeze_block_gives_the_same_domain_without_a_search(tmp_path, capsys):
    acfg, feats, gj, raw = _site(tmp_path)
    region, rep = align.align_site(acfg, feats, gj, CRS, tmp_path, 1.0, 4.0)
    frozen = json.loads(json.dumps(rep["freeze"]))              # as pasted from the report
    fr = frozen["frame"]
    assert fr["crs"] == CRS and frozen["resolution_m"] == RES
    assert fr["footprints_md5"] == align.md5_of(gj) and fr["fetch_box"] == rep["search_box"]
    assert fr["window_score"] == pytest.approx(rep["window_score"], abs=1e-4)
    capsys.readouterr()
    fcfg = align.parse_align(frozen, tmp_path)
    region2, rep2 = align.align_site(fcfg, feats, gj, CRS, tmp_path, 1.0, 4.0)
    assert "angles x" not in capsys.readouterr().out
    assert region2.frame == region.frame and region2.bbox == region.bbox
    assert rep2["mode"] == "frame" and rep2["iou_in_window"] == rep["iou_in_window"]
    assert rep2["score"] == pytest.approx(fr["window_score"], abs=1e-4) and rep2["warnings"] == []


def test_frozen_frame_checks_the_score_and_the_download(tmp_path):
    acfg, feats, gj, raw = _site(tmp_path)
    _, rep = align.align_site(acfg, feats, gj, CRS, tmp_path, 1.0, 4.0)
    frozen = json.loads(json.dumps(rep["freeze"]))
    # another download than the one the frame was found with: a warning, the frame is kept
    gj.write_text(gj.read_text() + "\n", encoding="utf-8")
    region, rep2 = align.align_site(align.parse_align(frozen, tmp_path), feats, gj, CRS,
                                    tmp_path, 1.0, 4.0)
    assert len(rep2["warnings"]) == 1 and "not the download this frame" in rep2["warnings"][0]
    # the frame pasted elsewhere (40 m off): its window scores far below the record -> stop
    moved = json.loads(json.dumps(frozen))
    moved["frame"]["origin"][0] += 40.0
    with pytest.raises(ConfigError, match="far below"):
        align.align_site(align.parse_align(moved, tmp_path), feats, gj, CRS, tmp_path, 1.0, 4.0)
    # without the recorded score only min_score is left to compare with: a warning
    del moved["frame"]["window_score"]
    _, rep3 = align.align_site(align.parse_align(moved, tmp_path), feats, gj, CRS, tmp_path,
                               1.0, 4.0)
    assert rep3["score"] < 0.5 and any("no window_score" in w for w in rep3["warnings"])
    # no footprint in the window at all
    far = [{"type": "Feature", "properties": {},
            "geometry": mapping(box(X0 + 5000, Y0, X0 + 5010, Y0 + 10))}]
    with pytest.raises(ConfigError, match="no building footprint lies in the frozen window"):
        align.align_site(align.parse_align(frozen, tmp_path), far, gj, CRS, tmp_path, 1.0, 4.0)


def test_window_outside_the_search_square_is_refused(tmp_path):
    # radius just over half the diagonal (156 m), the window centre 25 m off "near"
    acfg, feats, gj, raw = _site(tmp_path, shift=(25.0, 0.0), search_radius_m=160.0)
    with pytest.raises(ConfigError, match="outside the search square"):
        align.align_site(acfg, feats, gj, CRS, tmp_path, 1.0, 4.0)


def test_no_footprint_near_is_a_clear_error(tmp_path):
    acfg, feats, gj, raw = _site(tmp_path)
    far = [{"type": "Feature", "properties": {},
            "geometry": mapping(box(X0 + 5000, Y0, X0 + 5010, Y0 + 10))}]
    with pytest.raises(ConfigError, match="no building footprint within"):
        align.align_site(acfg, far, gj, CRS, tmp_path, 1.0, 4.0)


# --------------------------------------------------------------------------
# build_site with "align" (download and mesher stubbed: conftest.offline_build)
# --------------------------------------------------------------------------
def _align_config(tmp_path, align_block, **extra):
    site = {"name": "al", "crs": CRS, "snap_m": 4.0, "dx": 4.0, "output_dir": str(tmp_path / "site"),
            "align": align_block, "edge_buildings": "extend", "edge_clearance_m": 4.0,
            "buildings": {"height_source": "osm"}, "vegetation": {"enabled": False}, **extra}
    p = tmp_path / f"al_{len(list(tmp_path.glob('al_*.jsonc')))}.jsonc"
    p.write_text(json.dumps(site), encoding="utf-8")
    return p


def _search_block(tmp_path):
    acfg, feats, gj, raw = _site(tmp_path)
    gj.unlink()                                                 # build_site downloads its own
    return feats, {**raw, "reference": str(tmp_path / "ref.npy")}


def test_build_site_align_search_wiring(tmp_path, offline_build, capsys):
    from udales_geo.build import build_site
    offline_build.features, block = _search_block(tmp_path)
    rep = build_site(_align_config(tmp_path, block))
    (box_, crs), = offline_build.fetches                        # one download: the search square
    sq = rep["align"]["search_box"]
    assert crs == CRS and list(box_) == pytest.approx(sq, abs=1e-3)
    rec = json.loads((tmp_path / "site" / "buildings_osm.fetch.json").read_text())
    assert rec == {"crs": CRS, "bbox": [round(v, 3) for v in box_]}
    cfg = offline_build.cfgs[0]
    frame = cfg.raw["buildings"]["_frame"]
    assert frame.x_axis_deg == pytest.approx(337.0, abs=0.5)
    assert rep["region"]["frame"]["origin_world"] == [frame.x0, frame.y0]
    assert rep["edge_buildings"]["lattice_anchor"] == [0.0, 0.0]   # align: reference lattice
    assert rep["align"]["accepted_by"] == "score" and rep["align"]["warnings"] == []
    assert rep["fetch"]["n_in_search_square"] == len(offline_build.features)
    assert 0 < rep["fetch"]["n_buildings"] < rep["fetch"]["n_in_search_square"]
    out = capsys.readouterr().out
    assert "fetched for the search" in out and "WARNING   :" not in out


def test_build_site_align_counts_the_first_tile_only(tmp_path, offline_build):
    from udales_geo.build import build_site
    offline_build.features, block = _search_block(tmp_path)
    one = build_site(_align_config(tmp_path, block))
    two = build_site(_align_config(tmp_path, block, tile=[2, 1]))   # reuses download and search
    assert len(offline_build.fetches) == 1
    assert two["fetch"]["n_buildings"] == one["fetch"]["n_buildings"]
    assert two["region"]["xlen"] == 2 * one["region"]["xlen"]
    first = two["region"]["frame"]["first_tile"]
    assert first["corners_world"] == one["region"]["frame"]["corners_world"]
    assert "first_tile" not in one["region"]["frame"]


def test_build_site_freeze_keeps_the_download(tmp_path, offline_build, capsys):
    from udales_geo.build import build_site
    offline_build.features, block = _search_block(tmp_path)
    rep = build_site(_align_config(tmp_path, block))
    site_dir = tmp_path / "site"
    frozen = rep["align"]["freeze"]
    offline_build.allow_fetch = False                           # a download fails the test
    capsys.readouterr()
    rep2 = build_site(_align_config(tmp_path, frozen))
    assert rep2["align"]["mode"] == "frame" and rep2["align"]["warnings"] == []
    assert rep2["region"]["frame"]["origin_world"] == rep["region"]["frame"]["origin_world"]
    assert rep2["region"]["bbox"] == rep["region"]["bbox"]
    assert not (site_dir / "buildings_osm.previous.geojson").exists()
    assert "angles x" not in capsys.readouterr().out
    # a hand edit (another md5): kept, with a warning in the report and the summary
    gj = site_dir / "buildings_osm.geojson"
    gj.write_text(gj.read_text(encoding="utf-8") + "\n", encoding="utf-8")
    rep3 = build_site(_align_config(tmp_path, frozen))
    assert any("not the download this frame" in w for w in rep3["align"]["warnings"])
    assert "WARNING   : [align]" in capsys.readouterr().out
    assert not (site_dir / "buildings_osm.previous.geojson").exists()
    # a download for another area: refused, nothing moved, nothing downloaded
    rec = site_dir / "buildings_osm.fetch.json"
    rec.write_text(json.dumps({"crs": CRS, "bbox": [X0 + 2000, Y0, X0 + 2100, Y0 + 100]}))
    with pytest.raises(ConfigError, match="does not cover the frozen"):
        build_site(_align_config(tmp_path, frozen))
    assert not (site_dir / "buildings_osm.previous.geojson").exists()
    # --force-fetch downloads the search's box again (recorded in the frame)
    offline_build.allow_fetch = True
    build_site(_align_config(tmp_path, frozen), force_fetch=True)
    (box_, _), = offline_build.fetches[1:]
    assert list(box_) == frozen["frame"]["fetch_box"]
    # no record at all: the file is used as it is
    rec.unlink()
    offline_build.allow_fetch = False
    rep4 = build_site(_align_config(tmp_path, frozen))
    assert rep4["align"]["mode"] == "frame"
