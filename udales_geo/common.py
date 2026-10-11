"""Shared config / grid / frame helpers: the JSON-with-comments reader
(load_jsonc), ConfigError, the uDALES tools lookup, Grid and Config.

Frames:
  world  : the site's projected CRS (auto = UTM zone of the site), metres.
  sim    : uDALES frame [0,xlen] x [0,ylen], ground at z=0;
           sim = world - (bbox xmin, bbox ymin). The domain corner lies on the
           1 cm snap lattice (region.on_lattice).

Relative paths in a config are resolved against the config file's folder.
"""
from __future__ import annotations

import json
import sys
from dataclasses import dataclass, field
from pathlib import Path

def resolve_udales_tools(explicit=None) -> Path:
    """Locate uDALES's python tools (udgeom) for the geometry QA gate.
    Order: explicit config value > UDALES_TOOLS environment variable.
    No hardcoded machine paths — every machine states its own clone."""
    import os
    candidates = []
    if explicit:
        candidates.append(Path(explicit))
    env = os.environ.get("UDALES_TOOLS")
    if env:
        candidates.append(Path(env))
    for c in candidates:
        if c.exists():
            return c.resolve()
    raise ConfigError(
        "uDALES python tools (udgeom) not found. Set \"udales_tools\" in the "
        "config, or the UDALES_TOOLS environment variable, to "
        "<your u-dales clone>/tools/python."
        + (f" Tried: {[str(c) for c in candidates]}" if candidates else ""))


@dataclass
class Grid:
    itot: int
    jtot: int
    ktot: int
    xlen: float
    ylen: float
    zsize: float
    # optional explicit z cell-centre list (stretched grids); overrides uniform zsize/ktot
    zt_list: list | None = None

    @property
    def dx(self) -> float:
        return self.xlen / self.itot

    @property
    def dy(self) -> float:
        return self.ylen / self.jtot

    def zt(self):
        import numpy as np
        if self.zt_list is not None:
            zt = np.asarray(self.zt_list, dtype=float)
            if len(zt) != self.ktot:
                raise ValueError(f"zt_list has {len(zt)} entries, ktot={self.ktot}")
            return zt
        dz = self.zsize / self.ktot
        return (np.arange(self.ktot) + 0.5) * dz

    def dzt(self):
        import numpy as np
        zt = self.zt()
        if self.zt_list is not None and zt[-1] >= self.zsize:
            raise ValueError("zt_list top cell centre must lie below zsize")
        # cell edges: 0, midpoints, zsize — mirrors udgrid.z_grid_from_profile,
        # which is what udprep's sim.dzt uses (top edge is always zsize)
        zm = np.concatenate([[0.0], 0.5 * (zt[:-1] + zt[1:]), [self.zsize]])
        return np.diff(zm)

    def xt(self):
        import numpy as np
        return (np.arange(self.itot) + 0.5) * self.dx

    def yt(self):
        import numpy as np
        return (np.arange(self.jtot) + 0.5) * self.dy


@dataclass
class Config:
    path: Path                      # config file location (anchor for relative paths)
    raw: dict = field(repr=False, default_factory=dict)

    name: str = "unnamed"
    crs: str = "EPSG:27700"
    bbox: tuple = (0.0, 0.0, 0.0, 0.0)   # xmin, ymin, xmax, ymax in `crs`
    grid: Grid | None = None
    udales_tools: Path | None = None

    @property
    def origin(self) -> tuple:
        return (self.bbox[0], self.bbox[1])

    def resolve(self, p: str | Path) -> Path:
        p = Path(p)
        return p if p.is_absolute() else (self.path.parent / p).resolve()

    def section(self, name: str) -> dict:
        return self.raw.get(name, {})


def load_config(config_file: str | Path) -> Config:
    config_file = Path(config_file).resolve()
    raw = json.loads(config_file.read_text(encoding="utf-8-sig"))
    grid = Grid(**raw["grid"]) if "grid" in raw else None
    tools = raw.get("udales_tools")
    if tools is not None:
        tools = Path(tools)
        if not tools.is_absolute():
            tools = (config_file.parent / tools).resolve()
    tools = resolve_udales_tools(tools)
    cfg = Config(
        path=config_file,
        raw=raw,
        name=raw.get("name", config_file.stem),
        crs=raw.get("crs", "EPSG:27700"),
        bbox=tuple(raw["bbox"]),
        grid=grid,
        udales_tools=tools,
    )
    if grid is not None:
        bx = cfg.bbox[2] - cfg.bbox[0]
        by = cfg.bbox[3] - cfg.bbox[1]
        if abs(bx - grid.xlen) > 0.5 or abs(by - grid.ylen) > 0.5:
            print(f"[common] WARNING: bbox size ({bx:.1f} x {by:.1f}) != grid xlen/ylen "
                  f"({grid.xlen} x {grid.ylen}); grid values are authoritative, "
                  f"bbox min is the origin.")
    return cfg


def add_udgeom_to_path(cfg: Config) -> None:
    """Make the uDALES python tools (udgeom, udprep) importable."""
    tools = str(cfg.udales_tools)
    if tools not in sys.path:
        sys.path.insert(0, tools)


def world_to_sim_polygon(poly, origin):
    from shapely.affinity import translate
    return translate(poly, xoff=-origin[0], yoff=-origin[1])


class ConfigError(ValueError):
    """A problem in the user's config file; the CLI prints it without a traceback."""


def load_jsonc(path: str | Path) -> dict:
    """Read a JSON config that may contain // line and /* block */ comments
    (so every parameter can carry its explanation). Comment markers inside
    "quoted strings" (URLs, paths) are left alone."""
    text = Path(path).read_text(encoding="utf-8-sig")
    out, i, n, in_str = [], 0, len(text), False
    while i < n:
        c = text[i]
        if in_str:
            out.append(c)
            if c == "\\" and i + 1 < n:
                out.append(text[i + 1])
                i += 1
            elif c == '"':
                in_str = False
        elif c == '"':
            in_str = True
            out.append(c)
        elif text.startswith("//", i):
            while i < n and text[i] != "\n":
                i += 1
            continue
        elif text.startswith("/*", i):
            end = text.find("*/", i + 2)
            stop = n if end < 0 else end + 2
            out.append("\n" * text.count("\n", i, stop))   # keep error line numbers right
            i = stop
            continue
        else:
            out.append(c)
        i += 1
    try:
        return json.loads("".join(out))
    except json.JSONDecodeError as e:
        hint = ("write paths with / instead of \\ (e.g. \"C:/data/tile.laz\")"
                if "escape" in e.msg else "a missing or extra comma or quote?")
        raise ConfigError(f"{Path(path).name}, line {e.lineno}: {e.msg} — {hint}") from None


def write_report(path: Path, report: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, default=str), encoding="utf-8")
    print(f"[report] {path}")
