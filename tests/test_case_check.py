"""Post-udprep checker (case_check.py) on synthetic IBM masks. Offline."""
import numpy as np

from udales_geo.case_check import narrow_cavities


def test_one_cell_shaft_is_found_and_wide_courtyard_is_not():
    solid = np.zeros((20, 20, 6), dtype=bool)
    solid[5:12, 5:12, :4] = True          # block, 4 levels
    solid[8, 8, :4] = False               # one-cell shaft inside it
    solid[12:19, 2:18, :3] = True         # second block with a 5 x 5 courtyard
    solid[13:18, 5:10, :3] = False
    comp, cells, cols = narrow_cavities(solid)
    assert cols == 1 and comp == 4 and cells == 4      # the shaft, on each of 4 levels


def test_pocket_at_the_edge_joins_the_street_across_the_periodic_seam():
    solid = np.zeros((20, 20, 2), dtype=bool)
    solid[0:3, 4:9, :] = True             # buildings against the x = 0 edge ...
    solid[0:2, 6, :] = False              # ... with a slot that is open only across the seam
    solid[1, 4:9, :] = True
    solid[0, 5:8, :] = True
    solid[0, 6, :] = False
    comp, _, _ = narrow_cavities(solid)
    assert comp == 0


def test_y_seam_and_periodic_block_test():
    # a slot at the y = 0 edge, open only across the y seam
    solid = np.zeros((20, 20, 2), dtype=bool)
    solid[4:9, 0:3, :] = True
    solid[6, 0:2, :] = False
    solid[4:9, 1, :] = True
    solid[5:8, 0, :] = True
    solid[6, 0, :] = False
    assert narrow_cavities(solid)[0] == 0
    # an enclosed pocket that is 3 x 3 only when wrapped across the x seam: not narrow
    solid = np.ones((20, 20, 1), dtype=bool)
    solid[:, 0:8, 0] = False                       # open area (the street network)
    solid[18:20, 14:17, 0] = False                 # pocket rows 18, 19 ...
    solid[0:1, 14:17, 0] = False                   # ... and 0: 3 x 3 across the seam
    assert narrow_cavities(solid)[0] == 0


def test_non_periodic_boundary_counts_as_open():
    # a slot reaching the x = 0 boundary: open for inflow/outflow (BCxm != 1)
    solid = np.zeros((20, 20, 2), dtype=bool)
    solid[0:4, 5:10, :] = True
    solid[0:3, 7, :] = False
    solid[-1, 5:10, :] = True
    assert narrow_cavities(solid, periodic=(True, True))[0] > 0
    assert narrow_cavities(solid, periodic=(False, True))[0] == 0


def test_solid_mask_from_another_grid_is_unusable_input(tmp_path):
    # namoptions changed after udprep ran: exit code 2, not a crash or a pass
    from udales_geo.case_check import main
    case = tmp_path / "007"
    case.mkdir()
    (case / "namoptions.007").write_text("&DOMAIN\nitot = 8\njtot = 8\nktot = 4\n/\n")
    (case / "solid_c.txt").write_text("# i j k\n   1    1    9\n")
    assert main([str(case)]) == 2
    (case / "solid_c.txt").write_text("# i j k\n   1    1    1\n")
    assert main([str(case)]) == 0
