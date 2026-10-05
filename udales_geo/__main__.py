"""CLI: python -m udales_geo configs/<site>.jsonc [--force-fetch] [--redo-trees]"""
from __future__ import annotations

import argparse
import sys

from .build import build_site
from .common import ConfigError


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("config", help="site config JSON (see configs/ for examples)")
    ap.add_argument("--force-fetch", action="store_true",
                    help="download the OSM buildings again and rebuild the trees; "
                         "overwrites hand edits of buildings_osm.geojson / trees.geojson "
                         "(LiDAR and canopy tiles already in cache/ are reused)")
    ap.add_argument("--redo-trees", action="store_true",
                    help="rebuild chm.tif and trees.geojson only (keeps buildings)")
    args = ap.parse_args()
    try:
        build_site(args.config, force_fetch=args.force_fetch, redo_trees=args.redo_trees)
    except ConfigError as e:
        print(f"\nConfig problem in {args.config}:\n{e}", file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    main()
