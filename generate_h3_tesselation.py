"""
h3_surface_from_geojson.py

Create an H3 "surface" (cell coverage) over a GeoJSON-defined extent.

Inputs
- --in     Path to input GeoJSON (Geometry | Feature | FeatureCollection)
- --res    H3 resolution (0..15)
- --out    Path to output GeoJSON (FeatureCollection of H3 cell polygons)

Output
- GeoJSON FeatureCollection where each feature is a Polygon representing an H3 cell,
  with properties:
    - h3_id: str (cell index, e.g., "8928308280fffff")
    - res: int

Requires
- pip install h3

Notes
- Handles Polygon and MultiPolygon (and Feature/FeatureCollection wrappers).
- For FeatureCollection, the output is the union of H3 cells across all polygonal features.
- Cell boundaries are emitted as GeoJSON Polygons (lon/lat), closed rings.
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Set, Tuple, Union

import h3  # pip install h3


GeoJSON = Dict[str, Any]
LonLat = Tuple[float, float]


def _load_geojson(path: Path) -> GeoJSON:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception as e:
        raise RuntimeError(f"Failed to read/parse GeoJSON: {path}\n{e}") from e


def _iter_polygon_geometries(obj: GeoJSON) -> Iterable[GeoJSON]:
    """
    Yield Polygon/MultiPolygon geometries from a GeoJSON object that may be:
    - Geometry
    - Feature
    - FeatureCollection
    """
    t = obj.get("type")
    if t == "FeatureCollection":
        for feat in obj.get("features", []) or []:
            if not isinstance(feat, dict):
                continue
            geom = feat.get("geometry")
            if isinstance(geom, dict):
                yield from _iter_polygon_geometries(geom)
        return

    if t == "Feature":
        geom = obj.get("geometry")
        if isinstance(geom, dict):
            yield from _iter_polygon_geometries(geom)
        return

    # Geometry
    if t in ("Polygon", "MultiPolygon"):
        yield obj
        return

    # Ignore non-polygon geometries (Point/LineString/etc.)
    return


def _validate_res(res: int) -> None:
    if not isinstance(res, int) or not (0 <= res <= 15):
        raise ValueError("H3 resolution must be an integer in [0, 15].")


def _cells_from_geometry(geom: dict, res: int) -> set[str]:
    """
    Convert a GeoJSON Polygon or MultiPolygon dict into a set of H3 cell IDs (strings).

    Works with h3-py v4+ by:
      GeoJSON dict -> H3Shape via h3.geo_to_h3shape()
      H3Shape -> cells via h3.h3shape_to_cells()
    """
    gtype = geom.get("type")
    if gtype not in ("Polygon", "MultiPolygon"):
        return set()

    shape = h3.geo_to_h3shape(geom)
    return set(h3.h3shape_to_cells(shape, res=res))



def _close_ring(ring: List[LonLat]) -> List[LonLat]:
    if not ring:
        return ring
    if ring[0] != ring[-1]:
        ring = ring + [ring[0]]
    return ring


def _cell_to_geojson_feature(cell_id: str) -> GeoJSON:
    """
    Convert an H3 cell to a GeoJSON Polygon Feature with properties.
    h3.cell_to_boundary returns (lat, lon); GeoJSON needs (lon, lat).
    """
    boundary_latlon = h3.cell_to_boundary(cell_id)  # list[(lat, lon)]
    ring_lonlat: List[LonLat] = [(lon, lat) for (lat, lon) in boundary_latlon]
    ring_lonlat = _close_ring(ring_lonlat)

    return {
        "type": "Feature",
        "properties": {"GRID_ID": cell_id, "res": h3.get_resolution(cell_id)},
        "geometry": {"type": "Polygon", "coordinates": [ring_lonlat]},
    }


def build_h3_surface(in_geojson: GeoJSON, res: int) -> GeoJSON:
    _validate_res(res)

    geoms = list(_iter_polygon_geometries(in_geojson))
    if not geoms:
        raise ValueError(
            "No Polygon/MultiPolygon geometries found in input GeoJSON. "
            "Provide an extent as Polygon/MultiPolygon (Feature/FeatureCollection ok)."
        )

    cells: Set[str] = set()
    for geom in geoms:
        cells.update(_cells_from_geometry(geom, res))

    features = [_cell_to_geojson_feature(cid) for cid in sorted(cells)]

    return {"type": "FeatureCollection", "features": features}


def parse_args(argv: List[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Create an H3 surface (cells as polygons) from a GeoJSON extent."
    )
    p.add_argument("--in", dest="in_path", required=True, help="Path to input GeoJSON.")
    p.add_argument("--res", dest="res", required=True, type=int, help="H3 resolution (0..15).")
    p.add_argument("--out", dest="out_path", required=True, help="Path to output GeoJSON.")
    p.add_argument(
        "--indent",
        dest="indent",
        type=int,
        default=2,
        help="JSON indentation level for output (default: 2).",
    )
    return p.parse_args(argv)


def main(argv: List[str] | None = None) -> int:
    args = parse_args(argv)
    in_path = Path(args.in_path).expanduser().resolve()
    out_path = Path(args.out_path).expanduser().resolve()

    if not in_path.exists():
        print(f"[ERROR] Input file not found: {in_path}", file=sys.stderr)
        return 2

    try:
        in_geojson = _load_geojson(in_path)
        out_geojson = build_h3_surface(in_geojson, args.res)

        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(out_geojson, indent=args.indent), encoding="utf-8")

        print(f"[OK] Wrote {len(out_geojson['features'])} H3 cells to: {out_path}")
        return 0

    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
