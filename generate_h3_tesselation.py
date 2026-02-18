"""
generate_h3_tesselation.py

Create an H3 "surface" (cell coverage) over a defined extent.

Inputs
- --in     Path to input vector file (GeoJSON, Shapefile, or KML)
- --res    H3 resolution (0..15)
- --out    Path to output GeoJSON (FeatureCollection of H3 cell polygons)

Output
- GeoJSON FeatureCollection where each feature is a Polygon representing an H3 cell,
  with properties:
    - GRID_ID: str (cell index, e.g., "8928308280fffff")
    - res: int

Requires
- pip install h3 geopandas pyogrio

Notes
- Uses GeoPandas to load input as a GeoDataFrame (supports .geojson/.shp/.kml depending on GDAL drivers).
- Handles Polygon and MultiPolygon geometries.
- For multi-feature inputs, the output is the union of H3 cells across all polygonal features.
- Cell boundaries are emitted as GeoJSON Polygons (lon/lat), closed rings.
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Set, Tuple, Union

import geopandas as gpd  # pip install geopandas pyogrio
import h3  # pip install h3


GeoJSON = Dict[str, Any]
LonLat = Tuple[float, float]


def _validate_res(res: int) -> None:
    if not isinstance(res, int) or not (0 <= res <= 15):
        raise ValueError("H3 resolution must be an integer in [0, 15].")


def _iter_polygon_geojson_geometries_from_gdf(gdf: "gpd.GeoDataFrame") -> Iterable[GeoJSON]:
    """
    Yield GeoJSON geometry dicts (Polygon/MultiPolygon) from a GeoDataFrame.

    - Reprojects to EPSG:4326 if CRS is defined and not already 4326 (H3 expects lon/lat degrees).
    - Ignores non-polygon geometries.
    """
    if gdf is None:
        return

    if gdf.empty:
        return

    if "geometry" not in gdf.columns:
        return

    # Ensure lon/lat (EPSG:4326) for H3
    if gdf.crs is not None:
        try:
            epsg = gdf.crs.to_epsg()
        except Exception:
            epsg = None

        if epsg != 4326:
            gdf = gdf.to_crs(epsg=4326)

    # GeoSeries.to_json() returns a GeoJSON FeatureCollection (string)
    # We'll parse and iterate geometries from there.
    geojson_fc = json.loads(gdf.geometry.to_json())
    for feat in geojson_fc.get("features", []) or []:
        geom = feat.get("geometry")
        if not isinstance(geom, dict):
            continue
        gtype = geom.get("type")
        if gtype in ("Polygon", "MultiPolygon"):
            yield geom


def _cells_from_geometry(geom: GeoJSON, res: int) -> Set[str]:
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


def build_h3_surface(gdf: "gpd.GeoDataFrame", res: int) -> GeoJSON:
    """
    Build an H3 surface from an in-memory GeoDataFrame.
    """
    _validate_res(res)

    geoms = list(_iter_polygon_geojson_geometries_from_gdf(gdf))
    if not geoms:
        raise ValueError(
            "No Polygon/MultiPolygon geometries found in input dataset. "
            "Provide polygonal extent features."
        )

    cells: Set[str] = set()
    for geom in geoms:
        cells.update(_cells_from_geometry(geom, res))

    features = [_cell_to_geojson_feature(cid) for cid in sorted(cells)]

    return {"type": "FeatureCollection", "features": features}


def parse_args(argv: List[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Create an H3 surface (cells as polygons) from a vector extent (GeoJSON/Shapefile/KML)."
    )
    p.add_argument(
        "--in",
        dest="in_path",
        required=True,
        help="Path to input vector file (.geojson/.shp/.kml).",
    )
    p.add_argument("--res", dest="res", required=True, type=int, help="H3 resolution (0..15).")
    p.add_argument("--out", dest="out_path", required=True, help="Path to output GeoJSON.")
    p.add_argument(
        "--indent",
        dest="indent",
        type=int,
        default=2,
        help="JSON indentation level for output (default: 2).",
    )
    p.add_argument(
        "--layer",
        dest="layer",
        default=None,
        help="Optional: layer name (sometimes needed for KML or multi-layer sources).",
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
        # Load as GeoDataFrame (in-memory)
        gdf = gpd.read_file(in_path, layer=args.layer) if args.layer else gpd.read_file(in_path)

        out_geojson = build_h3_surface(gdf, args.res)

        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(out_geojson, indent=args.indent), encoding="utf-8")

        print(f"[OK] Wrote {len(out_geojson['features'])} H3 cells to: {out_path}")
        return 0

    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
