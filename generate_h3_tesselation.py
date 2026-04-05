"""
generate_h3_tesselation.py

Create an H3 "surface" (cell coverage) over a defined extent.

Inputs
- --in     Path to input vector file (GeoJSON, Shapefile, KML, or FGDB)
- --res    H3 resolution (0..15)
- --out    Path to output vector file (format inferred from extension)

Output
- Vector file where each feature is a Polygon representing an H3 cell,
  with attributes:
    - GRID_ID: str (cell index, e.g., "8928308280fffff")
    - res: int

Requires
- pip install h3 geopandas shapely
- geo_io.py (project I/O module)

Notes
- Uses geo_io for all file I/O (supports .geojson, .shp, .kml, .gdb).
- Handles Polygon and MultiPolygon geometries.
- For multi-feature inputs, the output is the union of H3 cells across all
  polygonal features.
- Cell boundaries are emitted as Polygons in EPSG:4326 (lon/lat).
"""

import argparse
import sys
import time
from pathlib import Path
from typing import Any, Dict, Iterable, List, Set, Tuple

import geopandas as gpd
import h3  # pip install h3
from shapely.geometry import Polygon

from geo_io import read_vector, write_vector


GeoJSON = Dict[str, Any]
LonLat = Tuple[float, float]


def _validate_res(res: int) -> None:
    if not isinstance(res, int) or not (0 <= res <= 15):
        raise ValueError("H3 resolution must be an integer in [0, 15].")


def _iter_polygon_geojson_geometries(gdf: gpd.GeoDataFrame) -> Iterable[GeoJSON]:
    """
    Yield GeoJSON geometry dicts (Polygon/MultiPolygon) from a GeoDataFrame.

    Uses Shapely's __geo_interface__ directly instead of serializing the entire
    GeoSeries to a JSON string and parsing it back.
    """
    if gdf is None or gdf.empty:
        return

    for geom in gdf.geometry:
        if geom is None or geom.is_empty:
            continue
        if geom.geom_type in ("Polygon", "MultiPolygon"):
            yield geom.__geo_interface__


def _cells_from_geometry(geom: GeoJSON, res: int) -> Set[str]:
    """
    Convert a GeoJSON Polygon or MultiPolygon dict into a set of H3 cell IDs.

    Works with h3-py v4+:
      GeoJSON dict -> H3Shape via h3.geo_to_h3shape()
      H3Shape -> cells via h3.h3shape_to_cells()
    """
    gtype = geom.get("type")
    if gtype not in ("Polygon", "MultiPolygon"):
        return set()

    shape = h3.geo_to_h3shape(geom)
    return set(h3.h3shape_to_cells(shape, res=res))


def _cell_to_polygon(cell_id: str) -> Polygon:
    """
    Convert an H3 cell ID to a Shapely Polygon.

    h3.cell_to_boundary returns (lat, lon); Shapely needs (lon, lat).
    """
    boundary_latlon = h3.cell_to_boundary(cell_id)
    ring_lonlat = [(lon, lat) for (lat, lon) in boundary_latlon]
    return Polygon(ring_lonlat)


def build_h3_surface(gdf: gpd.GeoDataFrame, res: int) -> gpd.GeoDataFrame:
    """
    Build an H3 surface GeoDataFrame from an in-memory GeoDataFrame extent.

    Parameters
    ----------
    gdf : GeoDataFrame
        Input extent polygons (must be in EPSG:4326 or will be reprojected).
    res : int
        H3 resolution (0..15).

    Returns
    -------
    GeoDataFrame
        One row per H3 cell with columns: GRID_ID, res, geometry.
        CRS is EPSG:4326.
    """
    _validate_res(res)

    # Ensure lon/lat for H3
    if gdf.crs is not None:
        try:
            epsg = gdf.crs.to_epsg()
        except Exception:
            epsg = None
        if epsg != 4326:
            gdf = gdf.to_crs(epsg=4326)

    geoms = list(_iter_polygon_geojson_geometries(gdf))
    if not geoms:
        raise ValueError(
            "No Polygon/MultiPolygon geometries found in input dataset. "
            "Provide polygonal extent features."
        )

    # Collect all cell IDs across all input polygons
    cells: Set[str] = set()
    for geom in geoms:
        cells.update(_cells_from_geometry(geom, res))

    if not cells:
        raise ValueError(
            f"No H3 cells generated at resolution {res}. "
            "Check that the input extent is valid and non-degenerate."
        )

    # Build GeoDataFrame directly from cell IDs
    sorted_cells = sorted(cells)
    polygons = [_cell_to_polygon(cid) for cid in sorted_cells]

    return gpd.GeoDataFrame(
        {
            "GRID_ID": sorted_cells,
            "res": [res] * len(sorted_cells),
        },
        geometry=polygons,
        crs="EPSG:4326",
    )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Create an H3 surface (cells as polygons) from a vector extent. "
            "Supports GeoJSON, Shapefile, KML, and FGDB input/output."
        )
    )
    p.add_argument(
        "--in",
        dest="in_path",
        required=True,
        help="Path to input vector file (.geojson, .shp, .kml, .gdb).",
    )
    p.add_argument(
        "--res",
        dest="res",
        required=True,
        type=int,
        help="H3 resolution (0..15).",
    )
    p.add_argument(
        "--out",
        dest="out_path",
        required=True,
        help="Path to output vector file (format inferred from extension).",
    )
    p.add_argument(
        "--layer",
        dest="layer",
        default=None,
        help="Layer name for multi-layer inputs (FGDB, KML).",
    )
    p.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing output file.",
    )
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    start = time.perf_counter()

    in_path = Path(args.in_path)
    out_path = Path(args.out_path)

    try:
        gdf = read_vector(in_path, layer=args.layer, epsg=4326)
        print(f"[OK] Read {len(gdf)} features from: {in_path}")

        surface = build_h3_surface(gdf, args.res)

        write_vector(surface, out_path, overwrite=args.overwrite)

        elapsed = time.perf_counter() - start
        print(f"[OK] Wrote {len(surface)} H3 cells to: {out_path}")
        print(f"[OK] Completed in {elapsed:.3f}s")
        return 0

    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())