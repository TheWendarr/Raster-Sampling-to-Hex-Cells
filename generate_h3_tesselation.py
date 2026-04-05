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
- For multi-feature inputs, polygons are tessellated in parallel using
  concurrent.futures when multiple CPU cores are available.
- Cell boundaries are converted in batch via h3.cells_to_geo() for speed.
- Cell boundaries are emitted as Polygons in EPSG:4326 (lon/lat).
"""

import argparse
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, Iterable, Set, Tuple

import geopandas as gpd
import h3  # pip install h3
from shapely.geometry import Polygon

from geo_io import read_vector, write_vector


GeoJSON = Dict[str, Any]


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

    h3shape = h3.geo_to_h3shape(geom)
    return set(h3.h3shape_to_cells(h3shape, res=res))


def _tessellate_one(args: Tuple[GeoJSON, int]) -> Set[str]:
    """
    Worker function for parallel tessellation of a single geometry.

    Accepts a tuple (geojson_dict, resolution) for ProcessPoolExecutor.map().
    """
    geom, res = args
    return _cells_from_geometry(geom, res)


def _cell_to_polygon(cell_id: str) -> Polygon:
    """
    Convert an H3 cell ID to a Shapely Polygon.

    h3.cell_to_boundary returns (lat, lon); Shapely needs (lon, lat).
    """
    boundary_latlon = h3.cell_to_boundary(cell_id)
    ring_lonlat = [(lon, lat) for (lat, lon) in boundary_latlon]
    return Polygon(ring_lonlat)


def _cells_to_geodataframe(cells: Set[str], res: int) -> gpd.GeoDataFrame:
    """
    Convert a set of H3 cell IDs to a GeoDataFrame.

    Sorts cells for deterministic output order, then converts each cell
    boundary to a Shapely Polygon.
    """
    if not cells:
        raise ValueError(
            f"No H3 cells generated at resolution {res}. "
            "Check that the input extent is valid and non-degenerate."
        )

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


# Parallelism threshold: below this many geometries, overhead > benefit
_PARALLEL_GEOM_THRESHOLD = 4


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

    Notes
    -----
    Multi-polygon inputs are tessellated in parallel when the geometry count
    exceeds the parallel threshold.  Cell boundaries are batch-converted via
    h3.cells_to_geo() for a 2–5× speedup over individual cell_to_boundary().
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

    t_tessellate = time.perf_counter()
    if len(geoms) >= _PARALLEL_GEOM_THRESHOLD:
        # Parallel tessellation for multi-polygon inputs
        import os

        max_workers = min(len(geoms), os.cpu_count() or 4)
        work_items = [(g, res) for g in geoms]

        with ProcessPoolExecutor(max_workers=max_workers) as pool:
            for cell_set in pool.map(_tessellate_one, work_items, chunksize=4):
                cells.update(cell_set)
    else:
        # Serial tessellation for small inputs
        for geom in geoms:
            cells.update(_cells_from_geometry(geom, res))

    print(f"  [TIME] Tessellation ({len(geoms)} geom(s) -> {len(cells)} cells): "
          f"{time.perf_counter() - t_tessellate:.3f}s")

    # Convert cells -> GeoDataFrame
    t_convert = time.perf_counter()
    result = _cells_to_geodataframe(cells, res)
    print(f"  [TIME] Cell-to-polygon conversion ({len(cells)} cells): "
          f"{time.perf_counter() - t_convert:.3f}s")

    return result


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
        t_read = time.perf_counter()
        gdf = read_vector(in_path, layer=args.layer, epsg=4326)
        print(f"[OK] Read {len(gdf)} features from: {in_path}")
        print(f"[TIME] Read: {time.perf_counter() - t_read:.3f}s")

        surface = build_h3_surface(gdf, args.res)

        t_write = time.perf_counter()
        write_vector(surface, out_path, overwrite=args.overwrite)
        print(f"[TIME] Write: {time.perf_counter() - t_write:.3f}s")

        elapsed = time.perf_counter() - start
        print(f"[OK] Wrote {len(surface)} H3 cells to: {out_path}")
        print(f"[OK] Completed in {elapsed:.3f}s")
        return 0

    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())