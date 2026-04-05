"""
build_surface.py

Manage Hex Surface Model FGDBs — the authoritative geospatial data fabric
built from H3 tessellations attributed with raster-derived terrain factors.

Subcommands:
init
    Create a new FGDB tagged as a Hex Surface Model.

commit
    Tessellate raster footprints into H3 cells, sample raster values,
    and write/merge into the target Surface Model FGDB.

    Supports two modes:
      1. Raster-only (original): tessellate from raster footprints
         --raster path:method:resolution

      2. Pre-built surface: skip tessellation, sample into provided cells
         --surface path.geojson --raster path:method[:resolution]

info
    Display model metadata, LOD summary, and commit history.

Usage:
    # Create a new surface model
    python build_surface.py init --target "path/to/model.gdb"

    # Commit rasters (tessellate from footprints)
    python build_surface.py commit --target "path/to/model.gdb" \\
        --raster path/to/elevation.tif:mean:9 \\
        --raster path/to/landcover.tif:majority:9

    # Commit rasters using a pre-built H3 surface
    python build_surface.py commit --target "path/to/model.gdb" \\
        --surface path/to/h3_cells.geojson \\
        --raster path/to/elevation.tif:mean \\
        --raster path/to/landcover.tif:majority

    # View model info
    python build_surface.py info --target "path/to/model.gdb"

Requires:
    geo_io.py, generate_h3_tesselation.py, sample_raster.py
    geopandas, pyogrio, rasterio, shapely, h3, numpy
"""


import argparse
import datetime
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import geopandas as gpd
import numpy as np
import rasterio
from shapely.geometry import Point, box

from geo_io import list_layers, read_vector
from generate_h3_tesselation import build_h3_surface
from sample_raster import RasterBinEnricher

try:
    import pyogrio  # type: ignore

    _HAS_PYOGRIO = True
except ImportError:
    pyogrio = None  # type: ignore
    _HAS_PYOGRIO = False


# Constants

META_TABLE = "__hex_surface_meta__"
LOG_TABLE = "__hex_surface_log__"
MAGIC_ID = "HEX_SURFACE_MODEL"
SCHEMA_VERSION = 1
LOD_PREFIX = "LOD_"

# Columns that belong to the H3 grid structure, not raster data
SYSTEM_COLUMNS = frozenset({"GRID_ID", "res", "geometry"})

# Columns injected by RasterBinEnricher that we strip before writing to the
# LOD feature class — provenance is tracked in the commit log instead.
META_COLUMN_PREFIX = "meta_"


# Data Structures

@dataclass(frozen=True)
class RasterConfig:
    """A single raster input with its sampling method and target resolution."""

    path: Path
    method: str  # "mean" or "majority"
    resolution: Optional[int]  # H3 resolution (0..15), None when using --surface


# LOD Naming

def lod_name(res: int) -> str:
    """H3 resolution -> LOD feature class name (e.g., 9 -> 'LOD_09')."""
    return f"{LOD_PREFIX}{res:02d}"


def parse_lod_name(name: str) -> Optional[int]:
    """LOD feature class name -> H3 resolution, or None if not a LOD layer."""
    if not name.startswith(LOD_PREFIX):
        return None
    try:
        return int(name[len(LOD_PREFIX) :])
    except ValueError:
        return None


# FGDB Layer Management

def _require_pyogrio() -> None:
    if not _HAS_PYOGRIO:
        raise RuntimeError(
            "pyogrio is required for FGDB operations. "
            "Install with: pip install pyogrio"
        )


def _write_to_gdb(
    gdf: gpd.GeoDataFrame,
    gdb_path: Path,
    layer_name: str,
) -> None:
    """
    Write a GeoDataFrame as a layer in an existing FGDB.

    If the layer already exists, it is overwritten (GDAL OpenFileGDB
    replaces the layer when writing without append=True).
    """
    _require_pyogrio()

    pyogrio.write_dataframe(
        gdf,
        str(gdb_path),
        layer=layer_name,
        driver="OpenFileGDB",
    )


def _read_from_gdb(gdb_path: Path, layer_name: str) -> gpd.GeoDataFrame:
    """Read a single layer from an FGDB via geo_io."""
    return read_vector(gdb_path, layer=layer_name)


def _gdb_layer_names(gdb_path: Path) -> Set[str]:
    """Return the set of all layer names in a FGDB."""
    return {rec["name"] for rec in list_layers(gdb_path) if rec["name"]}


# Metadata Table (FGDB tag + model info)

def _make_meta_gdf(
    created_at: Optional[str] = None,
) -> gpd.GeoDataFrame:
    """
    Build the metadata GeoDataFrame (single row).

    Uses a placeholder Point(0, 0) geometry because FGDB feature classes
    require a geometry column. The geometry is meaningless; the attributes
    are what matter.
    """
    now = created_at or datetime.datetime.now(datetime.timezone.utc).isoformat()
    return gpd.GeoDataFrame(
        {
            "magic_id": [MAGIC_ID],
            "schema_version": [SCHEMA_VERSION],
            "created_at": [now],
        },
        geometry=[Point(0, 0)],
        crs="EPSG:4326",
    )


def _make_log_entry(
    raster_path: str,
    method: str,
    resolution: int,
    lod: str,
    cell_count: int,
    action: str,
    source_surface: Optional[str] = None,
) -> Dict[str, Any]:
    """Build a single commit log row as a dict."""
    entry = {
        "committed_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "lod_name": lod,
        "resolution": resolution,
        "raster_source": str(raster_path),
        "method": method,
        "cell_count": cell_count,
        "action": action,
    }
    if source_surface is not None:
        entry["surface_source"] = source_surface
    return entry


def _read_log(gdb_path: Path) -> gpd.GeoDataFrame:
    """Read the commit log table, or return an empty GeoDataFrame if it doesn't exist."""
    existing = _gdb_layer_names(gdb_path)
    if LOG_TABLE not in existing:
        return gpd.GeoDataFrame(
            columns=[
                "committed_at",
                "lod_name",
                "resolution",
                "raster_source",
                "method",
                "cell_count",
                "action",
                "geometry",
            ],
            geometry="geometry",
            crs="EPSG:4326",
        )
    return _read_from_gdb(gdb_path, LOG_TABLE)


def _append_log_entries(
    gdb_path: Path, entries: List[Dict[str, Any]]
) -> None:
    """
    Append commit log entries to the log table.

    Reads the existing log, concatenates new entries, and overwrites
    the table (read-modify-write pattern).
    """
    import pandas as pd

    existing_log = _read_log(gdb_path)

    new_rows = gpd.GeoDataFrame(
        entries,
        geometry=[Point(0, 0)] * len(entries),
        crs="EPSG:4326",
    )

    updated_log = pd.concat([existing_log, new_rows], ignore_index=True)
    updated_log = gpd.GeoDataFrame(updated_log, geometry="geometry", crs="EPSG:4326")

    _write_to_gdb(updated_log, gdb_path, LOG_TABLE)


# Surface Model Validation

def is_surface_model(gdb_path: Path) -> bool:
    """
    Check whether a FGDB is tagged as a Hex Surface Model.

    Reads the metadata table and validates the magic identifier.
    """
    if not gdb_path.exists():
        return False

    existing = _gdb_layer_names(gdb_path)
    if META_TABLE not in existing:
        return False

    try:
        meta_gdf = _read_from_gdb(gdb_path, META_TABLE)
        if meta_gdf.empty:
            return False
        magic = str(meta_gdf.iloc[0].get("magic_id", ""))
        return magic == MAGIC_ID
    except Exception:
        return False


# Raster Footprint Extraction

def raster_footprint(raster_path: Path) -> gpd.GeoDataFrame:
    """
    Extract the bounding box of a raster as a single-polygon GeoDataFrame
    in EPSG:4326.
    """
    with rasterio.open(str(raster_path)) as src:
        bounds = src.bounds  # (left, bottom, right, top)
        crs = src.crs

    poly = box(bounds.left, bounds.bottom, bounds.right, bounds.top)
    gdf = gpd.GeoDataFrame(
        {"source": [str(raster_path)]},
        geometry=[poly],
        crs=crs,
    )

    # H3 needs EPSG:4326
    if gdf.crs is not None:
        try:
            epsg = gdf.crs.to_epsg()
        except Exception:
            epsg = None
        if epsg != 4326:
            gdf = gdf.to_crs(epsg=4326)

    return gdf


def union_raster_footprints(raster_paths: List[Path]) -> gpd.GeoDataFrame:
    """
    Union the bounding boxes of multiple rasters into a single extent polygon.

    All footprints are reprojected to EPSG:4326 before unioning.
    """
    import pandas as pd

    footprints = [raster_footprint(p) for p in raster_paths]
    combined = pd.concat(footprints, ignore_index=True)
    combined = gpd.GeoDataFrame(combined, geometry="geometry", crs="EPSG:4326")

    # Dissolve into a single polygon
    unioned = combined.dissolve()
    return unioned.reset_index(drop=True)


# Column Classification & Pruning

def _raster_data_columns(gdf: gpd.GeoDataFrame) -> List[str]:
    """
    Identify columns that contain raster-derived data (not system or metadata).
    """
    return [
        col
        for col in gdf.columns
        if col not in SYSTEM_COLUMNS and not col.startswith(META_COLUMN_PREFIX)
    ]


def _strip_meta_columns(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Remove meta_* columns from a GeoDataFrame.

    Provenance is tracked in the commit log, so per-cell meta columns
    from RasterBinEnricher are redundant in the LOD feature classes.
    """
    meta_cols = [c for c in gdf.columns if c.startswith(META_COLUMN_PREFIX)]
    if meta_cols:
        gdf = gdf.drop(columns=meta_cols)
    return gdf


def prune_empty_cells(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Drop rows where ALL raster-derived columns are null.

    Cells that landed entirely in NoData across every sampled raster
    are removed so we don't store empty hexagons.
    """
    data_cols = _raster_data_columns(gdf)
    if not data_cols:
        return gdf

    # Keep rows where at least one data column is non-null
    mask = gdf[data_cols].notna().any(axis=1)
    dropped = (~mask).sum()
    if dropped > 0:
        print(f"  Pruned {dropped} empty cell(s) (no raster data)")
    return gdf[mask].reset_index(drop=True)


# LOD Merge Logic

def merge_lod(
    existing_gdf: gpd.GeoDataFrame,
    new_gdf: gpd.GeoDataFrame,
) -> gpd.GeoDataFrame:
    """
    Merge new attributed cells into an existing LOD feature class.

    Join semantics (on GRID_ID):
      - Cells in both: new values overwrite old on shared columns.
        Where new is null, old values are preserved.
      - Cells only in existing: kept as-is; new-only columns filled with null.
      - Cells only in new: appended; existing-only columns filled with null.
    """
    crs = new_gdf.crs or existing_gdf.crs

    # Separate geometry from attributes — geometry doesn't play well with
    # combine_first, so we handle it manually.
    existing_geom = existing_gdf.set_index("GRID_ID")["geometry"]
    new_geom = new_gdf.set_index("GRID_ID")["geometry"]

    existing_attrs = existing_gdf.set_index("GRID_ID").drop(columns=["geometry"])
    new_attrs = new_gdf.set_index("GRID_ID").drop(columns=["geometry"])

    # combine_first: caller's non-null values take priority, then fills gaps
    # from the argument. So new_attrs.combine_first(existing_attrs) means
    # "use new where available, fill gaps from existing".
    merged_attrs = new_attrs.combine_first(existing_attrs)

    # Geometry merge: new wins, old fills cells not in new data.
    geom_lookup: Dict[str, Any] = {}
    for gid, geom in existing_geom.items():
        geom_lookup[gid] = geom
    for gid, geom in new_geom.items():
        geom_lookup[gid] = geom  # New overwrites old

    merged_geom = [geom_lookup[gid] for gid in merged_attrs.index]

    # Reconstruct GeoDataFrame
    result = gpd.GeoDataFrame(
        merged_attrs.reset_index(),
        geometry=merged_geom,
        crs=crs,
    )

    # Ensure 'res' is integer (combine_first can promote to float)
    if "res" in result.columns:
        result["res"] = result["res"].astype(int)

    return result


# Raster Config Parsing

def parse_raster_config(config_str: str, surface_mode: bool = False) -> RasterConfig:
    """
    Parse a raster config string into a RasterConfig.

    Formats:
      - Standard:      path:method:resolution  (e.g., dem.tif:mean:9)
      - Surface mode:  path:method             (e.g., dem.tif:mean)

    In surface mode, resolution is read from the surface's 'res' column
    instead of being specified per-raster.
    """
    if surface_mode:
        # path:method (no resolution)
        parts = config_str.rsplit(":", 1)
        if len(parts) == 3:
            # User provided resolution anyway — allow it as override
            return _parse_standard_config(config_str)
        if len(parts) != 2:
            raise ValueError(
                f"Invalid raster config: '{config_str}'. "
                "With --surface, expected format: path:method (e.g., dem.tif:mean)"
            )
        raw_path, method = parts
        path = Path(raw_path).expanduser().resolve()
        if not path.exists():
            raise FileNotFoundError(f"Raster not found: {path}")
        method = method.strip().lower()
        if method not in ("mean", "majority"):
            raise ValueError(
                f"Invalid method '{method}' for {path.name}. "
                "Must be 'mean' or 'majority'."
            )
        return RasterConfig(path=path, method=method, resolution=None)
    else:
        return _parse_standard_config(config_str)


def _parse_standard_config(config_str: str) -> RasterConfig:
    """Parse a 'path:method:resolution' string into a RasterConfig."""
    parts = config_str.rsplit(":", 2)
    if len(parts) != 3:
        raise ValueError(
            f"Invalid raster config: '{config_str}'. "
            "Expected format: path:method:resolution (e.g., dem.tif:mean:9)"
        )

    raw_path, method, raw_res = parts

    path = Path(raw_path).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(f"Raster not found: {path}")

    method = method.strip().lower()
    if method not in ("mean", "majority"):
        raise ValueError(
            f"Invalid method '{method}' for {path.name}. "
            "Must be 'mean' or 'majority'."
        )

    try:
        res = int(raw_res.strip())
    except ValueError:
        raise ValueError(
            f"Invalid resolution '{raw_res}' for {path.name}. "
            "Must be an integer 0..15."
        )
    if not (0 <= res <= 15):
        raise ValueError(
            f"Resolution {res} out of range for {path.name}. Must be 0..15."
        )

    return RasterConfig(path=path, method=method, resolution=res)


def group_by_resolution(
    configs: List[RasterConfig],
) -> Dict[int, List[RasterConfig]]:
    """Group raster configs by their target H3 resolution."""
    groups: Dict[int, List[RasterConfig]] = {}
    for cfg in configs:
        if cfg.resolution is not None:
            groups.setdefault(cfg.resolution, []).append(cfg)
    return groups


# Pre-built Surface Loading

def load_surface(
    surface_path: Path,
    layer: Optional[str] = None,
) -> Dict[int, gpd.GeoDataFrame]:
    """
    Load a pre-built H3 surface and group by resolution.

    The surface must have 'GRID_ID' and 'res' columns. Returns a dict
    mapping resolution -> GeoDataFrame subset.

    Parameters
    ----------
    surface_path : Path
        Path to vector file containing the H3 surface.
    layer : str, optional
        Layer name for multi-layer sources.

    Returns
    -------
    dict[int, GeoDataFrame]
        Resolution -> subset of the surface at that resolution.
    """
    gdf = read_vector(surface_path, layer=layer, epsg=4326)

    if "GRID_ID" not in gdf.columns:
        raise ValueError(
            f"Surface file must contain a 'GRID_ID' column. "
            f"Found columns: {list(gdf.columns)}"
        )
    if "res" not in gdf.columns:
        raise ValueError(
            f"Surface file must contain a 'res' column. "
            f"Found columns: {list(gdf.columns)}"
        )

    # Group by resolution
    groups: Dict[int, gpd.GeoDataFrame] = {}
    for res_val, subset in gdf.groupby("res"):
        groups[int(res_val)] = subset.reset_index(drop=True)

    return groups


# LOD Processing (parallelizable unit of work)

def _process_lod_rasters(
    surface: gpd.GeoDataFrame,
    raster_cfgs: List[RasterConfig],
) -> gpd.GeoDataFrame:
    """
    Sample all rasters for a single LOD into the provided surface.

    This is the hot path — enriches one surface with N rasters.
    Returns the enriched, pruned GeoDataFrame (no meta columns).
    """
    enricher = RasterBinEnricher(surface)
    for cfg in raster_cfgs:
        enricher.process_raster(
            raster_path=str(cfg.path),
            method=cfg.method,
        )

    result = enricher.result
    result = _strip_meta_columns(result)
    result = prune_empty_cells(result)
    return result


# Subcommand: init

def cmd_init(args: argparse.Namespace) -> int:
    """Create a new Hex Surface Model FGDB."""
    _require_pyogrio()

    target = Path(args.target).expanduser().resolve()

    if target.suffix.lower() != ".gdb":
        print(f"[ERROR] Target must end with .gdb: {target}", file=sys.stderr)
        return 2

    if target.exists():
        if not args.overwrite:
            print(
                f"[ERROR] {target} already exists.\n"
                "       Use --overwrite to delete and recreate it.",
                file=sys.stderr,
            )
            return 2
        import shutil

        shutil.rmtree(target)
        print(f"  Removed existing: {target}")

    # Create the FGDB by writing the metadata table
    target.parent.mkdir(parents=True, exist_ok=True)
    meta_gdf = _make_meta_gdf()

    pyogrio.write_dataframe(
        meta_gdf,
        str(target),
        layer=META_TABLE,
        driver="OpenFileGDB",
    )

    print(f"[OK] Initialized Hex Surface Model: {target}")
    print(f"     Schema version: {SCHEMA_VERSION}")
    return 0


# Subcommand: commit

def cmd_commit(args: argparse.Namespace) -> int:
    """Tessellate raster footprints, sample values, write to Surface Model."""
    _require_pyogrio()

    target = Path(args.target).expanduser().resolve()

    # Validate target
    if not target.exists():
        print(
            f"[ERROR] FGDB not found: {target}\n"
            f"       Run 'build_surface.py init --target {target}' first.",
            file=sys.stderr,
        )
        return 2

    if not is_surface_model(target):
        print(
            f"[ERROR] {target.name} is not a Hex Surface Model.\n"
            "       Only FGDBs created with 'build_surface.py init' can receive commits.",
            file=sys.stderr,
        )
        return 2

    # Parse raster configs
    if not args.raster:
        print("[ERROR] Provide at least one --raster path:method:resolution.", file=sys.stderr)
        return 2

    surface_mode = args.surface is not None

    try:
        configs = [parse_raster_config(r, surface_mode=surface_mode) for r in args.raster]
    except (ValueError, FileNotFoundError) as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        return 2

    existing_layers = _gdb_layer_names(target)
    log_entries: List[Dict[str, Any]] = []

    t_commit_start = time.perf_counter()
    surface_path_str = None

    if surface_mode:
        # ── Pre-built surface mode ───────────────────────────────────
        surface_path = Path(args.surface).expanduser().resolve()
        surface_path_str = str(surface_path)

        if not surface_path.exists():
            print(f"[ERROR] Surface file not found: {surface_path}", file=sys.stderr)
            return 2

        print(f"Loading pre-built surface: {surface_path}")
        t_load_surface = time.perf_counter()
        surface_groups = load_surface(surface_path, layer=args.surface_layer)
        print(f"[TIME] Load surface: {time.perf_counter() - t_load_surface:.3f}s")

        resolutions_in_surface = sorted(surface_groups.keys())
        print(f"Surface contains {sum(len(g) for g in surface_groups.values())} cells "
              f"across {len(resolutions_in_surface)} resolution(s): {resolutions_in_surface}")

        # Assign resolution to raster configs that don't have one
        if len(resolutions_in_surface) == 1:
            res = resolutions_in_surface[0]
            resolved_configs = []
            for cfg in configs:
                if cfg.resolution is None:
                    resolved_configs.append(
                        RasterConfig(path=cfg.path, method=cfg.method, resolution=res)
                    )
                else:
                    resolved_configs.append(cfg)
            configs = resolved_configs
        else:
            # Multiple resolutions in surface — rasters without explicit
            # resolution are applied to ALL resolutions
            resolved_configs = []
            for cfg in configs:
                if cfg.resolution is None:
                    for res in resolutions_in_surface:
                        resolved_configs.append(
                            RasterConfig(path=cfg.path, method=cfg.method, resolution=res)
                        )
                else:
                    resolved_configs.append(cfg)
            configs = resolved_configs

        resolution_groups = group_by_resolution(configs)

        print(
            f"\nTarget: {target}"
            f"\nCommitting {len(configs)} raster(s) across "
            f"{len(resolution_groups)} LOD(s) (surface mode)\n"
        )

        # Process each resolution group using the pre-built surface
        for res, raster_cfgs in sorted(resolution_groups.items()):
            lod = lod_name(res)
            raster_names = [cfg.path.name for cfg in raster_cfgs]

            t_lod = time.perf_counter()
            print(f"── {lod} (resolution {res}) ─────────────────────────────")
            print(f"  Rasters: {', '.join(raster_names)}")

            if res not in surface_groups:
                print(f"  [WARN] Resolution {res} not in surface — skipping")
                continue

            surface = surface_groups[res]
            print(f"  Using {len(surface)} pre-built H3 cells")

            # Sample rasters into the surface
            t_sample = time.perf_counter()
            result = _process_lod_rasters(surface, raster_cfgs)
            print(f"  [TIME] Raster sampling: {time.perf_counter() - t_sample:.3f}s")

            if result.empty:
                print(f"  [WARN] No cells with data for {lod} — skipping")
                continue

            # Merge or create the LOD feature class
            t_merge = time.perf_counter()
            if lod in existing_layers:
                print(f"  Merging into existing {lod}...")
                existing_gdf = _read_from_gdb(target, lod)
                merged = merge_lod(existing_gdf, result)
                action = "merge"
            else:
                print(f"  Creating new feature class: {lod}")
                merged = result
                action = "create"
            print(f"  [TIME] Merge/create: {time.perf_counter() - t_merge:.3f}s")

            t_write = time.perf_counter()
            _write_to_gdb(merged, target, lod)
            existing_layers.add(lod)
            print(f"  [TIME] FGDB write: {time.perf_counter() - t_write:.3f}s")
            print(f"  [OK] {lod}: {len(merged)} cells, "
                  f"{len(_raster_data_columns(merged))} data fields")
            print(f"  [TIME] LOD total: {time.perf_counter() - t_lod:.3f}s")

            for cfg in raster_cfgs:
                log_entries.append(
                    _make_log_entry(
                        raster_path=str(cfg.path),
                        method=cfg.method,
                        resolution=res,
                        lod=lod,
                        cell_count=len(merged),
                        action=action,
                        source_surface=surface_path_str,
                    )
                )

    else:
        # ── Standard mode (tessellate from raster footprints) ────────
        resolution_groups = group_by_resolution(configs)

        print(f"Target: {target}")
        print(
            f"Committing {len(configs)} raster(s) across "
            f"{len(resolution_groups)} LOD(s)\n"
        )

        for res, raster_cfgs in sorted(resolution_groups.items()):
            lod = lod_name(res)
            raster_paths = [cfg.path for cfg in raster_cfgs]
            raster_names = [cfg.path.name for cfg in raster_cfgs]

            t_lod = time.perf_counter()
            print(f"── {lod} (resolution {res}) ─────────────────────────────")
            print(f"  Rasters: {', '.join(raster_names)}")

            # 1. Extract union of raster footprints
            t_footprint = time.perf_counter()
            print("  Extracting raster footprints...")
            extent = union_raster_footprints(raster_paths)
            print(f"  [TIME] Footprint extraction: {time.perf_counter() - t_footprint:.3f}s")

            # 2. Build H3 surface from the extent
            print(f"  Tessellating at H3 resolution {res}...")
            surface = build_h3_surface(extent, res)
            print(f"  Generated {len(surface)} H3 cells")

            # 3. Sample rasters into the surface
            t_sample = time.perf_counter()
            result = _process_lod_rasters(surface, raster_cfgs)
            print(f"  [TIME] Raster sampling: {time.perf_counter() - t_sample:.3f}s")

            if result.empty:
                print(f"  [WARN] No cells with data for {lod} — skipping")
                continue

            # 4. Merge or create the LOD feature class
            t_merge = time.perf_counter()
            if lod in existing_layers:
                print(f"  Merging into existing {lod}...")
                existing_gdf = _read_from_gdb(target, lod)
                merged = merge_lod(existing_gdf, result)
                action = "merge"
            else:
                print(f"  Creating new feature class: {lod}")
                merged = result
                action = "create"
            print(f"  [TIME] Merge/create: {time.perf_counter() - t_merge:.3f}s")

            # 5. Write to FGDB
            t_write = time.perf_counter()
            _write_to_gdb(merged, target, lod)
            existing_layers.add(lod)
            print(f"  [TIME] FGDB write: {time.perf_counter() - t_write:.3f}s")
            print(f"  [OK] {lod}: {len(merged)} cells, "
                  f"{len(_raster_data_columns(merged))} data fields")
            print(f"  [TIME] LOD total: {time.perf_counter() - t_lod:.3f}s")

            # 6. Record commit log entries
            for cfg in raster_cfgs:
                log_entries.append(
                    _make_log_entry(
                        raster_path=str(cfg.path),
                        method=cfg.method,
                        resolution=res,
                        lod=lod,
                        cell_count=len(merged),
                        action=action,
                    )
                )

    # Write commit log
    if log_entries:
        t_log = time.perf_counter()
        _append_log_entries(target, log_entries)
        print(f"\n[TIME] Commit log write: {time.perf_counter() - t_log:.3f}s")
        print(f"[OK] Committed {len(log_entries)} raster(s) to {target.name}")
        print(f"[TIME] Total commit: {time.perf_counter() - t_commit_start:.3f}s")
    else:
        print("\n[WARN] Nothing was committed.")
        return 1

    return 0


# Subcommand: info

def cmd_info(args: argparse.Namespace) -> int:
    """Display Hex Surface Model metadata and commit history."""
    target = Path(args.target).expanduser().resolve()

    if not target.exists():
        print(f"[ERROR] Not found: {target}", file=sys.stderr)
        return 2

    if not is_surface_model(target):
        print(f"[ERROR] {target.name} is not a Hex Surface Model.", file=sys.stderr)
        return 2

    # Model header
    meta_gdf = _read_from_gdb(target, META_TABLE)
    meta = meta_gdf.iloc[0]
    created = meta.get("created_at", "unknown")
    version = meta.get("schema_version", "unknown")

    print(f"\nHex Surface Model: {target}")
    print(f"Created:           {created}")
    print(f"Schema version:    {version}")

    # LOD feature classes
    all_layers = list_layers(target)
    lod_layers = []
    for rec in all_layers:
        name = rec.get("name", "")
        res = parse_lod_name(name)
        if res is not None:
            lod_layers.append((name, res))

    lod_layers.sort(key=lambda x: x[1])

    print(f"\nFeature Classes: ({len(lod_layers)} LODs)")
    if not lod_layers:
        print("  (none)")
    else:
        print(f"  {'Name':<12} {'Res':>4} {'Cells':>10}  Data Fields")
        print(f"  {'─' * 12} {'─' * 4} {'─' * 10}  {'─' * 40}")
        for name, res in lod_layers:
            try:
                gdf = _read_from_gdb(target, name)
                n_cells = len(gdf)
                data_cols = _raster_data_columns(gdf)
                fields_str = ", ".join(data_cols) if data_cols else "(none)"
            except Exception:
                n_cells = "?"
                fields_str = "(read error)"
            print(f"  {name:<12} {res:>4} {n_cells:>10}  {fields_str}")

    # Commit history
    log_gdf = _read_log(target)
    print(f"\nCommit History: ({len(log_gdf)} entries)")
    if log_gdf.empty:
        print("  (none)")
    else:
        print(f"  {'Timestamp':<28} {'LOD':<8} {'Method':<10} {'Cells':>7}  {'Action':<7}  Source")
        print(f"  {'─' * 28} {'─' * 8} {'─' * 10} {'─' * 7}  {'─' * 7}  {'─' * 30}")
        for _, row in log_gdf.iterrows():
            ts = str(row.get("committed_at", ""))[:26]
            lod = str(row.get("lod_name", ""))
            method = str(row.get("method", ""))
            cells = row.get("cell_count", "?")
            action = str(row.get("action", ""))
            source = Path(str(row.get("raster_source", ""))).name
            print(f"  {ts:<28} {lod:<8} {method:<10} {cells:>7}  {action:<7}  {source}")

    print()
    return 0


# CLI

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="build_surface",
        description="Manage Hex Surface Model FGDBs.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # init
    p_init = sub.add_parser(
        "init",
        help="Create a new Hex Surface Model FGDB.",
    )
    p_init.add_argument(
        "--target",
        required=True,
        help="Path for the new .gdb (must end in .gdb).",
    )
    p_init.add_argument(
        "--overwrite",
        action="store_true",
        help="Delete and recreate if the target already exists.",
    )

    # commit
    p_commit = sub.add_parser(
        "commit",
        help="Commit raster data into an existing Surface Model.",
    )
    p_commit.add_argument(
        "--target",
        required=True,
        help="Path to the Surface Model .gdb.",
    )
    p_commit.add_argument(
        "--raster",
        action="append",
        required=True,
        help=(
            "Raster config as path:method:resolution (standard) or "
            "path:method (with --surface). Repeatable. "
            "Example: elevation.tif:mean:9 or elevation.tif:mean"
        ),
    )
    p_commit.add_argument(
        "--surface",
        default=None,
        help=(
            "Path to a pre-built H3 surface vector file (.geojson, .shp, .gdb). "
            "Must contain GRID_ID and res columns. When provided, tessellation "
            "is skipped and rasters are sampled directly into the provided cells. "
            "Raster configs may omit resolution (inferred from surface)."
        ),
    )
    p_commit.add_argument(
        "--surface-layer",
        default=None,
        help="Layer name within the surface file (for FGDB/KML sources).",
    )

    # info
    p_info = sub.add_parser(
        "info",
        help="Display Surface Model metadata and commit history.",
    )
    p_info.add_argument(
        "--target",
        required=True,
        help="Path to the Surface Model .gdb.",
    )

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    dispatch = {
        "init": cmd_init,
        "commit": cmd_commit,
        "info": cmd_info,
    }

    handler = dispatch.get(args.command)
    if handler is None:
        parser.print_help()
        return 2

    try:
        return handler(args)
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())