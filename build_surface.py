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
         --raster path:method:resolution[:factor]

      2. Pre-built surface: skip tessellation, sample into provided cells
         --surface path.geojson --raster path:method[:factor]

info
    Display model metadata, LOD summary, and commit history.

Factor Schema
-------------
    Factor 1 (F1) — Elevation
    Factor 2 (F2) — Slope
    Factor 3 (F3) — Vegetation
    Factor 4 (F4) — Soil

    When a factor tag is provided, the sampled column is renamed from the
    raster-derived name to F1, F2, F3, or F4. Raw values are stored — no
    normalization occurs at this stage.

    The commit log records which raster contributed to which factor for
    full provenance.

Usage:
    # Create a new surface model
    python build_surface.py init --target "path/to/model.gdb"

    # Commit rasters with factor tags
    python build_surface.py commit --target "path/to/model.gdb" \\
        --raster elevation.tif:mean:9:1 \\
        --raster slope.tif:mean:9:2 \\
        --raster landcover.tif:majority:9:3 \\
        --raster soil.tif:majority:9:4

    # Commit without factor tags (legacy; column named after raster)
    python build_surface.py commit --target "path/to/model.gdb" \\
        --raster elevation.tif:mean:9

    # Commit using a pre-built H3 surface
    python build_surface.py commit --target "path/to/model.gdb" \\
        --surface path/to/h3_cells.geojson \\
        --raster elevation.tif:mean:1

    # View model info
    python build_surface.py info --target "path/to/model.gdb"

Requires:
    geo_io.py, generate_h3_tesselation.py, sample_raster.py
    geopandas, pyogrio, rasterio, shapely, h3, numpy
"""


import argparse
import datetime
import os
import shutil
import sys
import time
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


# ═══════════════════════════════════════════════════════════════════════════
# Constants
# ═══════════════════════════════════════════════════════════════════════════

META_TABLE = "__hex_surface_meta__"
LOG_TABLE = "__hex_surface_log__"
MAGIC_ID = "HEX_SURFACE_MODEL"
SCHEMA_VERSION = 2  # Bumped for factor tag support
LOD_PREFIX = "LOD_"

# Factor definitions: number -> human-readable name
FACTOR_NAMES: Dict[int, str] = {
    1: "Elevation",
    2: "Slope",
    3: "Vegetation",
    4: "Soil",
}

# Valid factor numbers
VALID_FACTORS = frozenset(FACTOR_NAMES.keys())

# Columns that belong to the H3 grid structure, not raster data
SYSTEM_COLUMNS = frozenset({"GRID_ID", "res", "geometry"})

# Columns injected by RasterBinEnricher that we strip before writing to the
# LOD feature class — provenance is tracked in the commit log instead.
META_COLUMN_PREFIX = "meta_"


# Data Structures

@dataclass(frozen=True)
class RasterConfig:
    """A single raster input with its sampling method, target resolution, and factor."""

    path: Path
    method: str                  # "mean" or "majority"
    resolution: Optional[int]    # H3 resolution (0..15), None when using --surface
    factor: Optional[int]        # 1-4, or None for untagged (legacy)


# LOD Naming

def lod_name(res: int) -> str:
    """H3 resolution -> LOD feature class name (e.g., 9 -> 'LOD_09')."""
    return f"{LOD_PREFIX}{res:02d}"


def parse_lod_name(name: str) -> Optional[int]:
    """LOD feature class name -> H3 resolution, or None if not a LOD layer."""
    if not name.startswith(LOD_PREFIX):
        return None
    try:
        return int(name[len(LOD_PREFIX):])
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
    """Write a GeoDataFrame as a layer in an existing FGDB."""
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
    require a geometry column.
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
    factor: Optional[int] = None,
    source_surface: Optional[str] = None,
) -> Dict[str, Any]:
    """Build a single commit log row as a dict."""
    return {
        "committed_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "lod": lod,
        "resolution": resolution,
        "method": method,
        "cell_count": cell_count,
        "action": action,
        "raster_source": raster_path,
        "factor": factor if factor is not None else -1,
        "factor_name": FACTOR_NAMES.get(factor, "") if factor else "",
        "source_surface": source_surface or "",
    }


def _append_log_entries(
    gdb_path: Path,
    entries: List[Dict[str, Any]],
) -> None:
    """Append commit log entries to the log table in the FGDB."""
    if not entries:
        return

    new_log = gpd.GeoDataFrame(
        entries,
        geometry=[Point(0, 0)] * len(entries),
        crs="EPSG:4326",
    )

    existing = _gdb_layer_names(gdb_path)
    if LOG_TABLE in existing:
        old_log = _read_from_gdb(gdb_path, LOG_TABLE)
        import pandas as pd
        combined = gpd.GeoDataFrame(
            pd.concat([old_log, new_log], ignore_index=True),
            crs="EPSG:4326",
        )
    else:
        combined = new_log

    _write_to_gdb(combined, gdb_path, LOG_TABLE)


def is_surface_model(gdb_path: Path) -> bool:
    """Check whether a FGDB is tagged as a Hex Surface Model."""
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
    """Extract the bounding box of a raster as a GeoDataFrame in EPSG:4326."""
    with rasterio.open(str(raster_path)) as src:
        bounds = src.bounds
        crs = src.crs

    poly = box(bounds.left, bounds.bottom, bounds.right, bounds.top)
    gdf = gpd.GeoDataFrame(
        {"source": [str(raster_path)]},
        geometry=[poly],
        crs=crs,
    )

    if gdf.crs is not None:
        try:
            epsg = gdf.crs.to_epsg()
        except Exception:
            epsg = None
        if epsg != 4326:
            gdf = gdf.to_crs(epsg=4326)

    return gdf


def union_raster_footprints(raster_paths: List[Path]) -> gpd.GeoDataFrame:
    """Union the bounding boxes of multiple rasters into a single extent polygon."""
    import pandas as pd

    footprints = [raster_footprint(p) for p in raster_paths]
    combined = pd.concat(footprints, ignore_index=True)
    combined = gpd.GeoDataFrame(combined, geometry="geometry", crs="EPSG:4326")
    unioned = combined.dissolve()
    return unioned.reset_index(drop=True)


# Column Classification & Pruning

# Factor columns are F1, F2, F3, F4
FACTOR_COLUMNS = frozenset(f"F{i}" for i in VALID_FACTORS)


def _raster_data_columns(gdf: gpd.GeoDataFrame) -> List[str]:
    """Identify columns that contain raster-derived data."""
    return [
        col
        for col in gdf.columns
        if col not in SYSTEM_COLUMNS and not col.startswith(META_COLUMN_PREFIX)
    ]


def _strip_meta_columns(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Remove meta_* columns from a GeoDataFrame."""
    meta_cols = [c for c in gdf.columns if c.startswith(META_COLUMN_PREFIX)]
    if meta_cols:
        gdf = gdf.drop(columns=meta_cols)
    return gdf


def prune_empty_cells(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Drop rows where ALL raster-derived columns are null."""
    data_cols = _raster_data_columns(gdf)
    if not data_cols:
        return gdf

    mask = gdf[data_cols].notna().any(axis=1)
    dropped = (~mask).sum()
    if dropped > 0:
        print(f"  Pruned {dropped} empty cell(s) (no raster data)")
    return gdf[mask].reset_index(drop=True)


def _rename_sampled_to_factor(
    gdf: gpd.GeoDataFrame,
    raster_cfg: RasterConfig,
) -> gpd.GeoDataFrame:
    """
    Rename the raster-sampled column to its factor name (F1, F2, F3, F4).

    After RasterBinEnricher runs, columns are named like:
        {raster_stem}_{method}  (e.g., "elevation_mean", "landcover_majority")

    If a factor tag was provided, we rename that column to F{n}.
    If not, we leave it as-is (legacy behavior).
    """
    if raster_cfg.factor is None:
        return gdf

    # Determine what RasterBinEnricher named the column
    raster_stem = raster_cfg.path.stem
    sampled_col = f"{raster_stem}_{raster_cfg.method}"

    target_col = f"F{raster_cfg.factor}"

    if sampled_col in gdf.columns:
        # If target column already exists (from a previous raster in same batch),
        # the new one overwrites it
        if target_col in gdf.columns and target_col != sampled_col:
            gdf = gdf.drop(columns=[target_col])
        gdf = gdf.rename(columns={sampled_col: target_col})
    else:
        # Fallback: look for any new column that wasn't in SYSTEM_COLUMNS
        # (shouldn't normally reach here, but defensive)
        print(f"  [WARN] Expected column '{sampled_col}' not found after sampling")

    return gdf


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

    existing_geom = existing_gdf.set_index("GRID_ID")["geometry"]
    new_geom = new_gdf.set_index("GRID_ID")["geometry"]

    existing_attrs = existing_gdf.set_index("GRID_ID").drop(columns=["geometry"])
    new_attrs = new_gdf.set_index("GRID_ID").drop(columns=["geometry"])

    merged_attrs = new_attrs.combine_first(existing_attrs)

    geom_lookup: Dict[str, Any] = {}
    for gid, geom in existing_geom.items():
        geom_lookup[gid] = geom
    for gid, geom in new_geom.items():
        geom_lookup[gid] = geom

    merged_geom = [geom_lookup[gid] for gid in merged_attrs.index]

    result = gpd.GeoDataFrame(
        merged_attrs.reset_index(),
        geometry=merged_geom,
        crs=crs,
    )

    if "res" in result.columns:
        result["res"] = result["res"].astype(int)

    return result


# Raster Config Parsing

def parse_raster_config(config_str: str, surface_mode: bool = False) -> RasterConfig:
    """
    Parse a raster config string into a RasterConfig.

    Formats:
      Standard:      path:method:resolution[:factor]  (e.g., dem.tif:mean:9:1)
      Surface mode:  path:method[:factor]              (e.g., dem.tif:mean:1)

    Factor is optional. When provided (1-4), the sampled column is renamed
    to F1, F2, F3, or F4 in the output feature class.
    """
    if surface_mode:
        return _parse_surface_config(config_str)
    else:
        return _parse_standard_config(config_str)


def _parse_standard_config(config_str: str) -> RasterConfig:
    """Parse 'path:method:resolution[:factor]' into a RasterConfig."""
    parts = config_str.rsplit(":", 3)

    if len(parts) == 4:
        raw_path, method, raw_res, raw_factor = parts
        factor = _validate_factor(raw_factor, raw_path)
    elif len(parts) == 3:
        raw_path, method, raw_res = parts
        factor = None
    else:
        raise ValueError(
            f"Invalid raster config: '{config_str}'. "
            "Expected format: path:method:resolution[:factor] "
            "(e.g., dem.tif:mean:9 or dem.tif:mean:9:1)"
        )

    path = _validate_path(raw_path)
    method = _validate_method(method, path.name)
    res = _validate_resolution(raw_res, path.name)

    return RasterConfig(path=path, method=method, resolution=res, factor=factor)


def _parse_surface_config(config_str: str) -> RasterConfig:
    """Parse 'path:method[:factor]' (surface mode) into a RasterConfig."""
    parts = config_str.rsplit(":", 3)

    # Could be path:method:factor, path:method:resolution:factor,
    # path:method:resolution, or path:method
    if len(parts) == 4:
        # path:method:resolution:factor (user provided resolution + factor)
        raw_path, method, raw_res, raw_factor = parts
        path = _validate_path(raw_path)
        method = _validate_method(method, path.name)
        res = _validate_resolution(raw_res, path.name)
        factor = _validate_factor(raw_factor, raw_path)
        return RasterConfig(path=path, method=method, resolution=res, factor=factor)
    elif len(parts) == 3:
        raw_path, method, raw_third = parts
        path = _validate_path(raw_path)
        method = _validate_method(method, path.name)

        # Disambiguate: is the third part a resolution or a factor?
        # Factors are 1-4; resolutions are 0-15.
        # If it's 1-4, prefer interpreting as factor in surface mode.
        try:
            val = int(raw_third.strip())
        except ValueError:
            raise ValueError(
                f"Invalid config: '{config_str}'. Third component "
                f"'{raw_third}' is not a valid integer."
            )

        if 1 <= val <= 4:
            # Treat as factor in surface mode
            return RasterConfig(path=path, method=method, resolution=None, factor=val)
        elif 0 <= val <= 15:
            # Treat as resolution
            return RasterConfig(path=path, method=method, resolution=val, factor=None)
        else:
            raise ValueError(
                f"Invalid config: '{config_str}'. "
                f"Third component '{val}' is not a valid resolution (0-15) or factor (1-4)."
            )
    elif len(parts) == 2:
        raw_path, method = parts
        path = _validate_path(raw_path)
        method = _validate_method(method, path.name)
        return RasterConfig(path=path, method=method, resolution=None, factor=None)
    else:
        raise ValueError(
            f"Invalid raster config: '{config_str}'. "
            "With --surface, expected: path:method[:factor] (e.g., dem.tif:mean:1)"
        )


def _validate_path(raw_path: str) -> Path:
    """Validate and resolve a raster file path."""
    path = Path(raw_path).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(f"Raster not found: {path}")
    return path


def _validate_method(method: str, filename: str) -> str:
    """Validate the sampling method."""
    method = method.strip().lower()
    if method not in ("mean", "majority"):
        raise ValueError(
            f"Invalid method '{method}' for {filename}. "
            "Must be 'mean' or 'majority'."
        )
    return method


def _validate_resolution(raw_res: str, filename: str) -> int:
    """Validate the H3 resolution."""
    try:
        res = int(raw_res.strip())
    except ValueError:
        raise ValueError(
            f"Invalid resolution '{raw_res}' for {filename}. "
            "Must be an integer 0..15."
        )
    if not (0 <= res <= 15):
        raise ValueError(
            f"Resolution {res} out of range for {filename}. Must be 0..15."
        )
    return res


def _validate_factor(raw_factor: str, filename: str) -> int:
    """Validate a factor tag."""
    try:
        factor = int(raw_factor.strip())
    except ValueError:
        raise ValueError(
            f"Invalid factor '{raw_factor}' for {filename}. "
            "Must be an integer 1-4."
        )
    if factor not in VALID_FACTORS:
        raise ValueError(
            f"Factor {factor} out of range for {filename}. "
            f"Must be 1-4: {', '.join(f'{k}={v}' for k, v in FACTOR_NAMES.items())}"
        )
    return factor


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
    """Load a pre-built H3 surface and group by resolution."""
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

    If factor tags are present, renames sampled columns to F1-F4.
    Returns the enriched, pruned GeoDataFrame (no meta columns).
    """
    enricher = RasterBinEnricher(surface)

    for cfg in raster_cfgs:
        print(f"  Sampling: {cfg.path.name} ({cfg.method})"
              + (f" -> F{cfg.factor} ({FACTOR_NAMES[cfg.factor]})" if cfg.factor else ""))
        enricher.process_raster(
            raster_path=str(cfg.path),
            method=cfg.method,
        )

    result = enricher.result
    result = _strip_meta_columns(result)

    # Rename sampled columns to factor names where tagged
    for cfg in raster_cfgs:
        result = _rename_sampled_to_factor(result, cfg)

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
        shutil.rmtree(target)
        print(f"  Removed existing: {target}")

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

    if not args.raster:
        print("[ERROR] Provide at least one --raster.", file=sys.stderr)
        return 2

    surface_mode = args.surface is not None

    try:
        configs = [parse_raster_config(r, surface_mode=surface_mode) for r in args.raster]
    except (ValueError, FileNotFoundError) as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        return 2

    # Validate: no duplicate factors within the same LOD
    _check_factor_conflicts(configs)

    existing_layers = _gdb_layer_names(target)
    log_entries: List[Dict[str, Any]] = []
    surface_path_str = None

    if surface_mode:
        return _commit_surface_mode(
            args, target, configs, existing_layers, log_entries
        )
    else:
        return _commit_standard_mode(
            target, configs, existing_layers, log_entries
        )


def _check_factor_conflicts(configs: List[RasterConfig]) -> None:
    """Warn if multiple rasters target the same factor at the same resolution."""
    seen: Dict[Tuple[Optional[int], Optional[int]], str] = {}
    for cfg in configs:
        if cfg.factor is not None:
            key = (cfg.resolution, cfg.factor)
            if key in seen:
                print(
                    f"  [WARN] Multiple rasters target F{cfg.factor} at resolution "
                    f"{cfg.resolution}: {seen[key]} and {cfg.path.name}. "
                    f"Last-write-wins: {cfg.path.name} will overwrite."
                )
            seen[key] = cfg.path.name


def _commit_standard_mode(
    target: Path,
    configs: List[RasterConfig],
    existing_layers: Set[str],
    log_entries: List[Dict[str, Any]],
) -> int:
    """Standard commit: tessellate from raster footprints."""
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

        print(f"── {lod} (resolution {res}) ─────────────────────────────")
        print(f"  Rasters: {', '.join(raster_names)}")

        t_lod = time.perf_counter()

        # 1. Extract union of raster footprints
        print("  Extracting raster footprints...")
        extent = union_raster_footprints(raster_paths)

        # 2. Build H3 surface from the extent
        print(f"  Tessellating at H3 resolution {res}...")
        surface = build_h3_surface(extent, res)
        print(f"  Generated {len(surface)} H3 cells")

        # 3. Sample rasters + rename to factor columns
        result = _process_lod_rasters(surface, raster_cfgs)

        if result.empty:
            print(f"  [WARN] No cells with data for {lod} — skipping")
            continue

        # 4. Merge or create the LOD feature class
        if lod in existing_layers:
            print(f"  Merging into existing {lod}...")
            existing_gdf = _read_from_gdb(target, lod)
            merged = merge_lod(existing_gdf, result)
            action = "merge"
        else:
            print(f"  Creating new feature class: {lod}")
            merged = result
            action = "create"

        # 5. Write to FGDB
        _write_to_gdb(merged, target, lod)
        existing_layers.add(lod)

        elapsed = time.perf_counter() - t_lod
        print(f"  [OK] {lod}: {len(merged)} cells, "
              f"{len(_raster_data_columns(merged))} data fields ({elapsed:.2f}s)")

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
                    factor=cfg.factor,
                )
            )

    _append_log_entries(target, log_entries)
    print(f"\n[OK] Committed {len(log_entries)} raster(s) to {target.name}")
    return 0


def _commit_surface_mode(
    args: argparse.Namespace,
    target: Path,
    configs: List[RasterConfig],
    existing_layers: Set[str],
    log_entries: List[Dict[str, Any]],
) -> int:
    """Surface mode: skip tessellation, sample into pre-built cells."""
    surface_path = Path(args.surface).expanduser().resolve()

    if not surface_path.exists():
        print(f"[ERROR] Surface file not found: {surface_path}", file=sys.stderr)
        return 2

    print(f"Loading pre-built surface: {surface_path}")
    surface_groups = load_surface(surface_path, layer=args.surface_layer)

    resolutions_in_surface = sorted(surface_groups.keys())
    total_cells = sum(len(g) for g in surface_groups.values())
    print(f"Surface contains {total_cells} cells across "
          f"{len(resolutions_in_surface)} resolution(s): {resolutions_in_surface}")

    # Assign resolution to configs that don't have one
    if len(resolutions_in_surface) == 1:
        res = resolutions_in_surface[0]
        configs = [
            RasterConfig(path=c.path, method=c.method,
                         resolution=c.resolution or res, factor=c.factor)
            for c in configs
        ]
    else:
        unresolved = [c for c in configs if c.resolution is None]
        if unresolved:
            print(
                f"[ERROR] Surface has multiple resolutions {resolutions_in_surface}. "
                "Each raster must specify a resolution.",
                file=sys.stderr,
            )
            return 2

    resolution_groups = group_by_resolution(configs)

    print(f"Target: {target}")
    print(f"Committing {len(configs)} raster(s)\n")

    for res, raster_cfgs in sorted(resolution_groups.items()):
        lod = lod_name(res)
        if res not in surface_groups:
            print(f"  [WARN] No surface cells at resolution {res} — skipping")
            continue

        print(f"── {lod} (resolution {res}) ─────────────────────────────")

        t_lod = time.perf_counter()
        surface = surface_groups[res]
        print(f"  Using {len(surface)} pre-built H3 cells")

        result = _process_lod_rasters(surface, raster_cfgs)

        if result.empty:
            print(f"  [WARN] No cells with data for {lod} — skipping")
            continue

        if lod in existing_layers:
            print(f"  Merging into existing {lod}...")
            existing_gdf = _read_from_gdb(target, lod)
            merged = merge_lod(existing_gdf, result)
            action = "merge"
        else:
            print(f"  Creating new feature class: {lod}")
            merged = result
            action = "create"

        _write_to_gdb(merged, target, lod)
        existing_layers.add(lod)

        elapsed = time.perf_counter() - t_lod
        print(f"  [OK] {lod}: {len(merged)} cells, "
              f"{len(_raster_data_columns(merged))} data fields ({elapsed:.2f}s)")

        for cfg in raster_cfgs:
            log_entries.append(
                _make_log_entry(
                    raster_path=str(cfg.path),
                    method=cfg.method,
                    resolution=res,
                    lod=lod,
                    cell_count=len(merged),
                    action=action,
                    factor=cfg.factor,
                    source_surface=str(surface_path),
                )
            )

    _append_log_entries(target, log_entries)
    print(f"\n[OK] Committed {len(log_entries)} raster(s) to {target.name}")
    return 0


# Subcommand: info

def cmd_info(args: argparse.Namespace) -> int:
    """Display Surface Model metadata and commit history."""
    target = Path(args.target).expanduser().resolve()

    if not target.exists():
        print(f"[ERROR] Not found: {target}", file=sys.stderr)
        return 2

    if not is_surface_model(target):
        print(f"[ERROR] {target.name} is not a Hex Surface Model.", file=sys.stderr)
        return 2

    meta_gdf = _read_from_gdb(target, META_TABLE)
    row = meta_gdf.iloc[0]

    print(f"Hex Surface Model: {target.name}")
    print(f"Created:           {row.get('created_at', 'N/A')}")
    print(f"Schema version:    {row.get('schema_version', 'N/A')}")
    print()

    # LOD summary
    all_layers = _gdb_layer_names(target)
    lod_layers = sorted(
        [(name, parse_lod_name(name)) for name in all_layers if parse_lod_name(name) is not None],
        key=lambda x: x[1],
    )

    if lod_layers:
        print(f"Feature Classes: ({len(lod_layers)} LODs)")
        print(f"  {'Name':<14} {'Res':>4}  {'Cells':>8}  {'Data Fields'}")
        print(f"  {'─' * 13} {'─' * 4} {'─' * 9}  {'─' * 40}")

        for layer_name, res in lod_layers:
            gdf = _read_from_gdb(target, layer_name)
            data_cols = _raster_data_columns(gdf)
            factor_info = []
            for col in sorted(data_cols):
                if col.startswith("F") and col[1:].isdigit():
                    fnum = int(col[1:])
                    fname = FACTOR_NAMES.get(fnum, "?")
                    factor_info.append(f"{col}({fname})")
                else:
                    factor_info.append(col)
            print(f"  {layer_name:<14} {res:>4}  {len(gdf):>8}  {', '.join(factor_info)}")
    else:
        print("Feature Classes: (none)")

    # Commit history
    print()
    if LOG_TABLE in all_layers:
        log_gdf = _read_from_gdb(target, LOG_TABLE)
        if not log_gdf.empty:
            print(f"Commit History: ({len(log_gdf)} entries)")
            print(f"  {'Timestamp':<28} {'LOD':<8} {'Method':<10} {'Cells':>7}"
                  f"  {'Action':<7}  {'Factor':<16}  {'Source'}")
            print(f"  {'─' * 27} {'─' * 7} {'─' * 9} {'─' * 7}"
                  f"  {'─' * 6}  {'─' * 15}  {'─' * 26}")

            for _, row in log_gdf.iterrows():
                ts = str(row.get("committed_at", ""))[:27]
                lod = str(row.get("lod", ""))
                method = str(row.get("method", ""))
                cells = int(row.get("cell_count", 0))
                action = str(row.get("action", ""))
                source = Path(str(row.get("raster_source", ""))).name

                factor_val = row.get("factor", -1)
                if factor_val and factor_val > 0:
                    factor_str = f"F{int(factor_val)}({FACTOR_NAMES.get(int(factor_val), '?')})"
                else:
                    factor_str = "(none)"

                print(f"  {ts:<28} {lod:<8} {method:<10} {cells:>7}"
                      f"  {action:<7}  {factor_str:<16}  {source}")
    else:
        print("Commit History: (none)")

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
            "Raster config. Standard: path:method:resolution[:factor]. "
            "Surface mode: path:method[:factor]. Repeatable. "
            "Factor: 1=Elevation, 2=Slope, 3=Vegetation, 4=Soil. "
            "Example: elevation.tif:mean:9:1"
        ),
    )
    p_commit.add_argument(
        "--surface",
        default=None,
        help=(
            "Path to a pre-built H3 surface vector file. "
            "Must contain GRID_ID and res columns."
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


def main(argv: List[str] | None = None) -> int:
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

    return handler(args)


if __name__ == "__main__":
    raise SystemExit(main())