"""
find_path.py

Hierarchical least-cost pathfinding over a Hex Surface Model FGDB.

Finds the optimal path from point A to point B across an H3 hex grid
using the four terrain factors stored in the Surface Model:

    F1 — Elevation   (raw meters)
    F2 — Slope       (raw degrees)
    F3 — Vegetation  (raw classification value)
    F4 — Soil        (raw classification value)

    Mobility Factor = (F1 / F2) * F3 * F4

All factors are normalized to 0.00–1.00 before computing mobility,
where 1.0 = no restriction and 0.0 = impassable. Higher mobility
means faster traversal; the cost per edge is distance / speed,
where speed = max_speed * mobility.

Hierarchical Refinement
-----------------------
The pathfinder starts at the coarsest LOD available in the Surface
Model, finds a path using A*, then "zooms in" by expanding each cell
in the path into its H3 children at the next finer LOD (if data exists).
This repeats recursively until no finer data is available, producing
an increasingly detailed path.

Vehicle Profiles
----------------
Vehicle-specific normalization is stubbed. The current implementation
uses a generic normalization that treats all factors equally. The
architecture has clear hooks for per-vehicle factor mappings.

Usage:
    # Basic pathfinding (auto-selects farthest points)
    python find_path.py \\
        --surface model.gdb \\
        --output results.gdb \\
        --layer path_result

    # Specify start/end as lat,lon
    python find_path.py \\
        --surface model.gdb \\
        --start 39.75,-105.0 \\
        --end 39.60,-104.8 \\
        --output results.gdb \\
        --layer path_result

    # Specify start/end as H3 cell indices
    python find_path.py \\
        --surface model.gdb \\
        --start 892a1008003ffff \\
        --end 892a1008007ffff \\
        --output results.gdb \\
        --layer path_result

    # Export as GeoJSON instead
    python find_path.py \\
        --surface model.gdb \\
        --output path.geojson

Requires:
    geopandas, pyogrio, shapely, h3, numpy
"""

import argparse
import heapq
import math
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Set, Tuple

import geopandas as gpd
import h3
import numpy as np
from shapely.geometry import LineString, Point

from geo_io import list_layers, read_vector, write_vector
from profiles import ProfileEngine

try:
    import pyogrio

    _HAS_PYOGRIO = True
except ImportError:
    pyogrio = None
    _HAS_PYOGRIO = False


# Constants

META_TABLE = "__hex_surface_meta__"
MAGIC_ID = "HEX_SURFACE_MODEL"
LOD_PREFIX = "LOD_"

# Factor column names as stored in the Surface Model
FACTOR_COLS = ["F1", "F2", "F3", "F4"]

# Default max traversal speed (m/s) — ~30 mph, conservative cross-country
DEFAULT_MAX_SPEED_MS = 13.4

# Minimum speed floor to preserve connectivity (m/s)
# ~0.1 mph — extremely slow but not impassable
MIN_SPEED_MS = 0.05

# Factors below this threshold after normalization are treated as impassable
IMPASSABLE_THRESHOLD = 0.01


# Vehicle Profile (stub for future expansion)

@dataclass(frozen=True)
class VehicleProfile:
    """
    Vehicle platform parameters for factor normalization.

    Currently a stub — all factors use generic normalization.
    Future versions will define per-factor min/max mappings and
    threshold gates (e.g., "can't traverse slopes > 30°").
    """
    name: str = "Generic"
    max_speed_ms: float = DEFAULT_MAX_SPEED_MS

    # Factor normalization hooks (stubbed)
    # Each would be a callable: raw_value -> normalized 0.0-1.0
    # For now, normalization is computed from the data range.


DEFAULT_VEHICLE = VehicleProfile()


# Utilities

def haversine_m(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """Great-circle distance in meters between two WGS-84 points."""
    r = 6_371_000.0
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlmb = math.radians(lon2 - lon1)
    a = (math.sin(dphi / 2.0) ** 2
         + math.cos(p1) * math.cos(p2) * math.sin(dlmb / 2.0) ** 2)
    return 2.0 * r * math.asin(math.sqrt(a))


def parse_lod_name(name: str) -> Optional[int]:
    """LOD feature class name -> H3 resolution, or None."""
    if not name.startswith(LOD_PREFIX):
        return None
    try:
        return int(name[len(LOD_PREFIX):])
    except ValueError:
        return None


def is_surface_model(gdb_path: Path) -> bool:
    """Check whether a FGDB is tagged as a Hex Surface Model."""
    if not gdb_path.exists():
        return False
    layers = {rec["name"] for rec in list_layers(gdb_path) if rec["name"]}
    if META_TABLE not in layers:
        return False
    try:
        meta = read_vector(gdb_path, layer=META_TABLE)
        if meta.empty:
            return False
        return str(meta.iloc[0].get("magic_id", "")) == MAGIC_ID
    except Exception:
        return False


# Surface Model Loading

@dataclass
class LODData:
    """Data for a single Level of Detail."""
    resolution: int
    gdf: gpd.GeoDataFrame
    # Pre-computed lookups (populated during graph build)
    cell_data: Dict[str, Dict[str, float]] = field(default_factory=dict)
    coords: Dict[str, Tuple[float, float]] = field(default_factory=dict)
    cell_set: Set[str] = field(default_factory=set)


def load_surface_model(gdb_path: Path) -> Dict[int, LODData]:
    """
    Load all LOD feature classes from a Surface Model FGDB.

    Returns a dict mapping resolution -> LODData, sorted coarsest to finest.
    """
    layers = {rec["name"] for rec in list_layers(gdb_path) if rec["name"]}

    lod_map: Dict[int, LODData] = {}

    for layer_name in sorted(layers):
        res = parse_lod_name(layer_name)
        if res is None:
            continue

        gdf = read_vector(gdb_path, layer=layer_name)

        if "GRID_ID" not in gdf.columns:
            print(f"  [WARN] {layer_name} missing GRID_ID — skipping")
            continue

        # Check which factor columns are present
        present_factors = [c for c in FACTOR_COLS if c in gdf.columns]
        if not present_factors:
            print(f"  [WARN] {layer_name} has no factor columns (F1-F4) — skipping")
            continue

        print(f"  Loaded {layer_name}: {len(gdf)} cells, "
              f"factors: {', '.join(present_factors)}")

        lod_map[res] = LODData(resolution=res, gdf=gdf)

    return lod_map


# Factor Normalization

def normalize_factors(
    lod: LODData,
    vehicle: VehicleProfile = DEFAULT_VEHICLE,
    engine: Optional[ProfileEngine] = None,
) -> None:
    """
    Compute normalized factor values (0.0–1.0) for all cells in a LOD.

    Populates lod.cell_data with per-cell normalized factors and
    lod.coords / lod.cell_set for graph operations.

    If a ProfileEngine is provided, normalization uses the loaded factor
    schemas and vehicle profile (schema-based chain). Otherwise, falls
    back to generic min/max normalization from the data range.
    """
    gdf = lod.gdf

    if engine is not None:
        # ── Schema-based normalization via ProfileEngine ─────────────
        for _, row in gdf.iterrows():
            grid_id = str(row["GRID_ID"])
            lat, lon = h3.cell_to_latlng(grid_id)

            lod.coords[grid_id] = (lat, lon)
            lod.cell_set.add(grid_id)

            # Collect raw values
            raw = {}
            for col in FACTOR_COLS:
                if col in gdf.columns:
                    val = row.get(col)
                    # Keep raw value as-is (int, float, str, or NaN)
                    raw[col] = val
                # Missing columns are omitted — engine treats them as 1.0

            # Normalize through the engine chain
            lod.cell_data[grid_id] = engine.normalize_cell(raw)
        return

    # ── Fallback: generic min/max normalization ──────────────────────
    factor_ranges: Dict[str, Tuple[float, float]] = {}
    for col in FACTOR_COLS:
        if col not in gdf.columns:
            continue
        valid = gdf[col].dropna()
        if valid.empty:
            continue
        fmin, fmax = float(valid.min()), float(valid.max())
        factor_ranges[col] = (fmin, fmax)

    for _, row in gdf.iterrows():
        grid_id = str(row["GRID_ID"])
        lat, lon = h3.cell_to_latlng(grid_id)

        lod.coords[grid_id] = (lat, lon)
        lod.cell_set.add(grid_id)

        cell = {}
        for col in FACTOR_COLS:
            if col not in gdf.columns or col not in factor_ranges:
                cell[col] = 1.0
                continue

            raw = row.get(col)
            if raw is None or (isinstance(raw, float) and math.isnan(raw)):
                cell[col] = 0.0
                continue

            fmin, fmax = factor_ranges[col]
            if fmax == fmin:
                cell[col] = 1.0
            else:
                cell[col] = 1.0 - (float(raw) - fmin) / (fmax - fmin)

            cell[col] = max(0.0, min(1.0, cell[col]))

        lod.cell_data[grid_id] = cell


# Mobility Computation

def compute_mobility(cell: Dict[str, float]) -> float:
    """
    Compute the mobility factor for a single cell.

    Mobility = (F1 / F2) * F3 * F4

    All factors are normalized 0.0–1.0. To avoid division by zero,
    F2 is floored at IMPASSABLE_THRESHOLD.

    Returns a value in [0.0, ~1.0] range (can exceed 1.0 if F1 > F2).
    Values <= IMPASSABLE_THRESHOLD are treated as impassable.
    """
    f1 = cell.get("F1", 1.0)
    f2 = cell.get("F2", 1.0)
    f3 = cell.get("F3", 1.0)
    f4 = cell.get("F4", 1.0)

    # Gate: any factor at zero = impassable
    if f1 <= IMPASSABLE_THRESHOLD or f3 <= IMPASSABLE_THRESHOLD or f4 <= IMPASSABLE_THRESHOLD:
        return 0.0

    # Floor F2 to avoid division by zero
    f2_safe = max(f2, IMPASSABLE_THRESHOLD)

    mobility = (f1 / f2_safe) * f3 * f4

    return max(0.0, mobility)


def edge_cost(
    a_id: str,
    b_id: str,
    lod: LODData,
    vehicle: VehicleProfile = DEFAULT_VEHICLE,
) -> Optional[float]:
    """
    Compute traversal cost (time in seconds) for edge a->b.

    Returns None if either cell is impassable.
    """
    if a_id not in lod.cell_data or b_id not in lod.cell_data:
        return None

    mob_a = compute_mobility(lod.cell_data[a_id])
    mob_b = compute_mobility(lod.cell_data[b_id])

    # Average mobility across the edge
    avg_mob = (mob_a + mob_b) / 2.0

    if avg_mob <= IMPASSABLE_THRESHOLD:
        return None

    # Distance between cell centroids
    lat1, lon1 = lod.coords[a_id]
    lat2, lon2 = lod.coords[b_id]
    dist = haversine_m(lat1, lon1, lat2, lon2)

    if dist <= 0.0:
        return None

    # Effective speed: cap mobility at 1.0 for speed calculation
    speed = vehicle.max_speed_ms * min(avg_mob, 1.0)
    speed = max(speed, MIN_SPEED_MS)

    return dist / speed


# A* Pathfinding

def _heuristic(
    a_id: str,
    goal_id: str,
    coords: Dict[str, Tuple[float, float]],
    max_speed: float,
) -> float:
    """Admissible heuristic: straight-line time at max speed."""
    lat1, lon1 = coords[a_id]
    lat2, lon2 = coords[goal_id]
    return haversine_m(lat1, lon1, lat2, lon2) / max(max_speed, 0.1)


def astar(
    lod: LODData,
    start: str,
    goal: str,
    vehicle: VehicleProfile = DEFAULT_VEHICLE,
) -> Optional[List[str]]:
    """
    A* search over an H3 hex grid at a single LOD.

    Builds the neighbor graph lazily during search — no need to
    pre-compute all edges.

    Returns the path as a list of GRID_IDs, or None if no path exists.
    """
    if start not in lod.cell_set or goal not in lod.cell_set:
        return None

    coords = lod.coords
    max_speed = vehicle.max_speed_ms

    g_score: Dict[str, float] = {start: 0.0}
    came_from: Dict[str, str] = {}
    closed: Set[str] = set()

    f_start = _heuristic(start, goal, coords, max_speed)
    # Heap entries: (f_score, tie_breaker, node_id)
    counter = 0
    open_heap: List[Tuple[float, int, str]] = [(f_start, counter, start)]

    expansions = 0
    total = len(lod.cell_set)

    while open_heap:
        f_curr, _, current = heapq.heappop(open_heap)

        if current in closed:
            continue

        if current == goal:
            # Reconstruct path
            path = [current]
            while current in came_from:
                current = came_from[current]
                path.append(current)
            path.reverse()
            return path

        closed.add(current)
        expansions += 1

        if expansions % 500 == 0:
            print(f"    A* expanded {expansions}/{total} nodes...", end="\r")

        # Lazy neighbor discovery via H3 topology
        try:
            neighbors = h3.grid_disk(current, 1)
        except Exception:
            continue

        for nbr in neighbors:
            if nbr == current or nbr not in lod.cell_set or nbr in closed:
                continue

            cost = edge_cost(current, nbr, lod, vehicle)
            if cost is None:
                continue

            tentative = g_score[current] + cost
            if tentative < g_score.get(nbr, float("inf")):
                came_from[nbr] = current
                g_score[nbr] = tentative
                f_val = tentative + _heuristic(nbr, goal, coords, max_speed)
                counter += 1
                heapq.heappush(open_heap, (f_val, counter, nbr))

    return None  # No path found


# Hierarchical Refinement

def find_nearest_cell(
    lat: float,
    lon: float,
    lod: LODData,
) -> Optional[str]:
    """Find the cell in the LOD nearest to a given lat/lon."""
    target_cell = h3.latlng_to_cell(lat, lon, lod.resolution)
    if target_cell in lod.cell_set:
        return target_cell

    # If exact cell isn't in our data, find closest by distance
    best_id = None
    best_dist = float("inf")
    for cell_id in lod.cell_set:
        clat, clon = lod.coords[cell_id]
        d = haversine_m(lat, lon, clat, clon)
        if d < best_dist:
            best_dist = d
            best_id = cell_id

    return best_id


def find_farthest_pair(lod: LODData) -> Tuple[str, str]:
    """Find the two farthest cells in a LOD (by haversine distance)."""
    cells = list(lod.cell_set)
    if len(cells) <= 2:
        return (cells[0], cells[-1])

    # Sample if very large
    sample = cells if len(cells) <= 2000 else cells[::max(1, len(cells) // 2000)]

    best_dist = -1.0
    best_pair = (sample[0], sample[-1])

    for i, a in enumerate(sample):
        lat1, lon1 = lod.coords[a]
        for b in sample[i + 1:]:
            lat2, lon2 = lod.coords[b]
            d = haversine_m(lat1, lon1, lat2, lon2)
            if d > best_dist:
                best_dist = d
                best_pair = (a, b)

    return best_pair


def _build_scoped_lod(
    coarse_cells: List[str],
    fine_lod: LODData,
    fine_res: int,
) -> Optional[LODData]:
    """
    Build a lightweight LODData containing only the fine-resolution children
    of the given coarse cells that exist in the fine LOD.

    This scopes A* to a small neighborhood instead of the full dataset.
    """
    scoped_cells: Set[str] = set()

    for coarse_cell in coarse_cells:
        try:
            children = h3.cell_to_children(coarse_cell, fine_res)
        except Exception:
            continue
        scoped_cells.update(set(children) & fine_lod.cell_set)

    if not scoped_cells:
        return None

    # Build a mini LODData with only the scoped cells
    scoped = LODData(resolution=fine_res, gdf=fine_lod.gdf)  # gdf unused after normalization
    scoped.cell_set = scoped_cells
    scoped.coords = {cid: fine_lod.coords[cid] for cid in scoped_cells if cid in fine_lod.coords}
    scoped.cell_data = {cid: fine_lod.cell_data[cid] for cid in scoped_cells if cid in fine_lod.cell_data}

    return scoped


def refine_path(
    coarse_path: List[str],
    coarse_res: int,
    fine_lod: LODData,
    vehicle: VehicleProfile = DEFAULT_VEHICLE,
) -> List[str]:
    """
    Refine a coarse path by expanding each cell into its children at a finer LOD.

    For each consecutive pair of cells in the coarse path:
      1. Find all children of both cells that exist in the fine LOD
      2. Build a scoped mini-LOD containing only those children
      3. Identify entry/exit cells (children nearest to the coarse cell centroids)
      4. Run A* between them on the scoped LOD (fast — tens to hundreds of cells)
      5. Stitch the fine segments together

    Cells whose children don't exist in the fine LOD are kept at coarse
    resolution (their centroid is included in the path coordinates).
    """
    fine_res = fine_lod.resolution
    refined_segments: List[List[str]] = []
    segments_refined = 0
    segments_skipped = 0

    total_pairs = len(coarse_path) - 1
    print(f"    Refining {len(coarse_path)} coarse cells from res {coarse_res} -> {fine_res}")
    print(f"    ({total_pairs} segments to process)")

    for i in range(len(coarse_path)):
        coarse_cell = coarse_path[i]

        if i == 0:
            # First cell — check if it has children, add entry point
            try:
                children = set(h3.cell_to_children(coarse_cell, fine_res)) & fine_lod.cell_set
            except Exception:
                children = set()

            if children:
                lat, lon = h3.cell_to_latlng(coarse_cell)
                entry = _nearest_child(lat, lon, children, fine_lod)
                if entry:
                    refined_segments.append([entry])
            else:
                refined_segments.append([coarse_cell])
            continue

        # For each consecutive pair: prev_cell -> coarse_cell
        prev_cell = coarse_path[i - 1]

        # Get children of both cells
        try:
            prev_children = set(h3.cell_to_children(prev_cell, fine_res)) & fine_lod.cell_set
        except Exception:
            prev_children = set()

        try:
            curr_children = set(h3.cell_to_children(coarse_cell, fine_res)) & fine_lod.cell_set
        except Exception:
            curr_children = set()

        if not prev_children or not curr_children:
            # One side has no fine data — keep coarse cell
            refined_segments.append([coarse_cell])
            segments_skipped += 1
            continue

        # Build a scoped LOD with ONLY children of these two cells
        # (+ one ring of neighbors for connectivity)
        scope_parents = [prev_cell, coarse_cell]
        scoped_lod = _build_scoped_lod(scope_parents, fine_lod, fine_res)

        if scoped_lod is None or len(scoped_lod.cell_set) < 2:
            refined_segments.append([coarse_cell])
            segments_skipped += 1
            continue

        # Find entry point (child of prev nearest to the boundary between cells)
        # Use the midpoint between the two coarse centroids as the target
        prev_lat, prev_lon = h3.cell_to_latlng(prev_cell)
        curr_lat, curr_lon = h3.cell_to_latlng(coarse_cell)

        entry = _nearest_child(prev_lat, prev_lon, prev_children, fine_lod)
        exit_cell = _nearest_child(curr_lat, curr_lon, curr_children, fine_lod)

        if not entry or not exit_cell or entry == exit_cell:
            refined_segments.append([coarse_cell])
            segments_skipped += 1
            continue

        # Run A* on the SCOPED LOD (small! tens to hundreds of cells)
        segment = astar(scoped_lod, entry, exit_cell, vehicle)

        if segment:
            # Avoid duplicating the entry cell if it matches previous segment end
            if refined_segments and refined_segments[-1][-1] == segment[0]:
                segment = segment[1:]
            if segment:
                refined_segments.append(segment)
            segments_refined += 1
        else:
            # A* failed on scoped LOD — fall back to coarse
            refined_segments.append([coarse_cell])
            segments_skipped += 1

        # Progress
        if (i + 1) % 50 == 0 or i == total_pairs:
            print(f"    Processed {i + 1}/{total_pairs} segments "
                  f"(refined: {segments_refined}, skipped: {segments_skipped})", end="\r")

    print(f"    Processed {total_pairs}/{total_pairs} segments "
          f"(refined: {segments_refined}, skipped: {segments_skipped})")

    # Flatten segments into a single path
    flat: List[str] = []
    for seg in refined_segments:
        for cell_id in seg:
            if not flat or flat[-1] != cell_id:
                flat.append(cell_id)

    return flat if flat else coarse_path


def _nearest_child(
    lat: float,
    lon: float,
    children: Set[str],
    lod: LODData,
) -> Optional[str]:
    """Find the child cell nearest to a given point."""
    best_id = None
    best_dist = float("inf")
    for child in children:
        if child in lod.coords:
            clat, clon = lod.coords[child]
        else:
            clat, clon = h3.cell_to_latlng(child)
        d = haversine_m(lat, lon, clat, clon)
        if d < best_dist:
            best_dist = d
            best_id = child
    return best_id


# Main Pathfinding Pipeline

def find_path(
    surface_gdb: Path,
    start: Optional[str] = None,
    end: Optional[str] = None,
    vehicle: VehicleProfile = DEFAULT_VEHICLE,
    engine: Optional[ProfileEngine] = None,
) -> Tuple[List[str], Dict[str, Any]]:
    """
    Run hierarchical pathfinding over a Surface Model FGDB.

    Parameters
    ----------
    surface_gdb : Path
        Path to the Surface Model FGDB.
    start : str, optional
        Start point as "lat,lon" or an H3 cell index. If None, auto-selected.
    end : str, optional
        End point as "lat,lon" or an H3 cell index. If None, auto-selected.
    vehicle : VehicleProfile
        Vehicle parameters for normalization and speed.
    engine : ProfileEngine, optional
        If provided, uses schema-based normalization. Otherwise falls back
        to generic min/max normalization.

    Returns
    -------
    path : list[str]
        Ordered list of H3 cell IDs forming the path (mixed resolutions
        possible if refinement occurred).
    metadata : dict
        Path statistics: distance_m, time_s, cells, resolutions_used.
    """
    t_start = time.perf_counter()

    print(f"Loading Surface Model: {surface_gdb.name}")
    lod_map = load_surface_model(surface_gdb)

    if not lod_map:
        raise RuntimeError("No LOD feature classes with factor data found.")

    resolutions = sorted(lod_map.keys())
    print(f"Available LODs: {resolutions} (coarsest: {resolutions[0]}, finest: {resolutions[-1]})")

    # Normalize factors for all LODs
    if engine:
        print(f"Normalizing factors via {engine.vehicle.vehicle_name} profile...")
    else:
        print("Normalizing factors (generic min/max)...")
    for res in resolutions:
        normalize_factors(lod_map[res], vehicle, engine=engine)
        print(f"  LOD_{res:02d}: {len(lod_map[res].cell_set)} traversable cells")

    # Resolve start/end points
    coarsest = lod_map[resolutions[0]]

    start_id, end_id = _resolve_endpoints(start, end, coarsest)
    print(f"Start: {start_id}")
    print(f"End:   {end_id}")

    # Phase 1: Coarse pathfinding
    print(f"\nPhase 1: Pathfinding at LOD_{resolutions[0]:02d} "
          f"(resolution {resolutions[0]})...")
    t_phase1 = time.perf_counter()

    path = astar(coarsest, start_id, end_id, vehicle)

    if path is None:
        raise RuntimeError(
            f"No path found between {start_id} and {end_id} at "
            f"resolution {resolutions[0]}. The cells may not be connected."
        )

    print(f"  Found path: {len(path)} cells ({time.perf_counter() - t_phase1:.2f}s)")

    # Phase 2+: Hierarchical refinement
    current_res = resolutions[0]

    for finer_res in resolutions[1:]:
        print(f"\nPhase {resolutions.index(finer_res) + 1}: "
              f"Refining to LOD_{finer_res:02d} (resolution {finer_res})...")
        t_refine = time.perf_counter()

        refined = refine_path(path, current_res, lod_map[finer_res], vehicle)

        if len(refined) > len(path):
            print(f"  Refined: {len(path)} -> {len(refined)} cells "
                  f"({time.perf_counter() - t_refine:.2f}s)")
            path = refined
            current_res = finer_res
        else:
            print(f"  No refinement possible at resolution {finer_res} "
                  f"(no child data overlaps path)")

    # Compute path statistics
    total_dist = 0.0
    total_time = 0.0

    for i in range(len(path) - 1):
        a, b = path[i], path[i + 1]

        # Get coordinates — cells may be at different resolutions
        coord_a = _get_cell_coords(a, lod_map)
        coord_b = _get_cell_coords(b, lod_map)

        if coord_a and coord_b:
            d = haversine_m(*coord_a, *coord_b)
            total_dist += d

            # Estimate time using the finest LOD that has both cells
            cost = _get_edge_cost(a, b, lod_map, vehicle)
            total_time += cost if cost else d / MIN_SPEED_MS

    elapsed = time.perf_counter() - t_start

    metadata = {
        "distance_m": total_dist,
        "distance_km": total_dist / 1000.0,
        "time_s": total_time,
        "time_min": total_time / 60.0,
        "cells": len(path),
        "resolutions_used": resolutions,
        "elapsed_s": elapsed,
        "vehicle": engine.vehicle.vehicle_name if engine else vehicle.name,
    }

    print(f"\n{'═' * 50}")
    print(f"Path found: {len(path)} cells")
    print(f"Distance:   {total_dist / 1000.0:.2f} km")
    print(f"Est. time:  {total_time / 60.0:.1f} min")
    print(f"Computed in {elapsed:.2f}s")
    print(f"{'═' * 50}")

    return path, metadata


def _resolve_endpoints(
    start: Optional[str],
    end: Optional[str],
    lod: LODData,
) -> Tuple[str, str]:
    """Resolve start/end strings into H3 cell IDs within the given LOD."""
    start_id = _resolve_point(start, lod, "start") if start else None
    end_id = _resolve_point(end, lod, "end") if end else None

    if not start_id or not end_id:
        print("  Auto-selecting farthest points...")
        far_a, far_b = find_farthest_pair(lod)
        if not start_id:
            start_id = far_a
        if not end_id:
            end_id = far_b

    if start_id == end_id:
        raise ValueError("Start and end resolve to the same cell.")

    return start_id, end_id


def _resolve_point(spec: str, lod: LODData, label: str) -> Optional[str]:
    """
    Resolve a point spec to an H3 cell ID.

    Accepts:
      - H3 cell index string (e.g., "892a1008003ffff")
      - lat,lon string (e.g., "39.75,-105.0")
    """
    spec = spec.strip()

    # Try as H3 index first
    if h3.is_valid_cell(spec):
        # Check if it's in our data, or find the nearest
        if spec in lod.cell_set:
            return spec
        # It might be at a different resolution — find nearest
        lat, lon = h3.cell_to_latlng(spec)
        nearest = find_nearest_cell(lat, lon, lod)
        if nearest:
            print(f"  {label}: H3 index {spec} mapped to nearest cell {nearest}")
            return nearest
        return None

    # Try as lat,lon
    if "," in spec:
        parts = spec.split(",")
        if len(parts) == 2:
            try:
                lat = float(parts[0].strip())
                lon = float(parts[1].strip())
                nearest = find_nearest_cell(lat, lon, lod)
                if nearest:
                    print(f"  {label}: ({lat}, {lon}) mapped to cell {nearest}")
                    return nearest
            except ValueError:
                pass

    print(f"  [WARN] Could not resolve {label} '{spec}'", file=sys.stderr)
    return None


def _get_cell_coords(
    cell_id: str,
    lod_map: Dict[int, LODData],
) -> Optional[Tuple[float, float]]:
    """Get coordinates for a cell, checking all LODs."""
    for lod in lod_map.values():
        if cell_id in lod.coords:
            return lod.coords[cell_id]
    # Fallback: compute from H3
    try:
        return h3.cell_to_latlng(cell_id)
    except Exception:
        return None


def _get_edge_cost(
    a_id: str,
    b_id: str,
    lod_map: Dict[int, LODData],
    vehicle: VehicleProfile,
) -> Optional[float]:
    """Get edge cost from the finest LOD that contains both cells."""
    for res in sorted(lod_map.keys(), reverse=True):
        lod = lod_map[res]
        if a_id in lod.cell_data and b_id in lod.cell_data:
            return edge_cost(a_id, b_id, lod, vehicle)
    return None


# Output

def export_path(
    path: List[str],
    metadata: Dict[str, Any],
    output_path: Path,
    layer_name: Optional[str] = None,
    lod_map: Optional[Dict[int, LODData]] = None,
) -> None:
    """
    Export the path as a feature class.

    Creates two geometries:
      1. A LineString connecting all cell centroids (the route)
      2. Point features for each cell along the path with attributes

    The output format is inferred from the file extension.
    """
    # Build coordinate list
    coords_list: List[Tuple[float, float]] = []
    cell_records: List[Dict[str, Any]] = []

    for i, cell_id in enumerate(path):
        # Get coordinates
        coord = None
        if lod_map:
            coord = _get_cell_coords(cell_id, lod_map)
        if not coord:
            try:
                coord = h3.cell_to_latlng(cell_id)
            except Exception:
                continue

        lat, lon = coord
        coords_list.append((lon, lat))  # GeoJSON is lon,lat

        # Determine resolution of this cell
        try:
            cell_res = h3.get_resolution(cell_id)
        except Exception:
            cell_res = -1

        cell_records.append({
            "GRID_ID": cell_id,
            "sequence": i,
            "resolution": cell_res,
            "lat": lat,
            "lon": lon,
        })

    if len(coords_list) < 2:
        raise RuntimeError("Path has fewer than 2 valid coordinates; cannot export.")

    # Build line feature
    line = LineString(coords_list)
    line_gdf = gpd.GeoDataFrame(
        [{
            "path_type": "least_cost",
            "cells": len(path),
            "distance_km": round(metadata.get("distance_km", 0.0), 3),
            "time_min": round(metadata.get("time_min", 0.0), 1),
            "vehicle": metadata.get("vehicle", "Generic"),
        }],
        geometry=[line],
        crs="EPSG:4326",
    )

    # Build point features for cells along the path
    point_geoms = [Point(lon, lat) for lon, lat in coords_list]
    points_gdf = gpd.GeoDataFrame(
        cell_records,
        geometry=point_geoms,
        crs="EPSG:4326",
    )

    # Write output
    ext = output_path.suffix.lower()

    if ext == ".gdb":
        _require_pyogrio()
        line_layer = layer_name or "path_line"
        points_layer = f"{line_layer}_cells"

        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Check target is NOT a surface model
        if output_path.exists() and is_surface_model(output_path):
            raise ValueError(
                f"Cannot write path to a Surface Model FGDB: {output_path}\n"
                "Specify a different output GDB."
            )

        pyogrio.write_dataframe(
            line_gdf, str(output_path),
            layer=line_layer, driver="OpenFileGDB",
        )
        pyogrio.write_dataframe(
            points_gdf, str(output_path),
            layer=points_layer, driver="OpenFileGDB",
        )
        print(f"Exported to {output_path.name}:")
        print(f"  Line layer:   {line_layer}")
        print(f"  Points layer: {points_layer}")

    else:
        # GeoJSON, Shapefile, etc. — write line only
        write_vector(line_gdf, output_path, overwrite=True)
        print(f"Exported path to: {output_path}")

        # Also write points as a companion file
        points_path = output_path.with_stem(output_path.stem + "_cells")
        write_vector(points_gdf, points_path, overwrite=True)
        print(f"Exported cells to: {points_path}")


# CLI

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="find_path",
        description=(
            "Hierarchical least-cost pathfinding over a Hex Surface Model.\n\n"
            "Computes the optimal path from A to B using terrain factors\n"
            "(Elevation, Slope, Vegetation, Soil) with coarse-to-fine\n"
            "H3 hexagonal refinement."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--surface",
        required=True,
        help="Path to the Surface Model .gdb.",
    )
    parser.add_argument(
        "--start",
        default=None,
        help=(
            "Start point as 'lat,lon' or H3 cell index. "
            "If omitted, auto-selects the farthest pair."
        ),
    )
    parser.add_argument(
        "--end",
        default=None,
        help=(
            "End point as 'lat,lon' or H3 cell index. "
            "If omitted, auto-selects the farthest pair."
        ),
    )
    parser.add_argument(
        "--output",
        required=True,
        help=(
            "Output path. Use .gdb for File Geodatabase (writes line + cell layers), "
            "or .geojson/.shp for single-file export."
        ),
    )
    parser.add_argument(
        "--layer",
        default="path_result",
        help="Layer name for FGDB output (default: path_result).",
    )
    parser.add_argument(
        "--max-speed",
        type=float,
        default=DEFAULT_MAX_SPEED_MS,
        help=f"Max traversal speed in m/s (default: {DEFAULT_MAX_SPEED_MS:.1f}, ~30 mph). "
             "Overridden by --vehicle if provided.",
    )
    parser.add_argument(
        "--factor-schema",
        action="append",
        default=None,
        help=(
            "Path to a factor schema JSON (e.g., landfire_gap_2011.json, uscs_soils.json). "
            "Repeatable — one per categorical factor. "
            "Required when using --vehicle."
        ),
    )
    parser.add_argument(
        "--vehicle",
        default=None,
        help=(
            "Path to a vehicle profile JSON (e.g., light_wheeled.json). "
            "Defines per-factor normalization. Requires --factor-schema for "
            "categorical factors (F3 Vegetation, F4 Soil). "
            "If omitted, uses generic min/max normalization."
        ),
    )

    return parser


def main(argv: List[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    surface_path = Path(args.surface).expanduser().resolve()
    output_path = Path(args.output).expanduser().resolve()

    # Validate surface
    if not surface_path.exists():
        print(f"[ERROR] Surface Model not found: {surface_path}", file=sys.stderr)
        return 2

    if not is_surface_model(surface_path):
        print(
            f"[ERROR] {surface_path.name} is not a Hex Surface Model FGDB.",
            file=sys.stderr,
        )
        return 2

    # Validate output is not a surface model
    if output_path.exists() and output_path.suffix.lower() == ".gdb":
        if is_surface_model(output_path):
            print(
                f"[ERROR] Cannot write to Surface Model FGDB: {output_path}\n"
                "       Specify a different output GDB.",
                file=sys.stderr,
            )
            return 2

    vehicle = VehicleProfile(
        name="Generic",
        max_speed_ms=args.max_speed,
    )

    # Build ProfileEngine if schemas/vehicle are provided
    engine = None
    if args.vehicle or args.factor_schema:
        schema_paths = args.factor_schema or []

        # Validate schema files exist
        for sp in schema_paths:
            if not Path(sp).exists():
                print(f"[ERROR] Factor schema not found: {sp}", file=sys.stderr)
                return 2

        # Validate vehicle file exists
        vehicle_path = None
        if args.vehicle:
            vehicle_path = Path(args.vehicle)
            if not vehicle_path.exists():
                print(f"[ERROR] Vehicle profile not found: {vehicle_path}", file=sys.stderr)
                return 2

        print("Loading profiles...")
        engine = ProfileEngine(
            factor_schemas=schema_paths,
            vehicle_path=str(vehicle_path) if vehicle_path else None,
        )

        # Override vehicle max speed from the profile
        if engine.vehicle.max_speed_ms > 0:
            vehicle = VehicleProfile(
                name=engine.vehicle.vehicle_name,
                max_speed_ms=engine.vehicle.max_speed_ms,
            )

    try:
        path, metadata = find_path(
            surface_gdb=surface_path,
            start=args.start,
            end=args.end,
            vehicle=vehicle,
            engine=engine,
        )

        # Re-load LOD map for coordinate lookup during export
        lod_map = load_surface_model(surface_path)
        for lod in lod_map.values():
            normalize_factors(lod, vehicle, engine=engine)

        export_path(
            path=path,
            metadata=metadata,
            output_path=output_path,
            layer_name=args.layer,
            lod_map=lod_map,
        )

        return 0

    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())