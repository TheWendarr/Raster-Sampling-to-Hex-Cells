# Hex Surface Model

An open-source Python pipeline for translating raster data into H3 hexagonal cell surfaces and managing them as multi-resolution geospatial data fabrics. Built for Cross-Country Mobility (CCM) analysis and general-purpose terrain attribution.

## Overview

The pipeline takes raster inputs (elevation, slope, landcover, soil, etc.), tessellates their footprints into [H3 hexagonal cells](https://h3geo.org/) at a specified resolution, samples raster values into each cell, and stores the results in an Esri File Geodatabase (FGDB) organized by Level of Detail (LOD). Each LOD is a separate feature class containing only cells with valid data.

The end state is a **Hex Surface Model** — a multi-resolution attributed hex grid that can be incrementally updated as new data becomes available, with hierarchical pathfinding across the terrain.

## Factor Schema

Each raster committed to the model is tagged with a terrain factor:

| Factor | Column | Description |
|--------|--------|-------------|
| 1 | F1 | Elevation |
| 2 | F2 | Slope |
| 3 | F3 | Vegetation |
| 4 | F4 | Soil |

Raw values are stored as-is. Normalization to 0.0–1.0 is deferred to pathfinding time, enabling different vehicle profiles to interpret the same terrain data differently.

**Mobility Formula:** `Mobility = (F1 / F2) * F3 * F4`

Where each factor is normalized to 0.0 (impassable) – 1.0 (no restriction).

## Project Structure

```
geo_io.py                    # Geospatial I/O library (read/write any vector format)
generate_h3_tesselation.py   # H3 cell surface generation from vector extents
sample_raster.py             # Raster zonal statistics into polygon bins
build_surface.py             # Surface Model orchestrator (init / commit / info)
find_path.py                 # Hierarchical least-cost pathfinding
```

## Dependencies

```
pip install geopandas pyogrio h3 rasterio shapely numpy
```

Optional: `exactextract` (for faster raster sampling), `polars` (for attribute-only DataFrame operations via `geo_io.to_polars()`).

FGDB write support requires GDAL >= 3.6 (bundled with recent pyogrio releases or available via conda-forge).

## Quick Start

```bash
# Initialize a new Surface Model
python build_surface.py init --target model.gdb

# Commit rasters with factor tags
python build_surface.py commit --target model.gdb \
  --raster elevation.tif:mean:7:1 \
  --raster slope.tif:mean:7:2 \
  --raster landcover.tif:majority:7:3 \
  --raster soil.tif:majority:7:4

# Add higher-resolution data for a subarea
python build_surface.py commit --target model.gdb \
  --raster fine_elevation.tif:mean:9:1 \
  --raster fine_slope.tif:mean:9:2

# Check what's in the model
python build_surface.py info --target model.gdb

# Find path between two points
python find_path.py \
  --surface model.gdb \
  --start 39.75,-105.0 \
  --end 39.60,-104.8 \
  --output results.gdb \
  --layer best_route
```

## Scripts

### `build_surface.py`

The primary orchestrator. Manages Hex Surface Model FGDBs through three subcommands.

**`init`** — Create a new FGDB tagged as a Hex Surface Model.

```bash
python build_surface.py init --target model.gdb
python build_surface.py init --target model.gdb --overwrite
```

**`commit`** — Tessellate raster footprints, sample values, and write to the target FGDB. Each raster is specified as `path:method:resolution:factor`.

```bash
# Factor-tagged commit (new schema)
python build_surface.py commit --target model.gdb \
  --raster elevation.tif:mean:9:1 \
  --raster slope.tif:mean:9:2 \
  --raster landcover.tif:majority:9:3 \
  --raster soil.tif:majority:9:4

# Legacy (no factor tag — column named after raster)
python build_surface.py commit --target model.gdb \
  --raster elevation.tif:mean:9

# Using a pre-built H3 surface
python build_surface.py commit --target model.gdb \
  --surface h3_cells.geojson \
  --raster elevation.tif:mean:1 \
  --raster slope.tif:mean:2
```

The factor tag is the 4th component (1-4). When provided, the sampled column is named `F1`–`F4` instead of `{raster_name}_{method}`. This enables `find_path.py` to locate terrain factors by convention.

Merge semantics (last-write-wins): re-committing a raster with the same factor tag overwrites the previous values for overlapping cells.

**`info`** — Display model metadata, LOD summary, and full commit history.

```bash
python build_surface.py info --target model.gdb
```

```
Hex Surface Model: model.gdb
Created:           2026-04-05T14:32:00+00:00
Schema version:    2

Feature Classes: (2 LODs)
  Name          Res      Cells  Data Fields
  ──────────── ──── ──────────  ────────────────────────────────────────
  LOD_07          7          5  F1(Elevation), F2(Slope), F3(Vegetation), F4(Soil)
  LOD_09          9        211  F1(Elevation), F2(Slope)

Commit History: (6 entries)
  Timestamp                    LOD      Method       Cells  Action   Factor            Source
  ──────────────────────────── ──────── ────────── ───────  ───────  ───────────────── ──────
  2026-04-05 14:32             LOD_07   mean             5  create   F1(Elevation)     elevation.tif
  2026-04-05 14:32             LOD_07   mean             5  create   F2(Slope)         slope.tif
  2026-04-05 14:32             LOD_07   majority         5  create   F3(Vegetation)    landcover.tif
  2026-04-05 14:32             LOD_07   majority         5  create   F4(Soil)          soil.tif
  2026-04-05 14:35             LOD_09   mean           211  create   F1(Elevation)     fine_elev.tif
  2026-04-05 14:35             LOD_09   mean           211  create   F2(Slope)         fine_slope.tif
```

### `find_path.py`

Hierarchical least-cost pathfinding over the Hex Surface Model.

```bash
# Auto-select farthest points
python find_path.py --surface model.gdb --output path.geojson

# Specify start/end as lat,lon
python find_path.py \
  --surface model.gdb \
  --start 39.75,-105.0 \
  --end 39.60,-104.8 \
  --output results.gdb --layer route_1

# Specify start/end as H3 cell indices
python find_path.py \
  --surface model.gdb \
  --start 872a1008003ffff \
  --end 872a1008007ffff \
  --output results.gdb
```

**How hierarchical refinement works:**

1. Start at the coarsest LOD in the model (e.g., LOD_07)
2. Build an H3 neighbor graph, compute mobility per cell, run A*
3. For each cell in the coarse path, check if child cells exist at the next finer LOD (e.g., LOD_09)
4. Where child data exists, run A* at the finer resolution between entry/exit points
5. Repeat recursively until no finer data is available
6. Export the final mixed-resolution path as a line + cell points

**Output:** Writes to FGDB (line layer + cell points layer) or GeoJSON/Shapefile.

**Safety gate:** Refuses to write into a Surface Model FGDB to prevent contamination.

### `geo_io.py`

Unified geospatial I/O library. Reads and writes GeoJSON, Shapefile, KML, and FGDB through a single interface.

### `generate_h3_tesselation.py`

Generates H3 cell surfaces from vector extents at a specified resolution.

### `sample_raster.py`

Samples raster values into polygon bins using zonal statistics (mean or majority).

## Supported Formats

| Format | Read | Write | Extension | Notes |
|--------|------|-------|-----------|-------|
| GeoJSON | Yes | Yes | `.geojson`, `.json` | Default single-layer format |
| Shapefile | Yes | Yes | `.shp` | 10-char field name limit applies |
| KML | Yes | Yes | `.kml` | Auto-reprojects to EPSG:4326 |
| FGDB | Yes | Yes | `.gdb` | Requires GDAL >= 3.6 for write |

## Roadmap

- [x] H3 tessellation from vector extents
- [x] Raster zonal statistics into hex cells
- [x] Surface Model FGDB management (init/commit/info)
- [x] Factor tagging (F1-F4) with raw value storage
- [x] Hierarchical least-cost pathfinding (coarse-to-fine refinement)
- [ ] Vehicle platform profiles (per-factor normalization curves)
- [ ] Multi-resolution flattened authoritative model
- [ ] Vector data overlays (hydrology, bridges, roads)
- [ ] MCOO generation per vehicle platform category
- [ ] GUI wrapper for all CLI commands