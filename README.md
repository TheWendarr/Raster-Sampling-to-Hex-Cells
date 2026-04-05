# Hex Surface Model

An open-source Python pipeline for translating raster data into H3 hexagonal cell surfaces and managing them as multi-resolution geospatial data fabrics. Built for Cross-Country Mobility (CCM) analysis and general-purpose terrain attribution.

## Overview

The pipeline takes raster inputs (elevation, slope, landcover, soil, etc.), tessellates their footprints into [H3 hexagonal cells](https://h3geo.org/) at a specified resolution, samples raster values into each cell, and stores the results in an Esri File Geodatabase (FGDB) organized by Level of Detail (LOD). Each LOD is a separate feature class containing only cells with valid data.

The end state is a **Hex Surface Model** — a multi-resolution attributed hex grid that can be incrementally updated as new data becomes available.

## Project Structure

```
geo_io.py                    # Geospatial I/O library (read/write any vector format)
generate_h3_tesselation.py   # H3 cell surface generation from vector extents
sample_raster.py             # Raster zonal statistics into polygon bins
build_surface.py             # Surface Model orchestrator (init / commit / info)
```

## Dependencies

```
pip install geopandas pyogrio h3 rasterio shapely numpy exactextract scipy
```

`exactextract` provides 20–50× faster zonal statistics via a C++ engine. If not installed, the pipeline falls back to a batch-rasterize approach (still 10–20× faster than per-polygon loops). `scipy` is used by the fallback engine for grouped statistics.

Optional: `polars` (for attribute-only DataFrame operations via `geo_io.to_polars()`).

FGDB write support requires GDAL >= 3.6 (bundled with recent pyogrio releases or available via conda-forge).

## Quick Start

```bash
# Initialize a new Surface Model
python build_surface.py init --target model.gdb

# Commit rasters (format: path:method:resolution)
python build_surface.py commit --target model.gdb \
  --raster test_data/elevation.tif:mean:9 \
  --raster test_data/landcover.tif:majority:9

# Or use a pre-built H3 surface (skip tessellation)
python build_surface.py commit --target model.gdb \
  --surface my_surface.geojson \
  --raster test_data/elevation.tif:mean \
  --raster test_data/landcover.tif:majority

# Check what's in the model
python build_surface.py info --target model.gdb
```

## Scripts

### `build_surface.py`

The primary orchestrator. Manages Hex Surface Model FGDBs through three subcommands.

**`init`** — Create a new FGDB tagged as a Hex Surface Model. The tag is an internal metadata table (`__hex_surface_meta__`) that prevents accidental writes to unrelated FGDBs.

```bash
python build_surface.py init --target model.gdb
python build_surface.py init --target model.gdb --overwrite   # recreate from scratch
```

**`commit`** — Tessellate raster footprints, sample values, and write to the target FGDB. Each raster is specified as `path:method:resolution` where method is `mean` (continuous data) or `majority` (categorical data) and resolution is the H3 level (0–15). Multiple rasters can target the same or different LODs in a single command.

```bash
# Single LOD
python build_surface.py commit --target model.gdb \
  --raster elevation.tif:mean:9 \
  --raster slope.tif:mean:9 \
  --raster landcover.tif:majority:9

# Multiple LODs in one command
python build_surface.py commit --target model.gdb \
  --raster coarse_dem.tif:mean:7 \
  --raster fine_dem.tif:mean:12

# Add new data to an existing LOD (merges automatically)
python build_surface.py commit --target model.gdb \
  --raster soil_type.tif:majority:9
```

**`commit --surface`** — Use a pre-built H3 surface instead of tessellating from raster footprints. The surface file must contain `GRID_ID` and `res` columns. Resolution is read from the surface, so raster configs use the shorter `path:method` format.

```bash
# Generate a surface first
python generate_h3_tesselation.py \
  --in custom_aoi.shp --res 9 --out my_surface.geojson

# Commit rasters using the pre-built surface (resolution inferred)
python build_surface.py commit --target model.gdb \
  --surface my_surface.geojson \
  --raster elevation.tif:mean \
  --raster landcover.tif:majority

# You can still specify resolution explicitly to override
python build_surface.py commit --target model.gdb \
  --surface my_surface.geojson \
  --raster elevation.tif:mean:9

# Multi-resolution surfaces are supported — rasters without a resolution
# are applied to ALL resolutions in the surface
python build_surface.py commit --target model.gdb \
  --surface multi_res_surface.gdb --surface-layer cells \
  --raster elevation.tif:mean
```

When committing to an existing LOD, the merge logic joins on `GRID_ID`:
- **Shared cells, shared fields**: new values overwrite old values (where non-null).
- **Shared cells, new fields**: new columns are appended.
- **New cells**: appended as new rows.
- **Existing-only cells**: preserved with nulls for new-only columns.

**`info`** — Display model metadata, LOD summary, and full commit history.

```bash
python build_surface.py info --target model.gdb
```

```
Hex Surface Model: model.gdb
Created:           2026-04-04T14:32:00+00:00
Schema version:    1

Feature Classes: (2 LODs)
  Name          Res      Cells  Data Fields
  ──────────── ──── ──────────  ────────────────────────────────────────
  LOD_07          7          5  elevation_mean
  LOD_09          9        211  elevation_mean, landcover_majority, slope_mean

Commit History: (4 entries)
  Timestamp                    LOD      Method       Cells  Action   Source
  ──────────────────────────── ──────── ────────── ───────  ───────  ──────────────────────────
  2026-04-04 14:32             LOD_09   mean           211  create   elevation.tif
  2026-04-04 14:32             LOD_09   majority       211  create   landcover.tif
  2026-04-04 14:33             LOD_09   mean           211  merge    slope.tif
  2026-04-04 14:35             LOD_07   mean             5  create   coarse_dem.tif
```

### `geo_io.py`

Unified geospatial I/O library. Reads and writes GeoJSON, Shapefile, KML, and FGDB through a single interface backed by pyogrio with Arrow acceleration and a GeoPandas/Fiona fallback. All other scripts import from this module.

```python
from geo_io import read_vector, write_vector, list_layers

gdf = read_vector("input.shp", epsg=4326)
gdf = read_vector("data.gdb", layer="roads")

write_vector(gdf, "output.geojson", overwrite=True)
write_vector(gdf, "output.gdb", layer="results")

layers = list_layers("data.gdb", spatial_only=True)
```

CLI for quick conversions:

```bash
python geo_io.py --in input.shp --out output.geojson
python geo_io.py --in data.gdb --layer roads --out roads.shp --epsg 32613
python geo_io.py --in data.gdb --list-layers
```

### `generate_h3_tesselation.py`

Generates an H3 cell surface from any vector extent. Outputs a GeoDataFrame (or file) with one row per cell, carrying `GRID_ID` and `res` attributes.

Cell boundaries are batch-converted via `h3.cells_to_geo()` (single C call) instead of individual `cell_to_boundary()` calls. Multi-polygon inputs are tessellated in parallel using `concurrent.futures`.

```bash
python generate_h3_tesselation.py \
  --in extent.shp \
  --res 9 \
  --out h3_surface.geojson \
  --overwrite
```

Importable for in-memory use:

```python
from generate_h3_tesselation import build_h3_surface

surface = build_h3_surface(extent_gdf, res=9)
```

### `sample_raster.py`

Samples raster values into polygon bins using zonal statistics. Supports `mean` for continuous rasters and `majority` for categorical rasters with deterministic tie-breaking.

Uses a tiered engine strategy for performance:
1. **exactextract** (preferred) — C++ vectorized zonal stats, 20–50× faster
2. **Batch rasterize** (fallback) — single rasterize call + grouped NumPy ops, 10–20× faster
3. Both engines produce identical results; exactextract just does it faster

```bash
python sample_raster.py \
  --in-path h3_surface.geojson \
  --out-path attributed.geojson \
  --mean elevation.tif \
  --mean slope.tif \
  --majority landcover.tif \
  --threshold 0.1 \
  --overwrite
```

Importable for in-memory use:

```python
from sample_raster import RasterBinEnricher

enricher = RasterBinEnricher(surface_gdf)
enricher.process_raster("elevation.tif", method="mean")
enricher.process_raster("landcover.tif", method="majority")
result = enricher.result
```

## In-Memory Pipeline

All scripts are designed to chain without intermediate files:

```python
from geo_io import read_vector, write_vector
from generate_h3_tesselation import build_h3_surface
from sample_raster import RasterBinEnricher

extent = read_vector("aoi.shp", epsg=4326)
surface = build_h3_surface(extent, res=9)

enricher = RasterBinEnricher(surface)
enricher.process_raster("elevation.tif", method="mean")
enricher.process_raster("landcover.tif", method="majority")

write_vector(enricher.result, "attributed.gdb", layer="LOD_09")
```

Pre-built surfaces can skip tessellation entirely:

```python
from geo_io import read_vector, write_vector
from sample_raster import RasterBinEnricher

# Load a surface generated earlier or from an external source
surface = read_vector("my_h3_cells.geojson", epsg=4326)

enricher = RasterBinEnricher(surface)
enricher.process_raster("elevation.tif", method="mean")
enricher.process_raster("landcover.tif", method="majority")

write_vector(enricher.result, "attributed.gdb", layer="LOD_09")
```

## Supported Formats

| Format | Read | Write | Extension | Notes |
|--------|------|-------|-----------|-------|
| GeoJSON | Yes | Yes | `.geojson`, `.json` | Default single-layer format |
| Shapefile | Yes | Yes | `.shp` | 10-char field name limit applies |
| KML | Yes | Yes | `.kml` | Auto-reprojects to EPSG:4326; multi-layer supported |
| FGDB | Yes | Yes | `.gdb` | Requires GDAL >= 3.6 for write; multi-layer via `--layer` |

## Roadmap

- [ ] Vehicle platform mobility parameters (light/medium/heavy wheeled and tracked)
- [ ] Dijkstra / A* pathfinding over the Hex Surface Model
- [ ] Least-cost path generation with top-N branching from critical nodes
- [ ] Vector data overlays (hydrology, bridges, roads)
- [ ] MCOO generation per vehicle platform category
- [ ] GUI for the full pipeline