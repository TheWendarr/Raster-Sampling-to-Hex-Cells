# Raster-Sampling-to-Hex-Cells
An open source way to translate raster data to hexagon cells within the same CRS and extent using Python. Inputs are projected Rasters alongside GeoJSON of your Hex Surface.

# To-Do List:
- Modify Funtions to pass GeoDataFrames instead of geoJSON
- Create variables to determine highest priority of data
- Create pipeline for Flattened Authoritative Model
- Create method to overwrite lower fidelity / outdated Cells

- Incorporate Separate Repo (fgdb-tools) into workflow for import / export ability

- (Maybe) translate to C to compile a .exe file for quicker processing

## CLI Usage
Sample Rasters to Polygons
```bash
python sample_raster.py \
  --in-path path/to/HEXSURFACE.geojson \
  --out-path path/to/OUTPUT.geojson \
  --mean path/to/dem_elevation.tif \
  --mean path/to/slope_degrees.tif \
  --majority path/to/land_use_land_cover.tif \
  --threshold 0.1
```

Create H3 Polygons given input polygon (.geojson)
Will add the ability to use any file with extent later instead of limiting to geojson
```bash
python generate_h3_tesselation.py \
  --in path/to/bounding/polygon.geojson \
  --res int 0-15 for H3 Level of Detail \
  --out "path/to/output.geojson"
```

