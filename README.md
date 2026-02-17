# Raster-Sampling-to-Hex-Cells
An open source way to translate raster data to hexagon cells within the same CRS and extent using Python. Inputs are projected Rasters alongside GeoJSON of your Hex Surface.

## CLI Usage
```bash
python sample_raster.py \
  --in-path path/to/HEXSURFACE.geojson \
  --out-path path/to/OUTPUT.geojson \
  --mean path/to/dem_elevation.tif \
  --mean path/to/slope_degrees.tif \
  --majority path/to/land_use_land_cover.tif \
  --threshold 0.1
```
