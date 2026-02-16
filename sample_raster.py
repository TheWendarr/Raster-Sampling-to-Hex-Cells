import os
import datetime
import warnings
from typing import Union, Literal, List, Dict

import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio import features
from rasterio.windows import from_bounds
from shapely.geometry import box

# Suppress specific rasterio warnings
warnings.filterwarnings("ignore", category=rasterio.errors.NotGeoreferencedWarning)

class RasterBinEnricher:
    """
    A tool to enrich H3/Polygon bins with Raster statistics.
    Handles 'Mean' for continuous data and 'Majority' for discrete data.
    """

    VERSION = "0.1.0"

    def __init__(self, geojson_path: str):
        print(f"Loading Vector Data: {geojson_path}")
        self.gdf = gpd.read_file(geojson_path)
        self.original_crs = self.gdf.crs

    def process_raster(self, raster_path: str, method: Literal['mean', 'majority'], nodata_threshold: float = 0.1):
        """
        Enriches the internal GeoDataFrame with statistics from the provided raster.
        """
        if not os.path.exists(raster_path):
            print(f"Skipping: File not found -> {raster_path}")
            return

        raster_name = os.path.splitext(os.path.basename(raster_path))[0]
        output_col = f"{raster_name}_{method}"
        
        print(f"Processing Raster: {raster_name} | Method: {method.upper()}")

        results = []
        
        with rasterio.open(raster_path) as src:
            # Check CRS compatibility (warn only)
            if self.original_crs and src.crs != self.original_crs:
                print(f"   Warning: CRS mismatch. Vector: {self.original_crs}, Raster: {src.crs}")

            nodata_val = src.nodata

            for idx, row in self.gdf.iterrows():
                geom = row.geometry
                
                # 1. Spatial Windowing
                minx, miny, maxx, maxy = geom.bounds
                window = from_bounds(minx, miny, maxx, maxy, src.transform)
                
                try:
                    data = src.read(1, window=window)
                    win_transform = src.window_transform(window)
                except Exception:
                    results.append(None)
                    continue

                # 2. Rasterize Geometry (Mask)
                mask_shape = data.shape
                if mask_shape[0] == 0 or mask_shape[1] == 0:
                    results.append(None)
                    continue
                    
                mask_image = features.rasterize(
                    [(geom, 1)],
                    out_shape=mask_shape,
                    transform=win_transform,
                    fill=0,
                    dtype='uint8',
                    all_touched=False # Enforce Center of Pixel rule
                )
                
                # 3. Extract Valid Data
                flat_data = data.flatten()
                flat_mask = mask_image.flatten()
                
                pixels_in_poly = flat_data[flat_mask == 1]
                
                if pixels_in_poly.size == 0:
                    results.append(None)
                    continue

                if nodata_val is not None:
                    valid_pixels = pixels_in_poly[pixels_in_poly != nodata_val]
                else:
                    valid_pixels = pixels_in_poly
                
                # 4. Check NoData Threshold
                total_poly_pixels = pixels_in_poly.size
                valid_count = valid_pixels.size
                
                if total_poly_pixels == 0 or (valid_count / total_poly_pixels) < nodata_threshold:
                    results.append(None)
                    continue

                # 5. Calculate Statistics
                if method == 'mean':
                    val = self._calculate_mean(valid_pixels)
                    results.append(val)
                elif method == 'majority':
                    val = self._calculate_majority(valid_pixels)
                    results.append(val)

        # Assign results to GDF
        self.gdf[output_col] = results
        
        # Add Metadata
        self.gdf['meta_version'] = self.VERSION
        self.gdf['meta_processed_at'] = datetime.datetime.now().isoformat()
        self.gdf[f'meta_source_{raster_name}'] = raster_path

    def _calculate_mean(self, pixels):
        return float(np.mean(pixels))

    def _calculate_majority(self, pixels):
        """
        Calculates majority with Tie handling (returns 'Val1 / Val2').
        """
        values, counts = np.unique(pixels, return_counts=True)
        if len(counts) == 0:
            return None
            
        max_count = np.max(counts)
        candidates = values[counts == max_count]
        
        if len(candidates) == 1:
            return candidates[0].item()
        else:
            candidates.sort()
            return " / ".join([str(c) for c in candidates])

    def save_geojson(self, output_path: str):
        print(f"Saving to {output_path}...")
        self.gdf.to_file(output_path, driver='GeoJSON')
        print("Done.")


if __name__ == "__main__":
    
    # 1. Define Input/Output
    INPUT_GEOJSON = "PATH/TO/HEXSURFACE.geojson"
    OUTPUT_GEOJSON = "PATH/TO/OUTPUT.geojson"
    
    # 2. Configure Rasters List
    # Add as many dictionaries as you need. 
    # Keys: 'path', 'method' (mean/majority), 'threshold' (0.0 - 1.0)
    RASTER_CONFIG = [
        {
            "path": "path/to/dem_elevation.tif",
            "method": "mean",
            "threshold": 0.1
        },
        {
            "path": "path/to/slope_degrees.tif",
            "method": "mean",
            "threshold": 0.1
        },
        {
            "path": "path/to/land_use_land_cover.tif",
            "method": "majority",
            "threshold": 0.1
        }
    ]

    try:
        # Initialize the Enricher
        enricher = RasterBinEnricher(INPUT_GEOJSON)

        # Iterate through the config list
        print(f"Starting batch process for {len(RASTER_CONFIG)} rasters...")
        
        for config in RASTER_CONFIG:
            enricher.process_raster(
                raster_path=config["path"],
                method=config["method"],
                nodata_threshold=config.get("threshold", 0.1) # Default to 0.1 if missing
            )

        # Save Final Result
        enricher.save_geojson(OUTPUT_GEOJSON)
        
    except Exception as e:
        print(f"\nAn error occurred during processing: {e}")