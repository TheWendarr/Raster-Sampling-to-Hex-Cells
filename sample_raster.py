# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///

import os
import datetime
import argparse
import warnings
from typing import Literal

import numpy as np
import geopandas as gpd
import rasterio
from rasterio import features
from rasterio.windows import from_bounds

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


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Sample rasters into hex/polygon bins and save enriched GeoJSON."
    )
    parser.add_argument(
        "--in-path",
        required=True,
        help="Input hex/polygon GeoJSON path.",
    )
    parser.add_argument(
        "--out-path",
        required=True,
        help="Output GeoJSON path.",
    )
    parser.add_argument(
        "--mean",
        action="append",
        default=[],
        help="Raster path to process with MEAN. Use multiple times for multiple rasters.",
    )
    parser.add_argument(
        "--majority",
        action="append",
        default=[],
        help="Raster path to process with MAJORITY. Use multiple times for multiple rasters.",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.1,
        help="Minimum valid-pixel ratio per polygon (default: 0.1).",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if not args.mean and not args.majority:
        parser.error("Provide at least one raster via --mean or --majority.")

    if not 0.0 <= args.threshold <= 1.0:
        parser.error("--threshold must be between 0.0 and 1.0.")

    raster_config = (
        [{"path": path, "method": "mean", "threshold": args.threshold} for path in args.mean]
        + [{"path": path, "method": "majority", "threshold": args.threshold} for path in args.majority]
    )

    try:
        enricher = RasterBinEnricher(args.in_path)
        print(f"Starting batch process for {len(raster_config)} rasters...")

        for config in raster_config:
            enricher.process_raster(
                raster_path=config["path"],
                method=config["method"],
                nodata_threshold=config["threshold"],
            )

        enricher.save_geojson(args.out_path)
        return 0
    except Exception as exc:
        print(f"\nAn error occurred during processing: {exc}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
