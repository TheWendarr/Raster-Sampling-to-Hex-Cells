"""
sample_raster.py

Attribute polygon "bins" (e.g., H3 hex/pent cells or any polygon features) with
raster-derived statistics and write the result as a vector file.

Purpose
    This tool samples one or more rasters inside each input polygon and appends the
    computed value(s) as new attribute fields. It supports:
      - MEAN     (continuous rasters: DEM, slope, etc.)
      - MAJORITY (categorical rasters: landcover, soil class, etc.) with tie handling.

Inputs
    - --in-path (str)
        Path to input vector file containing polygon features (any format
        supported by geo_io: .geojson, .shp, .kml, .gdb).
    - --mean (str, repeatable)
        One or more raster paths to summarize with MEAN.
    - --majority (str, repeatable)
        One or more raster paths to summarize with MAJORITY.
    - --threshold (float, default 0.1)
        Minimum fraction of valid (non-NoData) pixels required within a polygon.
        If the valid fraction is below this value, the output for that polygon/raster
        is None.

Processing
    For each polygon and raster:
      1) Window the raster to the polygon bounds (fast subset read)
      2) Rasterize the polygon to a mask using a center-of-pixel rule (all_touched=False)
      3) Extract pixels under the mask and remove NoData values
      4) If valid pixel fraction >= threshold, compute MEAN or MAJORITY

Outputs
    - --out-path (str)
        Output vector file containing the original polygons plus:
          - <rasterName>_mean or <rasterName>_majority columns
          - meta_version, meta_processed_at, meta_source_<rasterName>

Notes
    - CRS mismatches between vector and raster will raise an error. Ensure the
      input polygons and rasters share a CRS, or reproject beforehand via
      geo_io.read_vector(..., epsg=<target>).
    - MAJORITY ties return the lowest tied value (deterministic). The count of
      tied classes is stored in a companion <rasterName>_majority_tie_count column.

Dependencies
    geopandas, rasterio, numpy, geo_io
"""

from __future__ import annotations

import datetime
import os
import argparse
import warnings
from pathlib import Path
from typing import Literal, Union

import numpy as np
import geopandas as gpd
import rasterio
from rasterio import features
from rasterio.windows import from_bounds

from geo_io import read_vector, write_vector

# Suppress specific rasterio warnings
warnings.filterwarnings("ignore", category=rasterio.errors.NotGeoreferencedWarning)


class RasterBinEnricher:
    """
    Enrich polygon bins (H3 cells or arbitrary polygons) with raster statistics.

    Handles 'mean' for continuous data and 'majority' for discrete/categorical data.

    Can be initialized from a file path (via geo_io) or directly from a GeoDataFrame
    for in-memory pipeline use.
    """

    VERSION = "0.2.0"

    def __init__(
        self,
        source: Union[str, Path, gpd.GeoDataFrame],
        layer: str | None = None,
        epsg: int | None = None,
    ):
        """
        Parameters
        ----------
        source : str, Path, or GeoDataFrame
            If a path, reads via geo_io.read_vector(). If a GeoDataFrame, uses
            it directly.
        layer : str, optional
            Layer name for multi-layer sources (FGDB, KML). Ignored if source
            is a GeoDataFrame.
        epsg : int, optional
            Reproject to this EPSG on read. Ignored if source is a GeoDataFrame.
        """
        if isinstance(source, gpd.GeoDataFrame):
            self.gdf = source.copy()
        else:
            path = Path(source)
            print(f"Loading vector data: {path}")
            self.gdf = read_vector(path, layer=layer, epsg=epsg)

        self.original_crs = self.gdf.crs

    def process_raster(
        self,
        raster_path: str,
        method: Literal["mean", "majority"],
        nodata_threshold: float = 0.1,
    ) -> None:
        """
        Enrich the internal GeoDataFrame with statistics from a raster.

        Parameters
        ----------
        raster_path : str
            Path to the raster file.
        method : 'mean' or 'majority'
            Aggregation method.
        nodata_threshold : float
            Minimum fraction of valid (non-NoData) pixels required within a
            polygon. Below this threshold the result is None.
        """
        if not os.path.exists(raster_path):
            print(f"Skipping: File not found -> {raster_path}")
            return

        raster_name = os.path.splitext(os.path.basename(raster_path))[0]
        output_col = f"{raster_name}_{method}"

        print(f"Processing raster: {raster_name} | Method: {method.upper()}")

        results = []
        tie_counts = [] if method == "majority" else None

        with rasterio.open(raster_path) as src:
            # CRS handling: if CRSs differ, temporarily reproject geometries
            # to the raster's CRS for accurate spatial sampling. The GeoDataFrame
            # itself stays in its original CRS.
            needs_reproject = False
            if self.original_crs is not None and src.crs is not None:
                if not self.original_crs.equals(src.crs):
                    needs_reproject = True
                    print(
                        f"  CRS mismatch (Vector: {self.original_crs}, "
                        f"Raster: {src.crs}) — reprojecting geometries for sampling"
                    )

            if needs_reproject:
                sample_geoms = self.gdf.geometry.to_crs(src.crs)
            else:
                sample_geoms = self.gdf.geometry

            nodata_val = src.nodata

            for idx, geom in enumerate(sample_geoms):

                if geom is None or geom.is_empty:
                    results.append(None)
                    if tie_counts is not None:
                        tie_counts.append(None)
                    continue

                # 1. Spatial windowing
                minx, miny, maxx, maxy = geom.bounds
                window = from_bounds(minx, miny, maxx, maxy, src.transform)

                try:
                    data = src.read(1, window=window)
                    win_transform = src.window_transform(window)
                except Exception:
                    results.append(None)
                    if tie_counts is not None:
                        tie_counts.append(None)
                    continue

                # 2. Rasterize geometry to mask
                mask_shape = data.shape
                if mask_shape[0] == 0 or mask_shape[1] == 0:
                    results.append(None)
                    if tie_counts is not None:
                        tie_counts.append(None)
                    continue

                mask_image = features.rasterize(
                    [(geom, 1)],
                    out_shape=mask_shape,
                    transform=win_transform,
                    fill=0,
                    dtype="uint8",
                    all_touched=False,  # Center-of-pixel rule
                )

                # 3. Extract valid data under mask
                pixels_in_poly = data[mask_image == 1]

                if pixels_in_poly.size == 0:
                    results.append(None)
                    if tie_counts is not None:
                        tie_counts.append(None)
                    continue

                valid_pixels = pixels_in_poly.copy()

                # Remove explicit NoData values (use np.isclose for float safety)
                if nodata_val is not None:
                    if np.issubdtype(valid_pixels.dtype, np.floating):
                        valid_pixels = valid_pixels[
                            ~np.isclose(valid_pixels, nodata_val)
                        ]
                    else:
                        valid_pixels = valid_pixels[valid_pixels != nodata_val]

                # Remove NaN values (float rasters may use NaN as effective nodata)
                if np.issubdtype(valid_pixels.dtype, np.floating):
                    valid_pixels = valid_pixels[~np.isnan(valid_pixels)]

                # 4. Check NoData threshold
                total_poly_pixels = pixels_in_poly.size
                valid_count = valid_pixels.size

                if total_poly_pixels == 0 or (valid_count / total_poly_pixels) < nodata_threshold:
                    results.append(None)
                    if tie_counts is not None:
                        tie_counts.append(None)
                    continue

                # 5. Calculate statistic
                if method == "mean":
                    results.append(self._calculate_mean(valid_pixels))
                elif method == "majority":
                    value, n_ties = self._calculate_majority(valid_pixels)
                    results.append(value)
                    tie_counts.append(n_ties)

        # Assign results to GDF
        self.gdf[output_col] = results
        if tie_counts is not None:
            self.gdf[f"{raster_name}_majority_tie_count"] = tie_counts

        # Add metadata
        self.gdf["meta_version"] = self.VERSION
        self.gdf["meta_processed_at"] = datetime.datetime.now().isoformat()
        self.gdf[f"meta_source_{raster_name}"] = raster_path

    @staticmethod
    def _calculate_mean(pixels: np.ndarray) -> float:
        return float(np.mean(pixels))

    @staticmethod
    def _calculate_majority(pixels: np.ndarray) -> tuple:
        """
        Calculate majority class with deterministic tie-breaking.

        Returns
        -------
        tuple of (value, tie_count)
            value: the lowest-valued class among those tied for first place.
            tie_count: number of classes sharing the max count (1 = no tie).
        """
        values, counts = np.unique(pixels, return_counts=True)
        if len(counts) == 0:
            return (None, None)

        max_count = np.max(counts)
        candidates = values[counts == max_count]
        candidates.sort()

        # Always return the lowest tied value (deterministic, consistent type)
        return (candidates[0].item(), int(len(candidates)))

    def save(
        self,
        output_path: str | Path,
        layer: str | None = None,
        epsg: int | None = None,
        overwrite: bool = False,
    ) -> Path:
        """
        Write the enriched GeoDataFrame via geo_io.

        Parameters
        ----------
        output_path : str or Path
            Output file path (format inferred from extension).
        layer : str, optional
            Layer name for multi-layer outputs (FGDB).
        epsg : int, optional
            Reproject before writing.
        overwrite : bool
            Overwrite existing output.

        Returns
        -------
        Path
            Resolved output path.
        """
        print(f"Saving to {output_path}...")
        result = write_vector(
            self.gdf, output_path, layer=layer, epsg=epsg, overwrite=overwrite
        )
        print(f"Done. Wrote {len(self.gdf)} features.")
        return result

    @property
    def result(self) -> gpd.GeoDataFrame:
        """Access the enriched GeoDataFrame directly for in-memory pipelines."""
        return self.gdf


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Sample rasters into polygon bins and save enriched vector output. "
            "Supports GeoJSON, Shapefile, KML, and FGDB via geo_io."
        )
    )
    parser.add_argument(
        "--in-path",
        required=True,
        help="Input vector file path (.geojson, .shp, .kml, .gdb).",
    )
    parser.add_argument(
        "--out-path",
        required=True,
        help="Output vector file path (format inferred from extension).",
    )
    parser.add_argument(
        "--layer",
        default=None,
        help="Layer name for multi-layer inputs (FGDB, KML).",
    )
    parser.add_argument(
        "--out-layer",
        default=None,
        help="Output layer name for FGDB writes.",
    )
    parser.add_argument(
        "--mean",
        action="append",
        default=[],
        help="Raster path to process with MEAN (repeatable).",
    )
    parser.add_argument(
        "--majority",
        action="append",
        default=[],
        help="Raster path to process with MAJORITY (repeatable).",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.1,
        help="Minimum valid-pixel ratio per polygon (default: 0.1).",
    )
    parser.add_argument(
        "--epsg",
        type=int,
        default=None,
        help="Reproject input to this EPSG before sampling.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing output file.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if not args.mean and not args.majority:
        parser.error("Provide at least one raster via --mean or --majority.")

    if not 0.0 <= args.threshold <= 1.0:
        parser.error("--threshold must be between 0.0 and 1.0.")

    raster_config = [
        {"path": path, "method": "mean", "threshold": args.threshold}
        for path in args.mean
    ] + [
        {"path": path, "method": "majority", "threshold": args.threshold}
        for path in args.majority
    ]

    try:
        enricher = RasterBinEnricher(
            args.in_path, layer=args.layer, epsg=args.epsg
        )
        print(f"Starting batch process for {len(raster_config)} raster(s)...")

        for config in raster_config:
            enricher.process_raster(
                raster_path=config["path"],
                method=config["method"],
                nodata_threshold=config["threshold"],
            )

        enricher.save(
            args.out_path,
            layer=args.out_layer,
            overwrite=args.overwrite,
        )
        return 0

    except Exception as exc:
        print(f"\n[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    import sys

    raise SystemExit(main())