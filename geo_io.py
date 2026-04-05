"""
geo_io.py

Reads and writes vector geospatial data (GeoJSON, Shapefile, KML, FGDB)
via a unified interface backed by pyogrio (with Arrow acceleration) and
a GeoPandas/Fiona fallback.

Usage (library):
    from geo_io import read_vector, write_vector, list_layers

    gdf = read_vector("extent.shp")
    gdf = read_vector("data.gdb", layer="roads", epsg=4326)

    write_vector(gdf, "output.geojson")
    write_vector(gdf, "output.gdb", layer="roads")

    layers = list_layers("data.gdb")

Usage (CLI):
    python geo_io.py --in extent.shp --out extent.geojson
    python geo_io.py --in data.gdb --layer roads --out roads.shp --epsg 4326
    python geo_io.py --in data.gdb --list-layers

Requires:
    geopandas, pyogrio (recommended), shapely

Optional:
    polars (for to_polars() convenience function)
"""


import argparse
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import geopandas as gpd

# Optional dependency: pyogrio (preferred I/O engine)
try:
    import pyogrio  # type: ignore

    _HAS_PYOGRIO = True
except ImportError:
    pyogrio = None  # type: ignore
    _HAS_PYOGRIO = False


# Constants & Format Registry

# Maps normalized format keys -> (GDAL driver name, file extension)
_FORMAT_REGISTRY: Dict[str, Tuple[str, str]] = {
    "geojson": ("GeoJSON", ".geojson"),
    "shp": ("ESRI Shapefile", ".shp"),
    "kml": ("KML", ".kml"),
    "fgdb": ("OpenFileGDB", ".gdb"),
}

# Extension -> normalized format key
_EXT_TO_FORMAT: Dict[str, str] = {
    ".geojson": "geojson",
    ".json": "geojson",
    ".shp": "shp",
    ".kml": "kml",
    ".gdb": "fgdb",
}

SUPPORTED_EXTENSIONS = set(_EXT_TO_FORMAT.keys())


# Internal Helpers

def _detect_format(path: Path) -> str:
    """Detect the normalized format key from a file path's extension."""
    ext = path.suffix.lower()

    # FGDB: path ends in .gdb and is (or will be) a directory
    if ext == ".gdb":
        return "fgdb"

    fmt = _EXT_TO_FORMAT.get(ext)
    if fmt is None:
        raise ValueError(
            f"Unsupported file extension: '{ext}'. "
            f"Supported: {sorted(SUPPORTED_EXTENSIONS)}"
        )
    return fmt


def _resolve_driver(fmt: str) -> str:
    """Return the GDAL driver name for a normalized format key."""
    if fmt == "kml":
        return _pick_kml_driver()
    entry = _FORMAT_REGISTRY.get(fmt)
    if entry is None:
        raise ValueError(f"Unknown format key: '{fmt}'")
    return entry[0]


def _pick_kml_driver() -> str:
    """
    GDAL exposes KML write support under 'KML' or 'LIBKML' depending on build.
    Probe for a writable driver and return the first match.
    """
    if not _HAS_PYOGRIO:
        return "KML"  # GeoPandas/Fiona will handle driver resolution

    drivers = pyogrio.list_drivers(write=True)
    for candidate in ("KML", "LIBKML"):
        mode = drivers.get(candidate)
        if mode is not None and "w" in str(mode):
            return candidate
    raise RuntimeError(
        "KML write is not supported by your GDAL/pyogrio build. "
        "No writable KML or LIBKML driver found."
    )


def _check_fgdb_write_support() -> None:
    """Raise early if OpenFileGDB write isn't available (requires GDAL >= 3.6)."""
    if not _HAS_PYOGRIO:
        raise RuntimeError(
            "FGDB write requires pyogrio with GDAL >= 3.6. "
            "Install with: pip install pyogrio"
        )
    drivers = pyogrio.list_drivers(write=True)
    mode = drivers.get("OpenFileGDB")
    if mode is None or "w" not in str(mode):
        gdal_ver = getattr(pyogrio, "__gdal_version__", "unknown")
        raise RuntimeError(
            f"FGDB write not supported (GDAL {gdal_ver}). "
            "OpenFileGDB write requires GDAL >= 3.6. "
            "A conda-forge environment is the most reliable way to get this."
        )


def _sanitize_layer_name(name: str, max_len: int = 120) -> str:
    """
    Produce a filesystem/GDB-safe layer name:
    letters, digits, underscores only; must start with a letter.
    """
    n = (name or "").strip()
    n = re.sub(r"\s+", "_", n)
    n = re.sub(r"[^A-Za-z0-9_]+", "_", n)
    n = re.sub(r"_+", "_", n).strip("_")
    if not n:
        n = "layer"
    if not re.match(r"^[A-Za-z]", n):
        n = f"lyr_{n}"
    return n[:max_len]


def _parse_layer_list(raw: Any) -> List[Dict[str, Any]]:
    """
    Normalize pyogrio.list_layers() output (which can be shaped (n,2) or (2,n))
    into a consistent list of dicts: [{"name": str, "geometry_type": str|None}, ...]
    """
    if raw is None:
        return []

    rows = raw.tolist() if hasattr(raw, "tolist") else list(raw)
    if not rows:
        return []

    # Shape (n, 2): each row is [name, geom_type]
    if isinstance(rows[0], (list, tuple)) and len(rows[0]) == 2:
        pass  # already in the right shape
    # Shape (2, n): needs transpose
    elif (
        len(rows) == 2
        and isinstance(rows[0], list)
        and isinstance(rows[1], list)
        and len(rows[0]) == len(rows[1])
    ):
        rows = list(zip(rows[0], rows[1]))
    # Flat list of names
    else:
        rows = [(x, None) for x in rows]

    return [
        {
            "name": None if name is None else str(name),
            "geometry_type": gtype,
        }
        for name, gtype in rows
    ]


# Public API: list_layers

def list_layers(
    path: Union[str, Path],
    spatial_only: bool = False,
) -> List[Dict[str, Any]]:
    """
    List layers/feature classes in a multi-layer source (FGDB, KML, etc.).

    Parameters
    ----------
    path : str or Path
        Path to the geospatial data source.
    spatial_only : bool
        If True, exclude non-spatial tables (geometry_type is None/empty).

    Returns
    -------
    list of dict
        Each dict has keys 'name' (str) and 'geometry_type' (str or None).
    """
    p = Path(path).expanduser().resolve()
    if not p.exists():
        raise FileNotFoundError(f"Path not found: {p}")

    if _HAS_PYOGRIO:
        raw = pyogrio.list_layers(str(p))
        layers = _parse_layer_list(raw)
    else:
        import fiona

        names = fiona.listlayers(str(p))
        layers = [{"name": n, "geometry_type": None} for n in names]

    if spatial_only:
        layers = [
            rec
            for rec in layers
            if rec.get("geometry_type") not in (None, "", "None")
        ]

    return layers


# Public API: read_vector

def read_vector(
    path: Union[str, Path],
    layer: Optional[str] = None,
    epsg: Optional[int] = None,
    force_2d: bool = False,
    use_arrow: bool = True,
) -> gpd.GeoDataFrame:
    """
    Read a vector geospatial file into a GeoDataFrame.

    Supports GeoJSON (.geojson/.json), Shapefile (.shp), KML (.kml),
    and Esri File Geodatabase (.gdb).

    Parameters
    ----------
    path : str or Path
        Path to the input file or .gdb directory.
    layer : str, optional
        Layer/feature class name. Required for multi-layer sources (FGDB, KML)
        unless the source has exactly one layer.
    epsg : int, optional
        If provided, reproject the output to this EPSG code.
    force_2d : bool
        If True, drop Z coordinates.
    use_arrow : bool
        If True (default), attempt pyogrio Arrow acceleration for faster reads.

    Returns
    -------
    geopandas.GeoDataFrame

    Raises
    ------
    FileNotFoundError
        If the input path does not exist.
    ValueError
        If a multi-layer source is opened without specifying a layer.
    RuntimeError
        If required dependencies are missing.
    """
    p = Path(path).expanduser().resolve()
    if not p.exists():
        raise FileNotFoundError(f"Input not found: {p}")

    fmt = _detect_format(p)

    # For multi-layer sources, resolve the layer name
    resolved_layer = _resolve_layer(p, fmt, layer)

    # Read via pyogrio (fast path) or GeoPandas/Fiona (fallback)
    gdf = _read_with_engine(p, resolved_layer, use_arrow, force_2d)

    # Reproject if requested
    if epsg is not None:
        if gdf.crs is None:
            raise ValueError(
                f"Cannot reproject to EPSG:{epsg} — input has no CRS defined. "
                "Set a CRS first with gdf.set_crs()."
            )
        gdf = gdf.to_crs(epsg=epsg)

    return gdf


def _resolve_layer(
    path: Path, fmt: str, layer: Optional[str]
) -> Optional[str]:
    """
    For multi-layer formats (FGDB, KML), validate or auto-select the layer.
    For single-layer formats (GeoJSON, Shapefile), return None.
    """
    if fmt in ("geojson", "shp"):
        # Single-layer formats: ignore the layer argument
        return None

    # Multi-layer format: FGDB or KML
    available = list_layers(path, spatial_only=True)

    if not available:
        raise ValueError(f"No spatial layers found in: {path}")

    available_names = [rec["name"] for rec in available if rec["name"]]

    if layer is not None:
        if layer not in available_names:
            raise ValueError(
                f"Layer '{layer}' not found in {path.name}. "
                f"Available: {available_names}"
            )
        return layer

    # No layer specified: auto-select if there's exactly one
    if len(available_names) == 1:
        return available_names[0]

    raise ValueError(
        f"Multiple layers in {path.name} — specify one with layer=. "
        f"Available: {available_names}"
    )


def _read_with_engine(
    path: Path,
    layer: Optional[str],
    use_arrow: bool,
    force_2d: bool,
) -> gpd.GeoDataFrame:
    """Attempt pyogrio read with Arrow, fall back to GeoPandas/Fiona."""
    kwargs: Dict[str, Any] = {}
    if layer is not None:
        kwargs["layer"] = layer
    if force_2d:
        kwargs["force_2d"] = True

    if _HAS_PYOGRIO and use_arrow:
        try:
            return pyogrio.read_dataframe(str(path), use_arrow=True, **kwargs)
        except TypeError:
            # Older pyogrio without use_arrow param
            return pyogrio.read_dataframe(str(path), **kwargs)
        except Exception:
            # Arrow acceleration failed; try without
            try:
                return pyogrio.read_dataframe(str(path), **kwargs)
            except Exception:
                pass  # Fall through to GeoPandas

    if _HAS_PYOGRIO:
        try:
            return pyogrio.read_dataframe(str(path), **kwargs)
        except Exception:
            pass  # Fall through to GeoPandas

    # GeoPandas/Fiona fallback
    read_kwargs: Dict[str, Any] = {}
    if layer is not None:
        read_kwargs["layer"] = layer
    return gpd.read_file(str(path), **read_kwargs)

# Public API: write_vector

def write_vector(
    gdf: gpd.GeoDataFrame,
    path: Union[str, Path],
    layer: Optional[str] = None,
    epsg: Optional[int] = None,
    overwrite: bool = False,
    use_arrow: bool = True,
) -> Path:
    """
    Write a GeoDataFrame to a vector geospatial file.

    Supports GeoJSON (.geojson/.json), Shapefile (.shp), KML (.kml),
    and Esri File Geodatabase (.gdb).

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        The data to write.
    path : str or Path
        Output file path. Format is inferred from the extension.
    layer : str, optional
        Layer name for multi-layer formats (FGDB). Ignored for single-layer
        formats. Auto-generated from the filename stem if not provided.
    epsg : int, optional
        If provided, reproject the data before writing.
    overwrite : bool
        If True, overwrite an existing file. If False and the file exists,
        raises FileExistsError.
    use_arrow : bool
        If True (default), attempt pyogrio Arrow acceleration for faster writes.

    Returns
    -------
    Path
        The resolved output path.

    Raises
    ------
    FileExistsError
        If the output exists and overwrite is False.
    ValueError
        If the GeoDataFrame is empty or has no geometry column.
    RuntimeError
        If required drivers are unavailable.
    """
    if gdf is None or not isinstance(gdf, gpd.GeoDataFrame):
        raise TypeError("Expected a GeoDataFrame, got: " + type(gdf).__name__)

    if "geometry" not in gdf.columns or gdf.geometry is None:
        raise ValueError("GeoDataFrame has no geometry column.")

    out = Path(path).expanduser().resolve()
    fmt = _detect_format(out)
    driver = _resolve_driver(fmt)

    # Reproject if requested (KML forces 4326)
    write_gdf = gdf
    effective_epsg = epsg
    if fmt == "kml" and epsg is None:
        effective_epsg = 4326

    if effective_epsg is not None:
        if write_gdf.crs is None:
            raise ValueError(
                f"Cannot reproject to EPSG:{effective_epsg} — GeoDataFrame has no CRS. "
                "Set a CRS first with gdf.set_crs()."
            )
        write_gdf = write_gdf.to_crs(epsg=effective_epsg)

    # Handle overwrite / existence checks
    if fmt == "fgdb":
        _check_fgdb_write_support()
        resolved_layer = _sanitize_layer_name(layer or out.stem)
    else:
        resolved_layer = None
        if out.exists() and not overwrite:
            raise FileExistsError(
                f"Output already exists: {out}. Use overwrite=True to replace."
            )

    # Ensure parent directory exists
    out.parent.mkdir(parents=True, exist_ok=True)

    # Write via pyogrio (fast path) or GeoPandas/Fiona (fallback)
    _write_with_engine(
        write_gdf, out, driver, resolved_layer, use_arrow, fmt
    )

    return out


def _write_with_engine(
    gdf: gpd.GeoDataFrame,
    path: Path,
    driver: str,
    layer: Optional[str],
    use_arrow: bool,
    fmt: str,
) -> None:
    """Attempt pyogrio write with Arrow, fall back to GeoPandas/Fiona."""
    write_kwargs: Dict[str, Any] = {"driver": driver}
    if layer is not None:
        write_kwargs["layer"] = layer
    if fmt == "fgdb":
        write_kwargs["append"] = False

    if _HAS_PYOGRIO and use_arrow:
        try:
            pyogrio.write_dataframe(gdf, str(path), use_arrow=True, **write_kwargs)
            return
        except TypeError:
            try:
                pyogrio.write_dataframe(gdf, str(path), **write_kwargs)
                return
            except Exception:
                pass
        except Exception:
            pass

    if _HAS_PYOGRIO:
        try:
            pyogrio.write_dataframe(gdf, str(path), **write_kwargs)
            return
        except Exception:
            pass

    # GeoPandas/Fiona fallback
    fallback_kwargs: Dict[str, Any] = {"driver": driver}
    if layer is not None:
        fallback_kwargs["layer"] = layer
    gdf.to_file(str(path), **fallback_kwargs)


# Public API: Polars convenience

def to_polars(
    gdf: gpd.GeoDataFrame,
    geometry_format: str = "wkt",
) -> Any:
    """
    Convert a GeoDataFrame to a Polars DataFrame.

    Geometry is serialized to a string column since Polars has no native
    geometry type. This is useful for fast attribute-only operations
    (groupby, join, filter) before converting back to GeoPandas.

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        Input data.
    geometry_format : str
        How to serialize geometry: 'wkt' (default) or 'wkb_hex'.

    Returns
    -------
    polars.DataFrame

    Raises
    ------
    ImportError
        If polars is not installed.
    """
    try:
        import polars as pl
    except ImportError:
        raise ImportError(
            "Polars is not installed. Install with: pip install polars"
        )

    import pandas as pd

    # Serialize geometry to a string column
    df = pd.DataFrame(gdf.drop(columns="geometry"))

    if geometry_format == "wkt":
        df["geometry_wkt"] = gdf.geometry.to_wkt()
    elif geometry_format == "wkb_hex":
        df["geometry_wkb_hex"] = gdf.geometry.to_wkb(hex=True)
    else:
        raise ValueError(f"geometry_format must be 'wkt' or 'wkb_hex', got '{geometry_format}'")

    # Preserve CRS as metadata in a constant column (recoverable on round-trip)
    crs_str = str(gdf.crs) if gdf.crs is not None else ""
    df["_crs"] = crs_str

    return pl.from_pandas(df)


def from_polars(
    pldf: Any,
    geometry_column: str = "geometry_wkt",
    crs_column: str = "_crs",
) -> gpd.GeoDataFrame:
    """
    Convert a Polars DataFrame back to a GeoDataFrame.

    Parameters
    ----------
    pldf : polars.DataFrame
        Input Polars DataFrame with a serialized geometry column.
    geometry_column : str
        Name of the column containing WKT or WKB hex geometry strings.
    crs_column : str
        Name of the column containing the CRS string. Dropped after use.

    Returns
    -------
    geopandas.GeoDataFrame
    """
    try:
        import polars as pl
    except ImportError:
        raise ImportError("Polars is not installed. Install with: pip install polars")

    from shapely import wkt, wkb

    pdf = pldf.to_pandas()

    # Recover CRS
    crs = None
    if crs_column in pdf.columns:
        crs_val = pdf[crs_column].iloc[0] if len(pdf) > 0 else ""
        if crs_val:
            crs = crs_val
        pdf = pdf.drop(columns=[crs_column])

    # Deserialize geometry
    if geometry_column not in pdf.columns:
        raise ValueError(f"Geometry column '{geometry_column}' not found in DataFrame.")

    geom_series = pdf[geometry_column]
    pdf = pdf.drop(columns=[geometry_column])

    # Detect WKT vs WKB hex
    if "wkb" in geometry_column.lower():
        geometries = gpd.GeoSeries.from_wkb(geom_series, crs=crs)
    else:
        geometries = gpd.GeoSeries.from_wkt(geom_series, crs=crs)

    return gpd.GeoDataFrame(pdf, geometry=geometries, crs=crs)


# CLI Interface

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="geo_io",
        description=(
            "Convert between geospatial vector formats: "
            "GeoJSON, Shapefile, KML, Esri File Geodatabase."
        ),
    )
    p.add_argument(
        "--in",
        dest="in_path",
        required=True,
        help="Input file path (.geojson, .shp, .kml, .gdb).",
    )
    p.add_argument(
        "--out",
        dest="out_path",
        default=None,
        help="Output file path. Format is inferred from extension.",
    )
    p.add_argument(
        "--layer",
        default=None,
        help="Layer/feature class name (required for multi-layer sources with >1 layer).",
    )
    p.add_argument(
        "--out-layer",
        default=None,
        help="Output layer name for FGDB writes. Defaults to output filename stem.",
    )
    p.add_argument(
        "--epsg",
        type=int,
        default=None,
        help="Reproject output to this EPSG code.",
    )
    p.add_argument(
        "--force-2d",
        action="store_true",
        help="Drop Z coordinates on read.",
    )
    p.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing output files.",
    )
    p.add_argument(
        "--list-layers",
        action="store_true",
        help="List layers in the input source and exit.",
    )
    p.add_argument(
        "--no-arrow",
        action="store_true",
        help="Disable pyogrio Arrow acceleration.",
    )
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    in_path = Path(args.in_path)
    use_arrow = not args.no_arrow

    # --list-layers mode
    if args.list_layers:
        try:
            layers = list_layers(in_path)
        except Exception as e:
            print(f"[ERROR] {e}", file=sys.stderr)
            return 1

        if not layers:
            print("(no layers found)")
            return 0

        print(f"{'Layer Name':<40} {'Geometry Type'}")
        print("-" * 60)
        for rec in layers:
            name = rec.get("name") or "(unnamed)"
            gtype = rec.get("geometry_type") or "(non-spatial)"
            print(f"{name:<40} {gtype}")
        return 0

    # Convert mode: --out is required
    if args.out_path is None:
        print("[ERROR] --out is required for conversion.", file=sys.stderr)
        return 2

    out_path = Path(args.out_path)

    try:
        print(f"Reading: {in_path}" + (f" [layer={args.layer}]" if args.layer else ""))
        gdf = read_vector(
            in_path,
            layer=args.layer,
            force_2d=args.force_2d,
            use_arrow=use_arrow,
        )
        print(f"  {len(gdf)} features, CRS={gdf.crs}")

        result = write_vector(
            gdf,
            out_path,
            layer=args.out_layer,
            epsg=args.epsg,
            overwrite=args.overwrite,
            use_arrow=use_arrow,
        )
        print(f"Wrote: {result}")
        return 0

    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())