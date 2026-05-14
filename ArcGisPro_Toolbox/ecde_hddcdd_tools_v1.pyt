# -*- coding: utf-8 -*-
"""
ECDE HDD/CDD Tools
==================

ArcGIS Pro Python toolbox for working with Copernicus ECDE
(sis-ecde-climate-indicators) Heating and Cooling Degree Days data.

Tools provided
--------------
1. Convert ECDE NetCDF to GeoTIFF
2. Inspect ECDE NetCDF (debug)
3. Compute Ensemble Statistics (projections only)

Changelog
---------
v0.4.1 (2026-05-14)
  - Fix: PolygonToRaster was being called with value_field="OID@", which
    is a SearchCursor/CalculateField token, not a real field name.
    arcpy.conversion.PolygonToRaster requires the actual OID column name
    ("FID" for shapefiles, "OBJECTID" for File GDB feature classes).
    Now resolved dynamically via arcpy.Describe(clip_fc).OIDFieldName.

v0.4 (2026-05-14)
  - Tool 1 conversion engine rewritten to use netCDF4 + osgeo.gdal
    directly, bypassing arcpy.md.MakeNetCDFRasterLayer + CopyRaster.
    Expected ~20-40x speedup. Both libraries ship with the default
    ArcGIS Pro 3.x Python environment.
  - Smart spatial slicing: when a clip polygon is provided, the engine
    (a) computes the polygon's bounding box, (b) reads ONLY the bbox
    region from the NetCDF (instead of the full European domain),
    (c) rasterises the polygon to the NetCDF grid ONCE per unique
    grid signature (cached across files), (d) applies the mask
    in-memory. ~99% reduction in data read when clipping to Portugal.
    Without clip, reads full extent (no mask, no slicing).
  - Per-file timing: log line now shows elapsed seconds and seconds
    per slice, e.g. "wrote 148 GeoTIFFs in 12.3s (0.08s/slice)".
  - Fixed latent bug: removed PROJECTIONS_DEFAULT_START_YEAR=2006
    fallback that would have produced wrong years if time values
    failed to parse (the projection files actually start in 1950 or
    1951 or 1970, not 2006). Time parsing now uses netCDF4.num2date
    which handles all CF calendars including 360-day (HadGEM2-ES).
  - Output GeoTIFFs now use DEFLATE compression natively (via GDAL),
    typically ~70% smaller than uncompressed.

v0.3 (2026-05-14): SetProgressor, verbose flag, arcpy.mp layer resolver,
  Tool 3 projections-only clarity, period validation.
v0.2 (2026-05-14): Time parsing fix via filename, GPRasterLayer param.
v0.1 (2026-05-14): initial version.

Filename convention parsed:
  Reanalysis:  03_heating_degree_days-reanalysis-yearly-grid-1940-2025-v2.0.nc
  Projections: 03_heating_degree_days-projections-yearly-rcp_4_5-<rcm>-<gcm>-<ens>-grid-v2.0.nc

Notes on the projection files
-----------------------------
Each *-projections-*.nc actually contains 1950-2100 (or 1951-2100 / 1970-2100
depending on the chain), with the historical experiment (1950/51/70-2005) and
the RCP scenario (2006-2100) merged and bias-adjusted with the same correction
function. This means the historical baseline is embedded; no separate
historical experiment download is needed for delta calculations within the
same chain.

Multiprocessing
---------------
NOT implemented. With the netCDF4+GDAL engine, multiprocessing is technically
feasible (workers don't touch arcpy), but the expected gain (~3-5 min -> ~1-2 min
for full 32-file run) does not justify the complexity of process spawning,
pickling, and cross-process progress synchronisation for a workflow that runs
once per data refresh.

Author: drafted for Pedro Goncalves, 41 Norte / LNEG, 2026.
"""

from __future__ import annotations

import csv
import datetime as dt
import os
import re
import time as _time

import arcpy
import numpy as np

# Required for the fast conversion path. Both ship with ArcGIS Pro 3.x.
try:
    import netCDF4
except ImportError as e:
    raise ImportError(
        "netCDF4 is required. It ships with ArcGIS Pro 3.x by default. "
        "If missing, run: conda install -c esri netCDF4 in your ArcGIS Pro env."
    ) from e

try:
    from osgeo import gdal, osr
    gdal.UseExceptions()
except ImportError as e:
    raise ImportError(
        "osgeo (GDAL) is required. It ships with ArcGIS Pro 3.x by default."
    ) from e


# ===========================================================================
# Configuration
# ===========================================================================

RCM_CANONICAL = {
    "cclm4_8_17": "CCLM4-8-17",
    "hirham5":    "HIRHAM5",
    "racmo22e":   "RACMO22E",
    "rca4":       "RCA4",
    "wrf381p":    "WRF381P",
}

GCM_CANONICAL = {
    "mpi_esm_lr":   "MPI-ESM-LR",
    "noresm1_m":    "NorESM1-M",
    "ec_earth":     "EC-EARTH",
    "hadgem2_es":   "HadGEM2-ES",
    "ipsl_cm5a_mr": "IPSL-CM5A-MR",
}

INDICATOR_TO_VAR = {
    "03_heating_degree_days": "HDD",
    "04_cooling_degree_days": "CDD",
}

NODATA_VAL = -9999.0

# GeoTIFF creation options for the fast (GDAL) writer
GTIFF_OPTIONS = [
    "COMPRESS=DEFLATE",
    "PREDICTOR=2",
    "TILED=YES",
    "BLOCKXSIZE=256",
    "BLOCKYSIZE=256",
    "BIGTIFF=IF_SAFER",
]


# ===========================================================================
# Filename parsing
# ===========================================================================

def parse_nc_filename(nc_path):
    """Parse an ECDE NetCDF filename into a metadata dict."""
    stem = os.path.splitext(os.path.basename(nc_path))[0]
    meta = {
        "indicator": "unknown",
        "variable": "VAR",
        "origin": "unknown",
        "agg": "unknown",
        "experiment": "",
        "rcm": "",
        "gcm": "",
        "ensemble": "",
        "version": "",
    }

    for ind_key, var_label in INDICATOR_TO_VAR.items():
        if stem.startswith(ind_key):
            meta["indicator"] = ind_key
            meta["variable"] = var_label
            break

    parts = stem.split("-")
    if len(parts) >= 2:
        meta["origin"] = parts[1]
    if len(parts) >= 3:
        meta["agg"] = parts[2]

    m = re.search(r"rcp_(\d)_(\d)", stem)
    if m:
        meta["experiment"] = "rcp" + m.group(1) + m.group(2)

    for token, canonical in RCM_CANONICAL.items():
        if token in stem:
            meta["rcm"] = canonical
            break
    for token, canonical in GCM_CANONICAL.items():
        if token in stem:
            meta["gcm"] = canonical
            break

    m = re.search(r"r\d+i\d+p\d+", stem)
    if m:
        meta["ensemble"] = m.group(0)
    m = re.search(r"v(\d+\.\d+)", stem)
    if m:
        meta["version"] = "v" + m.group(1)

    return meta


def infer_start_year_from_filename(nc_path):
    """Determine the time series start year from the ECDE filename.

    Used only as fallback when netCDF4.num2date fails. For reanalysis,
    the filename has explicit YYYY-YYYY range. For projections, returns
    None (the time dim values are the authoritative source).
    """
    name = os.path.basename(nc_path)
    m = re.search(r"-(\d{4})-(\d{4})-", name)
    if m:
        return int(m.group(1))
    return None


def build_output_name(meta, year, period_label=None):
    """Construct output GeoTIFF filename from metadata + year."""
    var = meta["variable"]
    suffix = "_{}".format(period_label) if period_label else ""
    if meta["origin"] == "reanalysis":
        return "{var}_reanalysis_{agg}_{year}{suffix}.tif".format(
            var=var, agg=meta["agg"], year=year, suffix=suffix
        )
    if meta["origin"] == "projections":
        return "{var}_{exp}_{rcm}_{gcm}_{ens}_{year}{suffix}.tif".format(
            var=var,
            exp=meta["experiment"] or "rcpXX",
            rcm=meta["rcm"] or "RCM",
            gcm=meta["gcm"] or "GCM",
            ens=meta["ensemble"] or "rXi1p1",
            year=year, suffix=suffix,
        )
    return "{var}_{stem}_{year}{suffix}.tif".format(
        var=var, stem=meta["indicator"], year=year, suffix=suffix
    )


def parse_output_filename(tif_path):
    """Parse a Tool 1 output filename back into metadata."""
    stem = os.path.splitext(os.path.basename(tif_path))[0]

    m = re.match(
        r"^(HDD|CDD)_(rcp\d{2})_([\w\-\.]+?)_([\w\-\.]+?)_(r\d+i\d+p\d+)_(\d{4})$",
        stem,
    )
    if m:
        return {
            "variable": m.group(1),
            "experiment": m.group(2),
            "rcm": m.group(3),
            "gcm": m.group(4),
            "ensemble": m.group(5),
            "year": int(m.group(6)),
            "origin": "projections",
        }

    m = re.match(r"^(HDD|CDD)_reanalysis_(\w+)_(\d{4})$", stem)
    if m:
        return {
            "variable": m.group(1),
            "experiment": "reanalysis",
            "agg": m.group(2),
            "year": int(m.group(3)),
            "origin": "reanalysis",
        }
    return None


# ===========================================================================
# NetCDF structure inspection (Tool 2 still uses arcpy.NetCDFFileProperties)
# ===========================================================================

def discover_structure(nc_props, user_variable=None):
    """Return (variable_name, x_dim, y_dim, t_dim) for an ECDE NetCDF."""
    dims = list(nc_props.getDimensions())
    dims_lower = {d.lower(): d for d in dims}

    x_dim = (dims_lower.get("lon") or dims_lower.get("longitude")
             or dims_lower.get("x"))
    y_dim = (dims_lower.get("lat") or dims_lower.get("latitude")
             or dims_lower.get("y"))
    t_dim = (dims_lower.get("time") or dims_lower.get("t")
             or dims_lower.get("stdtime"))

    missing = [n for n, v in (("x", x_dim), ("y", y_dim), ("t", t_dim)) if not v]
    if missing:
        raise ValueError(
            "Missing expected dimensions {} in NetCDF. Dims found: {}".format(
                missing, dims
            )
        )

    if user_variable:
        return user_variable, x_dim, y_dim, t_dim

    all_vars = list(nc_props.getVariables())
    coord_like = set(dims) | {"time_bnds", "lat_bnds", "lon_bnds",
                              "longitude_bnds", "latitude_bnds",
                              "crs", "spatial_ref"}
    data_vars = [v for v in all_vars if v not in coord_like]
    if not data_vars:
        raise ValueError("No data variable found in NetCDF.")

    priorities = [v for v in data_vars if v.lower() in ("hdd", "cdd")]
    if not priorities:
        priorities = [v for v in data_vars if "degree" in v.lower()]
    return (priorities[0] if priorities else data_vars[0]), x_dim, y_dim, t_dim


def extract_year_from_value(t_value):
    """Try to extract a year from various time representations (for Tool 2)."""
    if t_value is None:
        return None
    if not isinstance(t_value, str):
        try:
            y = int(getattr(t_value, "year", None))
            if 1800 < y < 2200:
                return y
        except (TypeError, ValueError):
            pass
    if isinstance(t_value, str):
        s = t_value.strip()
        m = re.match(r"^(\d{4})[/\-T\s]", s)
        if m:
            return int(m.group(1))
        m = re.match(r"^\d{1,2}[/\-]\d{1,2}[/\-](\d{4})", s)
        if m:
            return int(m.group(1))
        m = re.match(r"^(\d{4})$", s)
        if m:
            return int(m.group(1))
        return None
    try:
        v = int(t_value)
        if 1800 < v < 2200:
            return v
    except (TypeError, ValueError):
        pass
    return None


# ===========================================================================
# Layer resolution
# ===========================================================================

def get_layer_source_nc(layer_name):
    """Resolve a map layer name to its underlying .nc file path via arcpy.mp."""
    try:
        aprx = arcpy.mp.ArcGISProject("CURRENT")
    except Exception:
        return None

    candidate = None
    for m in aprx.listMaps():
        for lyr in m.listLayers():
            if lyr.name == layer_name:
                candidate = lyr
                break
            if getattr(lyr, "longName", None) == layer_name:
                candidate = lyr
                break
        if candidate:
            break

    if candidate is None:
        return None

    try:
        cp = candidate.connectionProperties
        if cp:
            ci = cp.get("connection_info", {}) or {}
            ds = cp.get("dataset") or ci.get("dataset")
            db = (cp.get("database") or ci.get("database")
                  or ci.get("server"))
            if ds:
                if os.path.isabs(str(ds)) and str(ds).lower().endswith(".nc"):
                    if os.path.isfile(ds):
                        return ds
                if str(ds).lower().endswith(".nc") and db and os.path.isdir(db):
                    full = os.path.join(db, ds)
                    if os.path.isfile(full):
                        return full
    except Exception:
        pass

    try:
        src = candidate.dataSource
        if src:
            if src.lower().endswith(".nc") and os.path.isfile(src):
                return src
            parts = re.split(r"[\\/]", src)
            for j in range(len(parts), 0, -1):
                cand = os.sep.join(parts[:j])
                if cand.lower().endswith(".nc") and os.path.isfile(cand):
                    return cand
    except Exception:
        pass

    m = re.match(r"^(.+?\.nc)[_:]", layer_name)
    if m:
        nc_name = m.group(1)
        search_dirs = []
        try:
            if aprx.homeFolder:
                search_dirs.append(aprx.homeFolder)
        except Exception:
            pass
        if arcpy.env.workspace:
            search_dirs.append(arcpy.env.workspace)
        for d in search_dirs:
            if d and os.path.isdir(d):
                full = os.path.join(d, nc_name)
                if os.path.isfile(full):
                    return full

    return None


def resolve_nc_inputs(file_param_raw, layer_param_raw, messages):
    """Combine file paths and layer references into a unified .nc path list."""
    nc_files = []

    if file_param_raw:
        for p in file_param_raw.split(";"):
            p = p.strip().strip("'\"")
            if p:
                nc_files.append(p)

    if layer_param_raw:
        for lyr in layer_param_raw.split(";"):
            lyr = lyr.strip().strip("'\"")
            if not lyr:
                continue
            path = get_layer_source_nc(lyr)
            if path:
                nc_files.append(path)
            else:
                messages.addWarningMessage(
                    "Could not resolve layer '{}' to a .nc file on disk. "
                    "Add the file via the file picker instead.".format(lyr)
                )

    seen = set()
    deduped = []
    for f in nc_files:
        key = os.path.normcase(os.path.abspath(f))
        if key not in seen:
            seen.add(key)
            deduped.append(f)
    return deduped


# ===========================================================================
# Fast conversion engine (netCDF4 + GDAL)
# ===========================================================================

def detect_variable_netcdf4(ds):
    """Find the HDD/CDD data variable in a netCDF4 Dataset."""
    candidates = []
    for vname, v in ds.variables.items():
        if vname in ds.dimensions:
            continue
        if vname in ("time_bnds", "lat_bnds", "lon_bnds", "crs"):
            continue
        if len(v.dimensions) >= 2:
            candidates.append(vname)

    if not candidates:
        raise ValueError("No data variable found in NetCDF.")

    priorities = [v for v in candidates if v.lower() in ("hdd", "cdd")]
    if not priorities:
        priorities = [v for v in candidates if "degree" in v.lower()]
    return priorities[0] if priorities else candidates[0]


def get_years_from_nc(ds, n_t, fallback_start_year):
    """Read time dimension and return a list of integer years.

    Uses netCDF4.num2date which handles CF calendars including 360-day.
    Falls back to start_year + i if num2date fails.
    """
    time_var = ds.variables["time"]
    try:
        units = time_var.units
        calendar = getattr(time_var, "calendar", "standard")
        times = netCDF4.num2date(
            time_var[:], units=units, calendar=calendar,
            only_use_cftime_datetimes=False,
            only_use_python_datetimes=False,
        )
        years = [int(t.year) for t in times]
        return years, "num2date"
    except Exception as e:
        if fallback_start_year is not None:
            return [fallback_start_year + i for i in range(n_t)], "filename"
        raise RuntimeError("Could not parse time axis: {}".format(e))


def compute_geotransform(lon, lat, lat_ascending):
    """GDAL geotransform tuple for north-up output: (origin_x, dx, 0, origin_y, 0, -dy)."""
    cs_x = float(abs(lon[1] - lon[0]))
    cs_y = float(abs(lat[1] - lat[0]))
    origin_x = float(lon[0]) - cs_x / 2.0
    if lat_ascending:
        origin_y = float(lat[-1]) + cs_y / 2.0
    else:
        origin_y = float(lat[0]) + cs_y / 2.0
    return (origin_x, cs_x, 0.0, origin_y, 0.0, -cs_y)


def write_geotiff_gdal(arr, out_path, geotransform, srs_wkt, nodata=NODATA_VAL):
    """Write a 2D NumPy array as a compressed GeoTIFF via GDAL."""
    driver = gdal.GetDriverByName("GTiff")
    nrows, ncols = arr.shape
    ds_out = driver.Create(
        out_path, ncols, nrows, 1, gdal.GDT_Float32, options=GTIFF_OPTIONS
    )
    ds_out.SetGeoTransform(geotransform)
    ds_out.SetProjection(srs_wkt)
    band = ds_out.GetRasterBand(1)
    band.WriteArray(arr.astype(np.float32))
    band.SetNoDataValue(float(nodata))
    band.FlushCache()
    ds_out = None  # close


def compute_bbox_slices(lat, lon, clip_extent, buffer_cells=1):
    """Compute (lat_slice, lon_slice) covering the clip extent bbox."""
    cs_x = abs(lon[1] - lon[0])
    cs_y = abs(lat[1] - lat[0])
    buf_x = buffer_cells * cs_x
    buf_y = buffer_cells * cs_y

    lon_mask = (lon >= clip_extent.XMin - buf_x) & (lon <= clip_extent.XMax + buf_x)
    lat_mask = (lat >= clip_extent.YMin - buf_y) & (lat <= clip_extent.YMax + buf_y)

    lon_idx = np.where(lon_mask)[0]
    lat_idx = np.where(lat_mask)[0]

    if len(lon_idx) == 0 or len(lat_idx) == 0:
        raise ValueError(
            "Clip extent does not intersect NetCDF grid. Clip extent: {} {} {} {}, "
            "NetCDF lon range: {}-{}, lat range: {}-{}".format(
                clip_extent.XMin, clip_extent.YMin,
                clip_extent.XMax, clip_extent.YMax,
                lon[0], lon[-1], lat[0], lat[-1]
            )
        )

    return (slice(lat_idx[0], lat_idx[-1] + 1),
            slice(lon_idx[0], lon_idx[-1] + 1))


# Cache: polygon mask per unique grid signature (avoids re-rasterising
# when many files share the same grid)
_MASK_CACHE = {}


def grid_signature(lat, lon):
    """Hashable key identifying a grid by extent + size."""
    return (len(lat), len(lon),
            round(float(lat[0]), 6), round(float(lat[-1]), 6),
            round(float(lon[0]), 6), round(float(lon[-1]), 6))


def rasterise_clip_mask(clip_fc, lat, lon, lat_ascending):
    """Rasterise the clip polygon to match the (sliced) lat/lon grid.

    Returns a boolean numpy array in the same orientation as the data
    arrays (lat-native: row 0 = lat[0]). True where the polygon covers.
    """
    cs_x = float(abs(lon[1] - lon[0]))
    cs_y = float(abs(lat[1] - lat[0]))

    left = float(lon[0]) - cs_x / 2.0
    right = float(lon[-1]) + cs_x / 2.0
    if lat_ascending:
        bottom = float(lat[0]) - cs_y / 2.0
        top = float(lat[-1]) + cs_y / 2.0
    else:
        bottom = float(lat[-1]) - cs_y / 2.0
        top = float(lat[0]) + cs_y / 2.0

    tmp = "in_memory/_ecde_mask_temp"
    if arcpy.Exists(tmp):
        arcpy.management.Delete(tmp)

    # PolygonToRaster requires a real field name, not the "OID@" token
    # (token resolution works only for SearchCursor / CalculateField).
    # Resolve dynamically: shapefile -> "FID", File GDB -> "OBJECTID".
    oid_field = arcpy.Describe(clip_fc).OIDFieldName

    with arcpy.EnvManager(
        extent="{} {} {} {}".format(left, bottom, right, top),
        cellSize=cs_x,
        snapRaster=None,
        outputCoordinateSystem=arcpy.SpatialReference(4326),
    ):
        arcpy.conversion.PolygonToRaster(
            in_features=clip_fc,
            value_field=oid_field,
            out_rasterdataset=tmp,
            cell_assignment="CELL_CENTER",
            cellsize=cs_x,
        )
        # arcpy returns array in north-up orientation (row 0 = highest lat)
        mask_arr = arcpy.RasterToNumPyArray(tmp, nodata_to_value=-1)
        arcpy.management.Delete(tmp)

    # Flip to lat-native orientation if needed
    if lat_ascending:
        mask_arr = mask_arr[::-1, :]

    # Trim or pad to exact grid size if there's a 1-row mismatch
    nrows_target, ncols_target = len(lat), len(lon)
    nrows_got, ncols_got = mask_arr.shape
    if (nrows_got, ncols_got) != (nrows_target, ncols_target):
        # Defensive: pad or crop to expected shape
        result = np.full((nrows_target, ncols_target), -1, dtype=mask_arr.dtype)
        r = min(nrows_got, nrows_target)
        c = min(ncols_got, ncols_target)
        result[:r, :c] = mask_arr[:r, :c]
        mask_arr = result

    return mask_arr != -1


def get_or_compute_mask(clip_fc, lat, lon, lat_ascending):
    if clip_fc is None:
        return None
    sig = grid_signature(lat, lon)
    if sig not in _MASK_CACHE:
        _MASK_CACHE[sig] = rasterise_clip_mask(clip_fc, lat, lon, lat_ascending)
    return _MASK_CACHE[sig]


def convert_nc_to_geotiffs(nc_path, out_folder, user_variable=None,
                           year_min=None, year_max=None,
                           clip_fc=None, messages=None,
                           verbose=False, advance_progressor=None):
    """Convert one NetCDF to per-year GeoTIFFs using netCDF4 + GDAL.

    With clip_fc:
      - Spatially slices the read to the bbox of the polygon (typically
        ~1% of the European grid for Portugal).
      - Applies a binary polygon mask (cached across files with the
        same grid) to preserve exact shape.
    Without clip_fc:
      - Reads the full extent, no mask.
    """
    def log(msg, level="msg"):
        if messages is None:
            arcpy.AddMessage(msg)
            return
        if level == "warn":
            messages.addWarningMessage(msg)
        elif level == "err":
            messages.addErrorMessage(msg)
        else:
            messages.addMessage(msg)

    t_start = _time.time()

    meta = parse_nc_filename(nc_path)
    log("  parsed: origin={}, agg={}, exp={}, rcm={}, gcm={}, ens={}".format(
        meta["origin"], meta["agg"], meta["experiment"] or "n/a",
        meta["rcm"] or "n/a", meta["gcm"] or "n/a", meta["ensemble"] or "n/a"
    ))

    ds = netCDF4.Dataset(nc_path)
    try:
        var_name = user_variable if user_variable else detect_variable_netcdf4(ds)
        var = ds.variables[var_name]
        lat = ds.variables["lat"][:]
        lon = ds.variables["lon"][:]
        n_t = var.shape[0]

        lat_ascending = bool(lat[1] > lat[0]) if len(lat) > 1 else True

        log("  variable='{}', shape={}, lat_ascending={}".format(
            var_name, var.shape, lat_ascending
        ))

        # Time parsing via netCDF4.num2date (handles all CF calendars)
        fallback_start = infer_start_year_from_filename(nc_path)
        years, year_source = get_years_from_nc(ds, n_t, fallback_start)
        log("  years: {} to {} ({} steps, source={})".format(
            years[0], years[-1], n_t, year_source
        ))

        # Spatial filtering
        if clip_fc:
            clip_extent = arcpy.Describe(clip_fc).extent
            lat_slice, lon_slice = compute_bbox_slices(lat, lon, clip_extent)
            sliced_lat = lat[lat_slice]
            sliced_lon = lon[lon_slice]
            log("  bbox slice: {} rows x {} cols (from {} x {} full grid)".format(
                len(sliced_lat), len(sliced_lon), len(lat), len(lon)
            ))
            mask = get_or_compute_mask(clip_fc, sliced_lat, sliced_lon, lat_ascending)
        else:
            lat_slice = slice(None)
            lon_slice = slice(None)
            sliced_lat = lat
            sliced_lon = lon
            mask = None

        # Get NoData / fill value from variable attributes
        fill_value = None
        for attr in ("_FillValue", "missing_value", "FillValue"):
            if hasattr(var, attr):
                fill_value = float(getattr(var, attr))
                break

        # Geotransform and CRS (computed once per file)
        gt = compute_geotransform(sliced_lon, sliced_lat, lat_ascending)
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)
        srs_wkt = srs.ExportToWkt()

        # Iterate time slices
        written = 0
        skipped = 0
        for i in range(n_t):
            year = years[i]
            if year_min is not None and year < year_min:
                skipped += 1
                if advance_progressor:
                    advance_progressor()
                continue
            if year_max is not None and year > year_max:
                skipped += 1
                if advance_progressor:
                    advance_progressor()
                continue

            try:
                arr = var[i, lat_slice, lon_slice]
                # Handle masked arrays from netCDF4
                if isinstance(arr, np.ma.MaskedArray):
                    arr = arr.filled(NODATA_VAL)
                arr = np.asarray(arr, dtype=np.float32)

                # Apply fill-value sentinel
                if fill_value is not None:
                    arr = np.where(arr == np.float32(fill_value), NODATA_VAL, arr)

                # Apply polygon mask
                if mask is not None:
                    arr = np.where(mask, arr, NODATA_VAL)

                # Flip to north-up for GeoTIFF output
                if lat_ascending:
                    arr_out = arr[::-1, :]
                else:
                    arr_out = arr

                out_name = build_output_name(meta, year)
                out_path = os.path.join(out_folder, out_name)
                write_geotiff_gdal(arr_out, out_path, gt, srs_wkt)

                written += 1
                if verbose:
                    log("    wrote {} ({}/{})".format(out_name, written, n_t - skipped))

            except Exception as e:
                log("  ERROR at index {} (year {}): {}".format(i, year, e), "warn")
            finally:
                if advance_progressor:
                    advance_progressor()

        elapsed = _time.time() - t_start
        per_slice = elapsed / written if written else 0
        log("  wrote {} GeoTIFFs in {:.1f}s ({:.2f}s/slice)".format(
            written, elapsed, per_slice
        ))

        # Build records for manifest
        records = []
        for i in range(n_t):
            year = years[i]
            if year_min is not None and year < year_min:
                continue
            if year_max is not None and year > year_max:
                continue
            out_name = build_output_name(meta, year)
            out_path = os.path.join(out_folder, out_name)
            if os.path.isfile(out_path):
                records.append({
                    "raster": out_path,
                    "variable": meta["variable"],
                    "origin": meta["origin"],
                    "agg": meta["agg"],
                    "experiment": meta["experiment"],
                    "gcm": meta["gcm"], "rcm": meta["rcm"],
                    "ensemble": meta["ensemble"],
                    "year": year,
                    "source_nc": os.path.basename(nc_path),
                })

        return records

    finally:
        ds.close()


def count_time_slices(nc_path):
    """Count time steps in a NetCDF for progressor sizing."""
    try:
        ds = netCDF4.Dataset(nc_path)
        try:
            return ds.dimensions["time"].size
        finally:
            ds.close()
    except Exception:
        return 0


# ===========================================================================
# Ensemble statistics helpers (Tool 3) -- unchanged from v0.3
# ===========================================================================

def find_tifs(folder):
    tifs = []
    for root, _dirs, files in os.walk(folder):
        for f in files:
            if f.lower().endswith(".tif"):
                tifs.append(os.path.join(root, f))
    return tifs


def categorise_tifs(tifs):
    projections, reanalysis, unparsable = [], [], []
    for r in tifs:
        meta = parse_output_filename(r)
        if not meta:
            unparsable.append(r)
        elif meta["origin"] == "projections":
            projections.append((r, meta))
        elif meta["origin"] == "reanalysis":
            reanalysis.append((r, meta))
        else:
            unparsable.append(r)
    return projections, reanalysis, unparsable


def group_by_year_scenario(projection_items):
    groups = {}
    for r, meta in projection_items:
        key = (meta["variable"], meta["experiment"], meta["year"])
        groups.setdefault(key, []).append(r)
    return groups


def raster_to_array_nan(raster_path):
    r = arcpy.Raster(raster_path)
    nd = r.noDataValue
    arr = arcpy.RasterToNumPyArray(
        r, nodata_to_value=NODATA_VAL
    ).astype(np.float32)
    arr = np.where(arr == NODATA_VAL, np.nan, arr)
    if nd is not None:
        try:
            arr = np.where(arr == float(nd), np.nan, arr)
        except (TypeError, ValueError):
            pass
    return arr


def array_to_raster(arr, ref_raster, out_path, sr):
    arr_clean = np.where(np.isnan(arr), NODATA_VAL, arr).astype(np.float32)
    out_r = arcpy.NumPyArrayToRaster(
        arr_clean,
        lower_left_corner=arcpy.Point(
            ref_raster.extent.XMin, ref_raster.extent.YMin
        ),
        x_cell_size=ref_raster.meanCellWidth,
        y_cell_size=ref_raster.meanCellHeight,
        value_to_nodata=NODATA_VAL,
    )
    out_r.save(out_path)
    arcpy.management.DefineProjection(out_path, sr)


STAT_FUNCS = {
    "mean":     lambda s: np.nanmean(s, axis=0),
    "median":   lambda s: np.nanmedian(s, axis=0),
    "p10":      lambda s: np.nanpercentile(s, 10, axis=0),
    "p90":      lambda s: np.nanpercentile(s, 90, axis=0),
    "std":      lambda s: np.nanstd(s, axis=0),
    "min":      lambda s: np.nanmin(s, axis=0),
    "max":      lambda s: np.nanmax(s, axis=0),
    "range":    lambda s: np.nanmax(s, axis=0) - np.nanmin(s, axis=0),
    "n_models": lambda s: np.sum(~np.isnan(s), axis=0).astype(np.float32),
}


def compute_year_ensemble(rasters, out_folder, var, exp, year, stats):
    stack = np.stack([raster_to_array_nan(r) for r in rasters], axis=0)
    ref = arcpy.Raster(rasters[0])
    sr = ref.spatialReference
    outputs = {}
    for stat in stats:
        if stat not in STAT_FUNCS:
            continue
        arr = STAT_FUNCS[stat](stack)
        fname = "{var}_{exp}_ensemble-{stat}_{year}.tif".format(
            var=var, exp=exp, stat=stat, year=year
        )
        out_path = os.path.join(out_folder, fname)
        array_to_raster(arr, ref, out_path, sr)
        outputs[stat] = out_path
    return outputs


def compute_period_mean(groups, var, exp, year_start, year_end,
                        out_folder, messages):
    yearly_means = []
    for (v, e, year), rasters in groups.items():
        if v != var or e != exp:
            continue
        if year < year_start or year > year_end:
            continue
        stack = np.stack([raster_to_array_nan(r) for r in rasters], axis=0)
        yearly_means.append(np.nanmean(stack, axis=0))

    if not yearly_means:
        messages.addWarningMessage(
            "  no data for {} {} period {}-{}".format(
                var, exp, year_start, year_end
            )
        )
        return None

    period_arr = np.nanmean(np.stack(yearly_means, axis=0), axis=0)
    any_key = next(iter(k for k in groups if k[0] == var and k[1] == exp))
    ref = arcpy.Raster(groups[any_key][0])
    sr = ref.spatialReference
    fname = "{var}_{exp}_periodmean_{ys}-{ye}.tif".format(
        var=var, exp=exp, ys=year_start, ye=year_end
    )
    out_path = os.path.join(out_folder, fname)
    array_to_raster(period_arr, ref, out_path, sr)
    return out_path


# ===========================================================================
# Toolbox
# ===========================================================================

class Toolbox(object):
    def __init__(self):
        self.label = "ECDE HDD/CDD Tools"
        self.alias = "ecde_hdd_cdd"
        self.tools = [
            ConvertECDEtoGeoTIFF,
            InspectECDENetCDF,
            ComputeEnsembleStatistics,
        ]


# ---------------------------------------------------------------------------
# Tool 1: Conversion (netCDF4 + GDAL fast path)
# ---------------------------------------------------------------------------

class ConvertECDEtoGeoTIFF(object):
    def __init__(self):
        self.label = "1. Convert ECDE NetCDF to GeoTIFF"
        self.description = (
            "Convert Copernicus ECDE HDD/CDD NetCDFs to yearly GeoTIFFs "
            "using netCDF4 + GDAL (fast path). When a clip polygon is "
            "provided, the engine reads only the bbox region and applies "
            "the exact polygon mask. Without clip, reads full extent."
        )
        self.canRunInBackground = False

    def getParameterInfo(self):
        p_files = arcpy.Parameter(
            displayName="Input NetCDF files (browse to disk)",
            name="in_nc_files", datatype="DEFile",
            parameterType="Optional", direction="Input", multiValue=True,
        )
        p_files.filter.list = ["nc"]

        p_layers = arcpy.Parameter(
            displayName="OR input multidimensional raster layers (from map)",
            name="in_nc_layers", datatype="GPRasterLayer",
            parameterType="Optional", direction="Input", multiValue=True,
        )

        p_out = arcpy.Parameter(
            displayName="Output folder",
            name="out_folder", datatype="DEFolder",
            parameterType="Required", direction="Input",
        )
        p_var = arcpy.Parameter(
            displayName="Variable name (blank = auto-detect: hdd/cdd/heating_degree_days)",
            name="variable", datatype="GPString",
            parameterType="Optional", direction="Input",
        )
        p_ymin = arcpy.Parameter(
            displayName="Year minimum (inclusive, optional)",
            name="year_min", datatype="GPLong",
            parameterType="Optional", direction="Input",
        )
        p_ymax = arcpy.Parameter(
            displayName="Year maximum (inclusive, optional)",
            name="year_max", datatype="GPLong",
            parameterType="Optional", direction="Input",
        )
        p_clip = arcpy.Parameter(
            displayName="Clip polygon feature class (optional - enables spatial slicing + masking)",
            name="clip_fc", datatype="GPFeatureLayer",
            parameterType="Optional", direction="Input",
        )
        p_clip.filter.list = ["Polygon"]

        p_manifest = arcpy.Parameter(
            displayName="Write manifest CSV",
            name="write_manifest", datatype="GPBoolean",
            parameterType="Optional", direction="Input",
        )
        p_manifest.value = True

        p_verbose = arcpy.Parameter(
            displayName="Verbose logging (one message per TIFF written)",
            name="verbose", datatype="GPBoolean",
            parameterType="Optional", direction="Input",
        )
        p_verbose.value = False

        return [p_files, p_layers, p_out, p_var, p_ymin, p_ymax,
                p_clip, p_manifest, p_verbose]

    def isLicensed(self):
        return True

    def updateMessages(self, parameters):
        p_files = parameters[0]
        p_layers = parameters[1]
        if not p_files.value and not p_layers.value:
            p_files.setErrorMessage(
                "Provide at least one input source: files via the picker "
                "or layers from the map."
            )

    def execute(self, parameters, messages):
        file_raw = parameters[0].valueAsText
        layer_raw = parameters[1].valueAsText
        out_folder = parameters[2].valueAsText
        user_variable = parameters[3].valueAsText or None
        year_min = parameters[4].value
        year_max = parameters[5].value
        clip_fc = parameters[6].valueAsText if parameters[6].value else None
        write_manifest = parameters[7].value
        verbose = bool(parameters[8].value)

        nc_files = resolve_nc_inputs(file_raw, layer_raw, messages)

        if not nc_files:
            messages.addErrorMessage("No valid NetCDF inputs resolved. Aborting.")
            return

        if not os.path.isdir(out_folder):
            os.makedirs(out_folder)

        # Pre-scan to size progressor
        slice_counts = [count_time_slices(nc) for nc in nc_files]
        total_slices = sum(slice_counts)

        messages.addMessage("Processing {} NetCDF file(s), {} total time slices".format(
            len(nc_files), total_slices
        ))
        messages.addMessage("Engine: netCDF4 + GDAL (fast path, v0.4)")
        messages.addMessage("Output: {}".format(out_folder))
        if year_min is not None or year_max is not None:
            messages.addMessage("Year filter: {} to {}".format(
                year_min if year_min is not None else "min",
                year_max if year_max is not None else "max",
            ))
        if clip_fc:
            messages.addMessage(
                "Clip: {} (will use bbox slicing + polygon mask)".format(clip_fc)
            )
        if verbose:
            messages.addMessage("Verbose logging enabled")

        # Reset mask cache between runs
        _MASK_CACHE.clear()

        progressor_state = {"pos": 0}

        def advance():
            progressor_state["pos"] += 1
            arcpy.SetProgressorPosition(progressor_state["pos"])

        arcpy.SetProgressor(
            type="step",
            message="Processing NetCDF time slices...",
            min_range=0,
            max_range=max(total_slices, 1),
            step_value=1,
        )

        all_records = []
        t_total_start = _time.time()
        try:
            for idx, nc in enumerate(nc_files):
                fname = os.path.basename(nc)
                arcpy.SetProgressorLabel(
                    "File {}/{}: {}".format(idx + 1, len(nc_files), fname)
                )
                messages.addMessage("--- [{}/{}] {}".format(
                    idx + 1, len(nc_files), fname
                ))
                try:
                    recs = convert_nc_to_geotiffs(
                        nc, out_folder,
                        user_variable=user_variable,
                        year_min=year_min, year_max=year_max,
                        clip_fc=clip_fc, messages=messages,
                        verbose=verbose,
                        advance_progressor=advance,
                    )
                    all_records.extend(recs)
                except Exception as e:
                    messages.addErrorMessage("FAILED on {}: {}".format(fname, e))
        finally:
            arcpy.ResetProgressor()

        if write_manifest and all_records:
            manifest_path = os.path.join(out_folder, "manifest.csv")
            fields = ["raster", "variable", "origin", "agg", "experiment",
                      "gcm", "rcm", "ensemble", "year", "source_nc"]
            with open(manifest_path, "w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(f, fieldnames=fields)
                writer.writeheader()
                writer.writerows(all_records)
            messages.addMessage("Manifest: {}".format(manifest_path))

        total_elapsed = _time.time() - t_total_start
        avg_per_slice = total_elapsed / len(all_records) if all_records else 0
        messages.addMessage(
            "Done. {} GeoTIFFs from {} NetCDF(s) in {:.1f}s ({:.2f}s/slice avg).".format(
                len(all_records), len(nc_files), total_elapsed, avg_per_slice
            )
        )


# ---------------------------------------------------------------------------
# Tool 2: Inspect (unchanged from v0.3)
# ---------------------------------------------------------------------------

class InspectECDENetCDF(object):
    def __init__(self):
        self.label = "2. Inspect ECDE NetCDF (debug)"
        self.description = (
            "Print dimensions, variables, parsed metadata and raw time "
            "value types for ECDE NetCDFs. Use before Tool 1 to verify "
            "structure. Resolves loaded multidim layers via arcpy.mp."
        )
        self.canRunInBackground = False

    def getParameterInfo(self):
        p_files = arcpy.Parameter(
            displayName="Input NetCDF files (browse to disk)",
            name="in_nc_files", datatype="DEFile",
            parameterType="Optional", direction="Input", multiValue=True,
        )
        p_files.filter.list = ["nc"]
        p_layers = arcpy.Parameter(
            displayName="OR input multidimensional raster layers (from map)",
            name="in_nc_layers", datatype="GPRasterLayer",
            parameterType="Optional", direction="Input", multiValue=True,
        )
        return [p_files, p_layers]

    def isLicensed(self):
        return True

    def updateMessages(self, parameters):
        p_files = parameters[0]
        p_layers = parameters[1]
        if not p_files.value and not p_layers.value:
            p_files.setErrorMessage(
                "Provide at least one input source: files via the picker "
                "or layers from the map."
            )

    def execute(self, parameters, messages):
        file_raw = parameters[0].valueAsText
        layer_raw = parameters[1].valueAsText
        nc_files = resolve_nc_inputs(file_raw, layer_raw, messages)

        for nc in nc_files:
            messages.addMessage("=" * 70)
            messages.addMessage(os.path.basename(nc))
            messages.addMessage("=" * 70)

            meta = parse_nc_filename(nc)
            for k, v in meta.items():
                messages.addMessage("  meta.{} = {}".format(k, v))

            inferred = infer_start_year_from_filename(nc)
            messages.addMessage(
                "  inferred start year (from filename): {}".format(inferred)
            )

            try:
                props = arcpy.NetCDFFileProperties(nc)
                dims = list(props.getDimensions())
                vars_ = list(props.getVariables())
                messages.addMessage("  dimensions: {}".format(dims))
                messages.addMessage("  variables : {}".format(vars_))
                for d in dims:
                    size = props.getDimensionSize(d)
                    messages.addMessage("    {}: size={}".format(d, size))
                    for k in range(min(3, size)):
                        try:
                            val = props.getDimensionValue(d, k)
                            messages.addMessage(
                                "      [{}]: type={}, repr={}".format(
                                    k, type(val).__name__, repr(val)[:80]
                                )
                            )
                            if d.lower() in ("time", "t", "stdtime"):
                                parsed = extract_year_from_value(val)
                                messages.addMessage(
                                    "          extract_year_from_value -> {}".format(
                                        parsed
                                    )
                                )
                        except Exception as e:
                            messages.addMessage(
                                "      [{}]: read failed: {}".format(k, e)
                            )
                try:
                    v, xd, yd, td = discover_structure(props)
                    messages.addMessage(
                        "  detected: variable='{}', x='{}', y='{}', t='{}'".format(
                            v, xd, yd, td
                        )
                    )
                except Exception as e:
                    messages.addWarningMessage(
                        "  auto-detect failed: {}".format(e)
                    )
            except Exception as e:
                messages.addErrorMessage("Inspection failed: {}".format(e))


# ---------------------------------------------------------------------------
# Tool 3: Ensemble Statistics (projections only) -- unchanged from v0.3
# ---------------------------------------------------------------------------

class ComputeEnsembleStatistics(object):
    def __init__(self):
        self.label = "3. Compute Ensemble Statistics (projections only)"
        self.description = (
            "Compute pixel-wise multi-model ensemble statistics from "
            "PROJECTION GeoTIFFs produced by Tool 1. Reanalysis TIFFs "
            "are excluded because ensemble stats require multiple models. "
            "Groups rasters by (variable, scenario, year). Optionally "
            "computes 30-year climatological period means of the ensemble mean."
        )
        self.canRunInBackground = False

    def getParameterInfo(self):
        p_in = arcpy.Parameter(
            displayName="Input folder (output of Tool 1, with projection TIFFs)",
            name="in_folder", datatype="DEFolder",
            parameterType="Required", direction="Input",
        )
        p_out = arcpy.Parameter(
            displayName="Output folder for ensemble statistics",
            name="out_folder", datatype="DEFolder",
            parameterType="Required", direction="Input",
        )
        p_stats = arcpy.Parameter(
            displayName="Statistics to compute",
            name="stats", datatype="GPString",
            parameterType="Required", direction="Input", multiValue=True,
        )
        p_stats.filter.type = "ValueList"
        p_stats.filter.list = list(STAT_FUNCS.keys())
        p_stats.values = ["mean", "median", "p10", "p90", "std"]

        p_ymin = arcpy.Parameter(
            displayName="Year minimum (inclusive, optional)",
            name="year_min", datatype="GPLong",
            parameterType="Optional", direction="Input",
        )
        p_ymax = arcpy.Parameter(
            displayName="Year maximum (inclusive, optional)",
            name="year_max", datatype="GPLong",
            parameterType="Optional", direction="Input",
        )
        p_periods = arcpy.Parameter(
            displayName="30-year period means (start-end pairs, e.g. '2041-2070;2071-2100')",
            name="periods", datatype="GPString",
            parameterType="Optional", direction="Input",
        )
        p_periods.value = "2011-2040;2041-2070;2071-2100"

        return [p_in, p_out, p_stats, p_ymin, p_ymax, p_periods]

    def isLicensed(self):
        return True

    def execute(self, parameters, messages):
        in_folder = parameters[0].valueAsText
        out_folder = parameters[1].valueAsText
        stats_raw = parameters[2].valueAsText
        year_min = parameters[3].value
        year_max = parameters[4].value
        periods_raw = parameters[5].valueAsText or ""

        stats = [s.strip().strip("'\"") for s in stats_raw.split(";")
                 if s.strip()]

        if not os.path.isdir(out_folder):
            os.makedirs(out_folder)

        tifs = find_tifs(in_folder)
        messages.addMessage("Found {} GeoTIFFs in input folder".format(len(tifs)))

        projection_items, reanalysis_items, unparsable = categorise_tifs(tifs)
        messages.addMessage("  - {} projection TIFFs".format(len(projection_items)))
        messages.addMessage(
            "  - {} reanalysis TIFFs (excluded - ensemble stats need multi-model)".format(
                len(reanalysis_items)
            )
        )
        if unparsable:
            messages.addWarningMessage(
                "  - {} TIFFs with unrecognised filename pattern (skipped)".format(
                    len(unparsable)
                )
            )

        if not projection_items:
            messages.addErrorMessage(
                "No projection TIFFs found. This tool requires output from "
                "Tool 1 run on the 32 *-projections-*.nc files (HDD/CDD x "
                "8 GCM-RCM x 2 RCPs). Point the Input folder at that output "
                "and re-run."
            )
            return

        groups = group_by_year_scenario(projection_items)
        if year_min is not None or year_max is not None:
            ymin = year_min if year_min is not None else -10000
            ymax = year_max if year_max is not None else 10000
            groups = {k: v for k, v in groups.items() if ymin <= k[2] <= ymax}

        messages.addMessage("Groups (var, scenario, year): {}".format(len(groups)))
        messages.addMessage("Statistics: {}".format(", ".join(stats)))

        total = 0
        decades_logged = set()
        for (var, exp, year), rasters in sorted(groups.items()):
            if len(rasters) < 2:
                messages.addWarningMessage(
                    "Skipping {} {} {} (only {} model in group)".format(
                        var, exp, year, len(rasters)
                    )
                )
                continue
            try:
                outs = compute_year_ensemble(
                    rasters, out_folder, var, exp, year, stats
                )
                total += len(outs)
                decade = (year // 10) * 10
                if (var, exp, decade) not in decades_logged:
                    messages.addMessage(
                        "  {} {} {}s: {} models, {} stats per year".format(
                            var, exp, decade, len(rasters), len(outs)
                        )
                    )
                    decades_logged.add((var, exp, decade))
            except Exception as e:
                messages.addWarningMessage(
                    "Failed {} {} {}: {}".format(var, exp, year, e)
                )

        messages.addMessage("Per-year stats: {} rasters written".format(total))

        if periods_raw:
            periods = []
            for p in periods_raw.split(";"):
                p = p.strip().strip("'\"")
                m = re.match(r"^(\d{4})-(\d{4})$", p)
                if not m:
                    messages.addWarningMessage("Bad period spec: '{}'".format(p))
                    continue
                ys, ye = int(m.group(1)), int(m.group(2))
                if ys >= ye:
                    messages.addWarningMessage(
                        "Invalid period {}-{}: start year must be before end. Skipping.".format(
                            ys, ye
                        )
                    )
                    continue
                length = ye - ys + 1
                if length != 30:
                    messages.addMessage(
                        "Note: period {}-{} spans {} years (WMO climatological "
                        "normal is 30 years)".format(ys, ye, length)
                    )
                periods.append((ys, ye))

            if periods:
                messages.addMessage("Period means: {}".format(periods))
                variables = sorted({k[0] for k in groups})
                experiments = sorted({k[1] for k in groups})
                period_total = 0
                for var in variables:
                    for exp in experiments:
                        for ys, ye in periods:
                            out = compute_period_mean(
                                groups, var, exp, ys, ye,
                                out_folder, messages
                            )
                            if out:
                                period_total += 1
                                messages.addMessage(
                                    "  {} {} {}-{}: {}".format(
                                        var, exp, ys, ye,
                                        os.path.basename(out)
                                    )
                                )
                messages.addMessage(
                    "Period means: {} rasters written".format(period_total)
                )

        messages.addMessage("Done.")
