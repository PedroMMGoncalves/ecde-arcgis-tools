# ecde-arcgis-tools

ArcGIS Pro Python toolbox (`.pyt`) for working with Copernicus
**European Climate Data Explorer (ECDE)** Heating Degree Days (HDD) and
Cooling Degree Days (CDD) NetCDF data. Converts the published NetCDFs into
GIS-ready GeoTIFFs and computes multi-model ensemble statistics.

Data source: [Copernicus Climate Data Store, dataset `sis-ecde-climate-indicators`](https://cds.climate.copernicus.eu/datasets/sis-ecde-climate-indicators).

---

## Contents

- [What it does](#what-it-does)
- [Requirements](#requirements)
- [Installation](#installation)
- [Tool 1 — Convert ECDE NetCDF to GeoTIFF](#tool-1--convert-ecde-netcdf-to-geotiff)
- [Tool 2 — Inspect ECDE NetCDF (debug)](#tool-2--inspect-ecde-netcdf-debug)
- [Tool 3 — Compute Ensemble Statistics](#tool-3--compute-ensemble-statistics-projections-only)
- [Output naming](#output-naming)
- [Manifest CSV](#manifest-csv)
- [Ensemble coverage caveat](#ensemble-coverage-caveat)
- [Method / citation](#method--citation)
- [Citing this toolbox](#citing-this-toolbox)
- [Limitations](#limitations)
- [License](#license)

---

## What it does

The toolbox `ecde_hddcdd_tools_v1.pyt` exposes three tools:

1. **Convert ECDE NetCDF to GeoTIFF** — one GeoTIFF per year, with optional
   polygon clip and a manifest CSV for downstream ingestion.
2. **Inspect ECDE NetCDF (debug)** — prints dimensions, variables, parsed
   filename metadata, and raw time-value diagnostics.
3. **Compute Ensemble Statistics (projections only)** — pixel-wise stats
   (mean, median, p10, p90, std, min, max, range, n_models) across the
   8 GCM-RCM chains, grouped by (variable, scenario, year). Optionally
   computes 30-year climatological period means of the ensemble mean.

### Supported data

- **Reanalysis** (ERA5-based), yearly grid, 1940-2025, single deterministic
  dataset.
- **Projections** (EURO-CORDEX bias-adjusted), yearly grid, 8 GCM-RCM chains,
  RCP4.5 and RCP8.5. Each file contains historical + scenario merged
  (1950/1951/1970 to 2100 depending on chain).
- **Grid:** regular 0.25 degree WGS84 lat/lon (EPSG:4326), 271 x 185 over
  Europe. Not rotated pole.

---

## Requirements

- **ArcGIS Pro 3.x** (Windows). The toolbox runs inside the default
  `arcgispro-py3` Python environment.
- No extra installs. Only libraries that ship with ArcGIS Pro by default
  are used:
  - `arcpy`
  - `numpy`
  - `netCDF4`
  - `osgeo.gdal`

Raw EURO-CORDEX at 0.11 degrees (rotated pole) is **not** supported. Use the
ECDE 0.25 degree bias-adjusted product only.

---

## Installation

1. Clone or download this repository:

   ```bash
   git clone https://github.com/PedroMMGoncalves/ecde-arcgis-tools.git
   ```

2. In ArcGIS Pro, open the **Catalog** pane.
3. Right-click **Toolboxes** -> **Add Toolbox**.
4. Browse to `ArcGisPro_Toolbox/ecde_hddcdd_tools_v1.pyt` and add it.
5. The three tools appear under the toolbox.

To download the source NetCDF data:

1. Register at [https://cds.climate.copernicus.eu/](https://cds.climate.copernicus.eu/).
2. Open the [`sis-ecde-climate-indicators`](https://cds.climate.copernicus.eu/datasets/sis-ecde-climate-indicators)
   dataset page.
3. Select:
   - **Variable:** Heating degree days and/or Cooling degree days
   - **Origin:** Reanalysis and/or Projections
   - **Temporal aggregation:** Yearly
   - For projections: all 8 GCM-RCM chains, RCP4.5 and RCP8.5
4. Submit the request and download. Save the `.nc` files to a local folder.

---

## Tool 1 — Convert ECDE NetCDF to GeoTIFF

Converts ECDE NetCDFs to one compressed GeoTIFF per time slice (year), using
`netCDF4` + GDAL directly. With a clip polygon, only the polygon's bounding
box is read from each NetCDF and a binary polygon mask is applied
(rasterised once per unique grid signature, cached across files).

### Parameters

| Name | Type | Required | Description |
|---|---|---|---|
| Input NetCDF files | Multi-file | One of the two | ECDE `.nc` files on disk. |
| Input multidimensional raster layers | Multi-layer | One of the two | Same data, but already added to the map as multidim layers. |
| Output folder | Folder | Yes | GeoTIFFs (and `manifest.csv`) are written here. |
| Variable name | String | No | Leave blank for auto-detect (`hdd`, `cdd`, `heating_degree_days`, `cooling_degree_days`). |
| Year minimum | Long | No | Inclusive lower bound. |
| Year maximum | Long | No | Inclusive upper bound. |
| Clip polygon feature class | Polygon FC | No | Enables bbox slicing + exact polygon mask. |
| Write manifest CSV | Boolean | No | Default on. Writes `manifest.csv` summarising every TIFF. |
| Verbose logging | Boolean | No | One message per TIFF written. Off by default. |

### Output

- One GeoTIFF per year per file (DEFLATE-compressed, tiled, EPSG:4326).
- `manifest.csv` (optional) with one row per TIFF: variable, origin,
  experiment, rcm, gcm, ensemble, year, source NetCDF, output path.

### Typical run

For Portugal-clipped output of the full ECDE catalogue (1 reanalysis file +
32 projection files, 4732 GeoTIFFs):

```text
Engine: netCDF4 + GDAL (fast path, v1.0.0)
... wrote 4732 GeoTIFFs in ~9 seconds (~0.002 s/slice)
```

---

## Tool 2 — Inspect ECDE NetCDF (debug)

Prints structural diagnostics for one or more ECDE NetCDFs:

- Dimensions and their sizes
- Variables (including coordinate variables)
- Parsed filename metadata (indicator, origin, experiment, RCM, GCM, ensemble)
- Raw time-axis types, units, calendar, and the first/last decoded values

Use this when Tool 1 misbehaves, to see exactly what arcpy and `netCDF4`
read from the file.

---

## Tool 3 — Compute Ensemble Statistics (projections only)

Reads the projection GeoTIFFs produced by Tool 1, groups them by
`(variable, scenario, year)`, and computes pixel-wise statistics across the
8 GCM-RCM model dimension.

### Parameters

| Name | Type | Required | Description |
|---|---|---|---|
| Input folder | Folder | Yes | Output of Tool 1, containing projection TIFFs. |
| Output folder | Folder | Yes | Ensemble statistics are written here. |
| Statistics to compute | Multi-value | Yes | Any of `mean`, `median`, `p10`, `p90`, `std`, `min`, `max`, `range`, `n_models`. Default: `mean, median, p10, p90, std`. |
| Year minimum / maximum | Long | No | Restrict to a year window. |
| 30-year period means | String | No | Semicolon-separated `start-end` pairs. Default `2011-2040;2041-2070;2071-2100` (IPCC AR6 near/mid/long-term). |

Reanalysis TIFFs in the input folder are silently filtered out (the count
is reported in the log), since ensemble statistics across one model are
undefined.

### Output

- Per-year stats: `HDD_rcp45_2050_mean.tif`, `HDD_rcp45_2050_p90.tif`, ...
- Per-period means of the ensemble mean (when requested):
  `HDD_rcp45_2041-2070_mean.tif`.
- An `n_models` raster pixel-wise documents how many models contributed
  (see caveat below).

---

## Output naming

```text
Reanalysis:  HDD_reanalysis_yearly_2020.tif
Projections: HDD_rcp45_RACMO22E_EC-EARTH_r1i1p1_2050.tif
Ensemble:    HDD_rcp45_2050_mean.tif
Periods:     HDD_rcp45_2041-2070_mean.tif
```

The same applies with `CDD` in place of `HDD`. Tool 3 parses these
filenames to group rasters by `(variable, scenario, year)`.

---

## Manifest CSV

Tool 1 writes `manifest.csv` in the output folder. Each row corresponds to
one GeoTIFF and includes:

- `variable` (HDD or CDD)
- `origin` (reanalysis or projections)
- `experiment` (rcp45, rcp85, or reanalysis)
- `rcm`, `gcm`, `ensemble`
- `year`
- `source_nc` (input file path)
- `tif_path` (output GeoTIFF path)

Use it for:

- Loading into a Mosaic Dataset.
- Driving Zonal Statistics over administrative boundaries in ArcGIS / Pandas.
- Filtering by experiment / chain / year in downstream analysis.

---

## Ensemble coverage caveat

The 8 GCM-RCM chains in the ECDE projections have different temporal
coverage:

| Period | Models available |
|---|---|
| 1950 | 1 (CCLM4-8-17 / MPI-ESM-LR only) |
| 1951-1969 | 5 |
| **1970-2097** | **8 (full coverage)** |
| 2098-2100 | 6 (HadGEM2-ES chains end at 2098) |

For a constant 8-model ensemble, restrict analysis to **1970-2097**.
Outside that window, `n_models` varies, and percentile / std estimates
degrade. The `n_models` statistic raster documents this pixel-by-pixel.

---

## Method / citation

HDD and CDD are computed in the ECDE source product using the **Spinoni et
al. (2018)** method, with a base temperature of **15.5 deg C**.

> Spinoni J, Vogt JV, Barbosa P, Dosio A, McCormick N, Bigano A, Fuessel HM
> (2018). Changes of heating and cooling degree-days in Europe from 1981 to
> 2100. *International Journal of Climatology* 38: e191-e208.
> DOI: [10.1002/joc.5362](https://doi.org/10.1002/joc.5362)

For a representative literature range to sanity-check outputs over the
Iberian Peninsula:

> Andrade et al. (2021). *Atmosphere* 12:715. Portugal: ~800-2500 HDD/year
> depending on region; CDD considerably lower.

If you use this toolbox in published work, please cite both the dataset
(`sis-ecde-climate-indicators` on Copernicus CDS) and Spinoni et al. (2018),
and link to this repository (see [Citing this toolbox](#citing-this-toolbox)).

---

## Citing this toolbox

If this toolbox supports your published work, please cite it. A
[`CITATION.cff`](CITATION.cff) file is included so GitHub renders a
"Cite this repository" widget and reference managers (Zotero, Mendeley)
import the metadata automatically.

**Software citation:**

> Goncalves, P. (2026). *ecde-arcgis-tools: ArcGIS Pro toolbox for
> Copernicus ECDE HDD/CDD data* (Version 1.0.0) [Software]. Zenodo.
> DOI: *to be assigned on first Zenodo release*

The Zenodo DOI is a permanent identifier that resolves to a snapshot of
the source code archived independently of GitHub. Once minted it will be
inserted here.

**Please also cite the source dataset and method:**

- Copernicus Climate Change Service (2024). *Climate indicators for
  Europe from 1940 to 2100 derived from reanalysis and climate
  projections.* Climate Data Store.
  <https://cds.climate.copernicus.eu/datasets/sis-ecde-climate-indicators>
- Spinoni et al. (2018). *International Journal of Climatology* 38:
  e191-e208. DOI: [10.1002/joc.5362](https://doi.org/10.1002/joc.5362)

---

## Limitations

These are intentional non-features. Open an issue if you have a strong
use case before implementing.

- No support for raw EURO-CORDEX at 0.11 degrees (rotated pole grid
  requires reprojection logic that is out of scope).
- No QGIS or CLI front-end. Toolbox runs inside ArcGIS Pro only.
- No multiprocessing. Sequential processing is already fast enough
  (~9 seconds for 4732 GeoTIFFs at 0.25 deg, clipped to Portugal).
- No delta / anomaly calculation (future minus reference period).
- No automatic Zonal Statistics over administrative boundaries. Do this
  downstream with ArcGIS or Pandas using the manifest.
- No Mosaic Dataset auto-creation. Create the mosaic and ingest GeoTIFFs
  manually if desired.

---

## License

Apache License 2.0. See [LICENSE](LICENSE).

Copyright 2026 Pedro Goncalves, LNEG (Laboratorio Nacional de Energia e
Geologia).
