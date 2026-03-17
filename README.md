# SpotiPy: Solar Active Region Analysis Framework

[![JOSS Status](https://joss.theoj.org/papers/a6818174542289139f72787834526017/status.svg)](https://joss.theoj.org/papers/a6818174542289139f72787834526017)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

SpotiPy is a modular Python framework for the analysis of sunspot center-to-limb variations (CLV). It provides an interactive pipeline for downloading, aligning, segmenting, and extracting physical quantities from SDO/HMI and AIA observations of solar active regions across five observables:

1. **Intensity ($I_c$)** — SDO/HMI Continuum
2. **Magnetogram ($M$)** — SDO/HMI LOS Magnetic Field
3. **Dopplergram ($V$)** — SDO/HMI Velocity
4. **Line Depth ($L_d$)** — SDO/HMI
5. **Line Width ($L_w$)** — SDO/HMI

---

## Key Features

SpotiPy is designed to be fully modular. You can run the complete pipeline end-to-end or import individual components into your own analysis scripts.

* **Multi-Observable Pipeline**
  Interactive processing of co-aligned intensity, magnetic field, Doppler velocity, line depth, and line width data.

* **CLV Scatter Plots with Polynomial Fit**
  Visualizes center-to-limb variation trends across heliocentric angles.
  
* **Time-Summed Strip Visualization**
  Generates a co-aligned, time-summed strip image of the solar disk showing the sunspot trajectory over the full observation window.

* **Consolidated Data Output**
  All cropped frames, segmentation masks, timestamps, and heliographic coordinates are saved into a single compressed `.npz` archive for easy use.

* **Modular Architecture**

  * `spotipy.downloading` — Batch data retrieval with automatic JSOC authentication
  * `spotipy.aligning` — Sub-pixel co-alignment using WCS reprojection
  * `spotipy.limbdarkening_removal` — Radial profile-based limb darkening correction for AIA data
  * `spotipy.segmentation` — Intensity-based umbra, penumbra, plage, and network masking
  * `spotipy.tracking` — Differential rotation modeling, region tracking, and strip visualization
  * `spotipy.gui_tools` — Interactive GUI for manual selection of the starting coordinates of a solar feature

---

## Dependencies

SpotiPy requires Python 3.9 or later. The following packages are installed automatically:

| Package | Purpose |
|---|---|
| `numpy` | Array operations and data handling |
| `matplotlib` | Plotting and visualization |
| `astropy` | FITS I/O, WCS, and time handling |
| `sunpy` | Solar coordinate transforms and differential rotation |
| `scipy` | Image processing and polynomial fitting |
| `reproject` | WCS-based image reprojection for AIA-to-HMI alignment |

---

## Installation

Clone the repository and install it in **editable mode**. This is recommended for research use, as changes to the code are applied immediately.

```bash
git clone https://github.com/Emily-Joe/SpotiPy.git
cd SpotiPy
pip install -e .
```

---

## Quick Start (Library Mode)

You can import SpotiPy tools directly into your own Python scripts:

```python
from spotipy.segmentation import get_masks
import numpy as np

# Generate segmentation masks from image data
# (Assuming `image_data` is a NumPy array you have already loaded)
masks = get_masks(
    image_data,
    umbra_range=(10, 55),
    penumbra_range=(75, 120)
)

print(f"Detected {np.sum(masks['umbra'])} pixels of Umbra.")
```

---

## Usage (Pipeline Mode)

For full end-to-end analysis (**Download → Align → Segment → Analyze → Plot**), place `params.txt` in your working directory and run:

```bash
python3 run_analysis.py
```

The pipeline will interactively prompt you for each step (download, coordinate selection, mask generation, etc.).

> **Note:** Downloading data requires a JSOC-registered email address. Visit [JSOC Email Registration](http://jsoc.stanford.edu/ajax/register_email.html) to register before running the pipeline. An example `params.txt` configuration file is included in the repository.

---

## Configuration File (`params.txt`)

The pipeline is controlled through a simple text-based configuration file. All parameters are listed below:

 Parameter | Description | Example |
|---|---|---|
| `NOAA_NUMBER` | Active Region ID | `12673` |
| `START_DATE` | ISO format start time | `2017-09-06T12:00:00` |
| `EMAIL` | JSOC-registered email for data access | `you@example.com` |
| `DAYS` | Duration of observation window in days | `11` |
| `CADENCE` | Time step between frames in hours | `6` |
| `X_START_ARCSEC` | Initial x-position of the region in arcsec | `-921` |
| `Y_START_ARCSEC` | Initial y-position of the region in arcsec | `130` |
| `FRAME_SIZE` | Crop size in pixels (square) | `400` |
| `UMBRA_MIN` / `UMBRA_MAX` | Intensity threshold range for umbra | `10` / `55` |
| `PENUMBRA_MIN` / `PENUMBRA_MAX` | Intensity threshold range for penumbra | `75` / `120` |
| `PLAGE_EXCESS_PCT` | AIA plage detection threshold (% above quiet Sun) | `20.0` |
| `QUIET_SUN_TOL_PCT` | AIA network detection threshold (% above quiet Sun) | `15.0` |
| `MIN_PLAGE_AREA` | Minimum connected area in pixels for plage detection | `450` |
| `HMI_POLY_ORDER` | Polynomial order for HMI CLV fit | `2` |
| `AIA_POLY_ORDER` | Polynomial order for AIA CLV fit | `5` |
| `MIN_MU_FOR_FIT` | Minimum mu value included in CLV fit | `0.15` |

---

## Available Data for Download

The downloading module enables the search and download of data from public archives of Solar Dynamics Observatory's HMI and AIA instruments. It provides interface built on `sunpy`'s Fido client. Download of HMI and AIA data requires prior registration via email with the JSOC data center. Visit [JSOC Email Registration](http://jsoc.stanford.edu/ajax/register_email.html) and register your email address. This email address must also be specified in `params.txt`.

The following datasets are available for download:

### 1. HMI Observables

HMI generally provides data series with cadences of 45 or 720 seconds. It offers continuum intensity measurements in the region of the Fe I absorption line at 6173 Å, as well as a version with the effects of limb darkening removed. The limb-darkening-corrected dataset is only available at a cadence of 12 minutes. Additional observables include the line depth and line width of a Gaussian absorption line profile for the 6173.3 Å line, as well as line-of-sight (LOS) magnetograms and Dopplergrams.

| Continuum Intensity (Ic)   | Linewidth (Lw) | Linedepth (Ld) | LOS Magnetograms (M) | LOS Dopplergrams (V) |
|----------------------------|----------------|----------------|----------------------|----------------------|
| hmi.Ic_45s                 | hmi.Lw_45s     | hmi.Ld_45s     | hmi.M_45s            | hmi.V_45s            |
| hmi.Ic_720s                | hmi.Lw_720s    | hmi.Ld_720s    | hmi.M_720s           | hmi.V_720s           |
| hmi.Ic_nolimbDark_720s     |                |                |                      |                      |
| hmi.Ic_nolimbDark_720s_nrt |                |                |                      |                      |

### 2. AIA Observables

At Level 1, AIA observables cover three wavelength ranges (soft X-ray/extreme UV, UV, and visible), with ten total wavelength channels. The wavelengths are given in angstroms at the end of each series name. Each wavelength range is available at different cadences: every 12 seconds for EUV, every 24 seconds for UV, and once per hour for the visible range.

| Soft X-ray – Extreme UV (euv), 12s | UV (uv), 24s          | Visible (vis), 1h     |
|------------------------------------|----------------------|-----------------------|
| aia.lev1_euv_12s[94]               | aia.lev1_uv_24s[1600] | aia.lev1_vis_1h[4500] |
| aia.lev1_euv_12s[131]              | aia.lev1_uv_24s[1700] |                       |
| aia.lev1_euv_12s[193]              |                       |                       |
| aia.lev1_euv_12s[211]              |                       |                       |
| aia.lev1_euv_12s[304]              |                       |                       |
| aia.lev1_euv_12s[335]              |                       |                       |

---

## Output Structure

SpotiPy automatically organizes all outputs under `results_NOAA_{noaa}/`:

| Folder / File             | Description                                                           |
|---------------------------|-----------------------------------------------------------------------|
| `AIA_Aligned/`            | AIA images reprojected onto the HMI grid                             |
| `AIA_NoLimbDark/`         | AIA images with limb darkening removed                               |
| `Masks_HMI/`              | HMI umbra/penumbra segmentation masks (FITS) and overlay PNGs        |
| `Masks_AIA/`              | AIA segmentation maps (FITS) and overlay PNGs                        |
| `CLV_Plots_HMI/`          | CLV scatter + fit plots for each HMI observable (PNG)                |
| `CLV_Plots_AIA/`          | CLV scatter + fit plots for AIA masks per HMI observable (PNG)       |
| `Results_HMI/{obs}/`      | Extracted mu and intensity arrays per category (.txt)                |
| `Results_AIA/{obs}/`      | Extracted mu and intensity arrays per AIA category (.txt)            |
| `spotipy_NOAA_{noaa}.npz` | Consolidated archive of all frames, masks, timestamps, and coordinates |

The `.npz` archive contains the following arrays:

| Key           | Description                                      |
|---------------|--------------------------------------------------|
| `hmi_context` | All HMI cropped frames, shape `(N, F, F)`        |
| `hmi_seg`     | HMI segmentation masks, shape `(N, F, F)`        |
| `aia_context` | All AIA cropped frames, shape `(N, F, F)`        |
| `aia_seg`     | AIA segmentation maps, shape `(N, F, F)`         |
| `times`       | ISO timestamps for each frame                    |
| `x_arcsec`    | Helioprojective x-coordinate of crop center      |
| `y_arcsec`    | Helioprojective y-coordinate of crop center      |
| `longitude`   | Heliographic Stonyhurst longitude in degrees     |
| `latitude`    | Heliographic Stonyhurst latitude in degrees      |

where `N` is the number of timesteps and `F` is `FRAME_SIZE`.

---

## Citation

If you use this software in your research, please cite the associated JOSS paper:

> Sumra, M.S. & Lößnitz, E.J. (2026). *SpotiPy: A Modular Python Pipeline for Solar Spot Analysis*. Journal of Open Source Software, X(X), XXXX.

---

## License

This project is licensed under the **MIT License**. See the `LICENSE` file for details.
