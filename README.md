# SpotiPy: Solar Active Region Analysis Framework

[![JOSS Status](https://joss.theoj.org/papers/a6818174542289139f72787834526017/status.svg)](https://joss.theoj.org/papers/a6818174542289139f72787834526017)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**SpotiPy** is a modular Python framework for automated analysis of sunspot center-to-limb variations (CLV). It extends standard solar pipelines to support multi-wavelength and multi-observable analysis, with a focus on studying physical asymmetries between the **Leading (West)** and **Following (East)** hemispheres of solar active regions.

The package automatically downloads, aligns, segments, and analyzes data for five solar observables:

1. **Intensity ($I_c$)** — SDO/HMI Continuum
2. **Magnetogram ($M$)** — SDO/HMI LOS Magnetic Field
3. **Dopplergram ($V$)** — SDO/HMI Velocity
4. **Line Depth ($L_d$)** — SDO/HMI
5. **Line Width ($L_w$)** — SDO/HMI

---

## Key Features

SpotiPy is designed to be fully modular. You can run the complete pipeline end-to-end or import individual components into your own analysis scripts.

* **Multi-Observable Pipeline**
  Automated processing of co-aligned intensity, magnetic field, Doppler velocity, line depth, and line width data.

* **Spatial Asymmetry Analysis**
  Automatic separation of sunspot data into Leading and Following hemispheres to study rotational and evolutionary effects.

* **Statistical CLV “Candle” Plots**
  Box-plot–style distributions of physical quantities as a function of $\mu = \cos(\theta)$.

* **Modular Architecture**

  * `spotipy.downloading` — Batch data retrieval with automatic JSOC authentication
  * `spotipy.aligning` — Sub-pixel co-alignment using WCS reprojection
  * `spotipy.segmentation` — Intensity-based umbra and penumbra masking
  * `spotipy.tracking` — Differential rotation modeling and region tracking

---

## Installation

Clone the repository and install it in **editable mode**. This is recommended for research use, as changes to the code are applied immediately.

```bash
git clone https://github.com/Emily-Joe/SpotiPy.git
cd SpotiPy
pip install -e .
```

This will automatically install required dependencies, including `sunpy`, `astropy`, `scipy`, and `reproject`.

---

## Quick Start (Library Mode)

You can import SpotiPy tools directly into your own Python scripts:

```python
import spotipy
from spotipy.segmentation import get_masks
import numpy as np

# Example: download data (optional)
# files = spotipy.download_data("2014-05-20", duration_days=0.5)

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

For full end-to-end analysis (**Download → Align → Segment → Analyze → Plot**), use the example scripts provided in the `examples/` directory.

```bash
# Run the full analysis using a configuration file
python3 examples/run_analysis.py --config params.txt
```

## Configuration File (`params.txt`)

The pipeline is controlled through a simple text-based configuration file defining the physical and observational parameters:

* **NOAA_NUMBER** — Active Region ID (e.g. `12673`)
* **START_DATE** — ISO timestamp (e.g. `2017-09-06T12:00:00`)
* **EMAIL** — Required for JSOC data access
* **CADENCE** — Time step in hours between observations
---

## Avilable Data for Download

The downloading module enables the search and retrieval of data from public archives such as the Solar Dynamics Observatory’s HMI and AIA instruments. It provides a streamlined interface built on `sunpy`'s Fido client. Retrieval of HMI and AIA data requires prior registration via email with the JSOC data center hosting the SDO data. Visit [JSOC Email Registration](http://jsoc.stanford.edu/ajax/register_email.html) and register your email adress. This email address must also be specified in the configuration file `params.txt`. 

The following datasets are available for download: 

1. HMI Observables:
HMI generally provides data series with cadences of 45 or 720 seconds. It offers continuum intensity measurements in the region of the Fe I absorption line at 6173 Å, as well as a version with the effects of limb darkening removed. The limb-darkening–corrected dataset is only available at a cadence of 12 minutes. Additional observables include the line depth and line width of a Gaussian absorption line profile for the 6173.3 Å line, as well as line-of-sight (LOS) magnetograms and Dopplergrams.

| Continuum Intensity (Ic)   | Linewidth (Lw) | Linedepth (Ld) | LOS Magnetograms (M) | LOS Dopplergrams (V) |
| ---------------------------|----------------|----------------|----------------------|----------------------|
| hmi.Ic_45s                 | hmi.Lw_45s     | hmi.Ld_45s     | hmi.M_45s            | hmi.V_45s            |
| hmi.Ic_720s                | hmi.Lw_720s    | hmi.Ld_720s    | hmi.M_720s           | hmi.V_720s           |
| hmi.Ic_nolimbDark_720s     |                |                |                      |                      |
| hmi.Ic_nolimbDark_720s_nrt |                |                |                      |                      |

2. AIA Observables:
At Level 1, AIA observables cover three wavelength ranges (soft X-ray/extreme UV, UV, and visible), with ten total wavelength channels. The wavelengths are given in angstroms at the end of each series name. Each wavelength range is available at different cadences: every 12 seconds for EUV, every 24 seconds for UV, and once per hour for the visible range.

| Soft Xray - Extreme UV (euv), 12s | UV (uv), 24s            | visible (vis), 1h     |
|-----------------------------------|-------------------------|-----------------------|
| aia.lev1_euv_12s[94]              | aia.lev1_uv_12s[1600]   |aia.lev1_vis_1h[4500]  |
| aia.lev1_euv_12s[131]             | aia.lev1_uv_12s[1700]   |                       |
| aia.lev1_euv_12s[193]             |                         |                       |
| aia.lev1_euv_12s[211]             |                         |                       |
| aia.lev1_euv_12s[304]             |                         |                       |
| aia.lev1_euv_12s[335]             |                         |                       |


## Output Structure

SpotiPy automatically organizes outputs into a structured directory tree:

| Folder              | Description                                                     |
| ------------------- | --------------------------------------------------------------- |
| `FITS_files_.../`   | Raw downloaded FITS data (separated by instrument)              |
| `masks_.../`        | Generated umbra/penumbra segmentation masks (FITS format)       |
| `Results_.../`      | Extracted numerical data (text files containing $\mu, I, x, y$) |
| `Post_CLV_candles/` | Final CLV “Candle” plots (PNG)                                  |

---

## Citation

If you use this software in your research, please cite the associated JOSS paper:

> Sumra, S. & Lößnitz, E.J. (2026). *SpotiPy: A Modular Python Pipeline for Solar Spot Analysis*. Journal of Open Source Software, X(X), XXXX.

---

## License

This project is licensed under the **MIT License**. See the `LICENSE` file for details.
