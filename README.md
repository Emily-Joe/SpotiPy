# SpotiPy: Solar Active Region Analysis Framework

[![JOSS Status](https://joss.theoj.org/papers/a6818174542289139f72787834526017/status.svg)](https://joss.theoj.org/papers/a6818174542289139f72787834526017)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**SpotiPy** is a modular Python framework for automated analysis of sunspot center-to-limb variations (CLV). It extends standard solar pipelines to support multi-wavelength and multi-observable analysis, with a focus on studying physical asymmetries between the **Leading (West)** and **Following (East)** hemispheres of solar active regions.

The package automatically downloads, aligns, segments, and analyzes data for **five solar observables**:

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

### Configuration File (`params.txt`)

The pipeline is controlled through a simple text-based configuration file defining the physical and observational parameters:

* **NOAA_NUMBER** — Active Region ID (e.g. `12673`)
* **START_DATE** — ISO timestamp (e.g. `2017-09-06T12:00:00`)
* **EMAIL** — Required for JSOC data access
* **CADENCE** — Time step in hours between observations

---

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

> Lößnitz, E.J., & Sumra, S. (2026). *SpotiPy: A Modular Python Pipeline for Solar Spot Analysis*. Journal of Open Source Software, X(X), XXXX.

---

## License

This project is licensed under the **MIT License**. See the `LICENSE` file for details.
