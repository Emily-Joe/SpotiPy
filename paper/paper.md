---
title: 'SpotiPy: A Modular Python Pipeline for Solar Spot Analysis and Tracking'
tags:
  - Python
  - astronomy
  - solar physics
  - sunspots
  - image processing
authors:
  - name: Saqib Sumra
    orcid: 0000-0000-0000-0000
    affiliation: 1, 2
  - name: Emily Joe
    orbit: 0000-0000-0000-0000
    affiliation: 1, 2
  - Alex Pietrow
    orcid: 0000-0000-0000-0000
    affiliation: 1
affiliations:
 - name: AIP
   index: 1
 - name: University of Potsdam, German 
   index: 2
date: 15 January 2026
bibliography: paper.bib
---

# Summary

Sunspots are fundamental features of the solar photosphere, manifesting as dark regions with intense magnetic fields. Analyzing their evolution is critical for understanding the solar cycle and space weather. `SpotiPy` is a Python package designed to automate the retrieval, processing, and analysis of sunspot data. It provides a modular pipeline that handles data downloading from public archives, image co-alignment, limb darkening removal, feature segmentation, and temporal tracking of sunspots.

# Statement of need

While general-purpose solar libraries like `SunPy` [@sunpy_community2020] provide excellent tools for handling coordinate systems and file I/O, there is a lack of high-level, cohesive tools specifically designed for sunspot lifecycle analysis. Researchers often write ad-hoc scripts to stitch together downloading, alignment, and segmentation steps, leading to non-reproducible workflows.

`SpotiPy` bridges this gap by providing a standardized, tested, and modular framework. It allows researchers to:
1.  **Download** large datasets from HMI and AIA instruments with minimal code.
2.  **Align** multi-instrument observations (e.g., UV and Visible light) with sub-pixel accuracy.
3.  **Process** images to remove effects like limb darkening.
4.  **Segment** features (umbra, penumbra) using physically motivated thresholds.
5.  **Track** these features over time to study their movement and evolution.

This modularity satisfies the need for "sustained development" and reusability, allowing users to incorporate specific components (e.g., just the segmentation algorithm) into their own pipelines.

# Key Features

`SpotiPy` is organized into five core modules:

* **Downloading:** A wrapper around `Fido` that handles JSOC email authentication and large-batch retrieval of HMI and AIA data.
* **Aligning:** Automates the co-alignment of images from different telescopes using World Coordinate Systems (WCS) and `reproject` [@reproject].
* **Processing:** Implements radial-profile subtraction to remove limb darkening, crucial for accurate intensity thresholding across the solar disk.
* **Segmentation:** Uses intensity thresholding combined with morphological operations to generate masks for the umbra, penumbra, and pores.
* **Tracking:** Calculates the expected position of solar features over time using differential rotation models, enabling consistent tracking of specific active regions.

# Usage Example

The following example demonstrates how to download data for a specific date and generate a sunspot mask:

```python
import spotipy
from spotipy.segmentation import get_masks

# 1. Download data (Feature 1)
files = spotipy.download_data("2014-05-20", duration_days=0.5, instrument='hmi')

# 2. Process the first image
# (Assuming file I/O handling here for brevity)
raw_image = ... 
clean_image = spotipy.remove_limb_darkening(raw_image, center=(2048, 2048), radius_pix=1900)

# 3. Generate Masks (Feature 4)
# Users can define their own physical thresholds
masks = get_masks(clean_image, umbra_range=(10, 55), penumbra_range=(75, 120))

print(f"Detected {masks['umbra'].sum()} pixels of Umbra.")
``` 

# Usage Example

 
