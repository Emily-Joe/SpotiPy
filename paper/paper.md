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
    orcid: 0009-0008-9890-4230
    affiliation: 1, 2
  - name: Emily Joe Lößnitz
    orcid: 0009-0006-3049-4875
    affiliation: 1, 2
  - name: Alexander G.M. Pietrow
    orcid: 0000-0002-0484-7634
    affiliation: 1
affiliations:
 - name: Leibniz Institute for Astrophysics Potsdam (AIP), Germany
   index: 1
 - name: University of Potsdam, Germany
   index: 2
date: 15 January 2026
bibliography: paper.bib
---

# Summary

Solar activity manifests itself in many different ways, from the near-omnipresent
network to the more confined active regions, which appear as groups of faculae,
spots, and pores in the photosphere, and as plage, filaments, and prominences in
the chromosphere. Our ever-improving physical understanding and astronomical
instrumentation have made these features the subject of study for several
centuries [@solanki2003; @Cretignier2024; @palumbo2024].

In recent years, more focus has been put on the consistent tracking and
observations of active regions with space-based facilities like the Solar
Dynamics Observatory (SDO) [@pesnell2012], the Chinese H$\alpha$ Solar
Explorer (CHASE) [@li2022], and the Interface Region Imaging Spectrograph
(IRIS) [@depontieu2014], but also in high resolution from the ground with
facilities such as the Swedish Solar Telescope (SST) [@scharmer2003] and
GREGOR [@schmidt2012], resulting in strong progression of our understanding
of these features [@pietrow2023; @morosin2022; @loessnitz2025]].

`SpotiPy` is a Python package developed for feature identification and tracking
in full-disk solar observations [@kontogiannis2024]. The framework has been
primarily designed for use with data from the Helioseismic and Magnetic Imager
(HMI) [@scherrer2012] and the Atmospheric Imaging Assembly (AIA) [@lemen2012]
aboard SDO, although the methods are not instrument-specific. The pipeline
enables the identification and tracking of arbitrary fields of view over time.
While initial applications have focused on sunspot studies, the implemented
procedures are equally valid for other types of solar structures, such as plage,
network, and filaments.

# Statement of Need

While general-purpose solar physics libraries such as `SunPy` [@sunpy_community2020]
provide functionality for data access, coordinate system handling, and file I/O,
they do not aim to offer an integrated pipeline for the tracking of solar active
regions. In particular, tasks such as time-series download, cross-instrument
alignment between observations, and consistent image preprocessing are typically
handled through custom scripts by individual researchers, which can be difficult
to reproduce and maintain across projects. While foundational libraries like astropy
[@astropy2022] handle the underlying coordinate transformations, and instrument-specific 
packages like aiapy [@aiapy2020] provide excellent calibration routines for SDO/AIA data, 
they lack the specific, multi-instrument feature-tracking integration that `SpotiPy` provides

`SpotiPy` addresses this practical gap by providing a modular framework that
supports coordinated retrieval of large datasets from HMI and AIA, sub-pixel
alignment of multi-instrument observations, and preprocessing steps such as limb
darkening correction. While the package was designed with feature tracing and
active region analysis in mind, its modular components can be applied to other
full-disk studies or as standalone functions within existing pipelines, enabling
reuse and sustained development.

# Functionality

`SpotiPy` is organized into five independent modules — (1) Downloading,
(2) Aligning, (3) Limb Darkening Correction, (4) Segmentation, and (5) Tracking —
that can be used together as a complete pipeline or individually within existing
workflows.

## 1. Downloading

`SpotiPy` provides a wrapper around `Fido` [@sunpy_community2020] that handles
JSOC email authentication and batch retrieval of HMI and AIA time series.
Because the module queries JSOC directly, it can retrieve any available data
series. The most common SDO data products used for feature tracking are
listed in Table 1.

| Instrument | Observable | Supported JSOC Series |
| :--- | :--- | :--- |
| HMI | Continuum Intensity | `hmi.Ic_45s`, `hmi.Ic_720s` |
| HMI | Continuum Intensity (Flattened) | `hmi.Ic_nolimbDark_720s`, `hmi.Ic_nolimbDark_720s_nrt` |
| HMI | Linewidth | `hmi.Lw_45s`, `hmi.Lw_720s` |
| HMI | Linedepth | `hmi.Ld_45s`, `hmi.Ld_720s` |
| HMI | LOS Magnetograms | `hmi.M_45s`, `hmi.M_720s` |
| HMI | Dopplergrams | `hmi.V_45s`, `hmi.V_720s` |
| AIA | Visible (4500 \AA) | `aia.lev1_vis_1h` |
| AIA | Ultraviolet (1600, 1700 \AA) | `aia.lev1_uv_24s` |
| AIA | Extreme UV (94, 131, 193, 211, 304, 335 \AA) | `aia.lev1_euv_12s` |

Table 1: HMI and AIA full-disk data products and passbands accessible via the downloading module.

To retrieve a time series, users must configure the following core parameters:

* **`series`**: the HMI or AIA data product identifier (e.g., `"hmi.Ic_720s"`)
* **`start_date`**: start of the time series in ISO 8601 format (e.g., `"2019-04-07T18:00:00"`)
* **`days`**: duration of the time series in days
* **`cadence_h`**: time step between downloaded frames in hours
* **`email`**: a JSOC-registered email address required for data export

While `download_series` is a general-purpose wrapper capable of retrieving any
JSOC data series, the following example demonstrates fetching the specific HMI
continuum and AIA 1700 \AA data required for the active region pipeline:

```python
# Download HMI continuum intensity and AIA 1700 A time series
from spotipy.downloading import download_series

hmi_files = download_series(
    series="hmi.Ic_720s",
    start_date="2019-04-07T18:00:00",
    days=11,
    cadence_h=6,
    email="user@institute.edu",
    out_dir="data/HMI"
)

aia_files = download_series(
    series="aia.lev1_uv_24s",
    start_date="2019-04-07T18:00:00",
    days=11,
    cadence_h=6,
    email="user@institute.edu",
    out_dir="data/AIA",
    wavelength_angstrom=1700  # Required to avoid downloading all 10 AIA channels
)
```

## 2. Aligning

`SpotiPy` automates the co-alignment of multi-instrument images using World
Coordinate Systems (WCS) and fast interpolation via the `reproject` package 
[@reproject]. The alignment module automatically handles the inconsistent
FITS file structures (e.g., data stored in primary versus secondary HDUs) across 
different SDO data products.This ensures that HMI and AIA observations share a common
pixel grid, as shown in \autoref{fig:alignment}.

```python
# Reproject AIA onto the HMI pixel grid and orientation
from spotipy.aligning import align_images

align_images(aia_path, hmi_reference_path, output_path)
```

![Reprojection of AIA 1700 \AA  onto the HMI pixel grid. Left: original AIA grid.
Center: HMI reference. Right: AIA reprojected to HMI
grid.\label{fig:alignment}](fig2_alignment.pdf){ width=100% }

## 3. Limb Darkening Correction

This function implements radial-profile subtraction to remove the center-to-limb
intensity gradient, which is essential for reliable intensity thresholding across
the full solar disk. Limb darkening is a well-known optical effect in which the
solar disk appears brighter at its center than at its edges, due to the viewing
geometry through different layers of the solar atmosphere. Removing this gradient
is therefore necessary before applying uniform intensity thresholds across the disk.

The correction is illustrated in \autoref{fig:limbdark} and \autoref{fig:radialprofile},
where a radial profile is generated from the data using the header values, and then
a limb is defined either based on the header or the maximum derivative taken around
the limb. After removal of this profile a mostly flat disk is left, occasionally
with artifacts due to large active regions (see \autoref{fig:radialprofile}).
However, these have not been shown to affect thresholding methods.

```python
# Remove limb darkening using radial profile subtraction
from spotipy.limbdarkening_removal import remove_limb_darkening, get_header_geometry

(cx, cy), r_pix = get_header_geometry(header)
aia_corrected = remove_limb_darkening(aia_raw, center=(cx, cy), radius_pix=r_pix)
```

![AIA 1700 \AA before (left) and after (right) limb darkening
correction.\label{fig:limbdark}](fig1_limb_darkening_aia.pdf){ width=100% }

![Radial intensity profile of the solar disk (left) and its derivative (right).
The dashed vertical line marks the detected limb radius (1567 pixels). The
intensity profile shows the characteristic center-to-limb gradient, which is
subtracted during correction. The derivative profile shows the sharp transition
at the limb boundary, used to determine the limb position when header values
are unavailable.\label{fig:radialprofile}](fig1b_radial_profile.pdf){ width=100% }

## 4. Segmentation

`SpotiPy` uses intensity thresholding combined with morphological operations to
generate binary masks for umbra, penumbra, and extended spot regions from HMI
continuum intensity, and for plage, network, and quiet Sun regions from AIA
1700 \AA. The thresholding approach identifies pixels within defined intensity
ranges normalized to the quiet Sun level. To optimize for OpenCV morphological 
operations, the HMI continuum arrays are mapped to an 8-bit scale `(0-255)` where 
the quiet Sun median corresponds to a pixel value of approximately 128. Umbral pixels
are identified using an 8-bit intensity range of 10 to 55 (corresponding to ~8% to 43% 
of the quiet Sun level), and penumbral pixels use a range of 75 to 120 (corresponding to 
~59% to 94% of the quiet Sun level). In contrast, plage regions in AIA are identified directly 
from the normalized data as pixels exceeding the quiet Sun level by more than 20% (`plage_excess_pct=20.0`) 
over an area of at least 600 pixels.


The segmentation masks can be applied across all five HMI observables
simultaneously — continuum intensity (`hmi.Ic_720s`), Doppler velocity
(`hmi.V_720s`), magnetic flux density (`hmi.M_720s`), line depth (`hmi.Ld_720s`),
and line width (`hmi.Lw_720s`) — enabling studies of center-to-limb variation
(CLV) within individual solar features.

\autoref{fig:fulldisk} shows the full-disk segmentation map, while
\autoref{fig:spotcrop} shows a zoomed view of an individual sunspot with
umbral and penumbral masks.

```python
# Generate segmentation masks using intensity thresholding
from spotipy.segmentation import get_masks, get_aia_masks

hmi_masks = get_masks(
    hmi_norm,                 # normalized HMI continuum intensity map
    disk_mask=ondisk,         # boolean mask selecting the solar disk
    umbra_range=(10, 55),     # umbra: 8-bit scale (~8-43% of quiet Sun)
    penumbra_range=(75, 120), # penumbra: 8-bit scale (~59-94% of quiet Sun)
    cleanup=True              # applies connected-component filter
)

aia_masks = get_aia_masks(
    aia_flat,                            # limb-darkening-corrected AIA image
    spot_mask=hmi_masks["spot"],         # main sunspot mask to exclude
    all_dark_mask=hmi_masks["all_dark"], # includes smaller pores for exclusion
    plage_excess_pct=20.0,
    qs_tol_pct=15.0,
    solar_center=(cx, cy),               # required for polar-slice QS reference
    r_pix=r_pix,
    full_disk_aia=aia_flat
)
```

![Full-disk segmentation map derived from HMI continuum and AIA 1700 \AA. Left:
HMI continuum intensity. Center: AIA 1700 \AA. Right: segmentation map showing
quiet Sun (gray), network (blue), plage (gold), and sunspot (red). Excluded
regions such as pores and the limb are shown in
black.\label{fig:fulldisk}](fig3_segmentation_20240522_181700.pdf){ width=100% }

![Sunspot segmentation of NOAA 12738 at disk center passage. (a) HMI continuum
crop. (b) HMI sunspot segmentation with umbra (red) and penumbra (orange).
(c) AIA 1700 \AA crop. (d) Combined segmentation map showing all region
types.\label{fig:spotcrop}](fig4_spot_segmentation.pdf){ width=100% }

## 5. Tracking

`SpotiPy` computes the expected pixel position of a solar feature over time 
using the differential rotation model derived by [@loessnitz2025]. This enables
consistent cropping of a tracked region across a multi-day time series. The 
tracking module first applies this physical model to predict where a feature at 
a given heliographic position will appear in subsequent frames. To account for 
local proper motions or small model deviations, the pipeline utilizes the 
`refine_centering` function to "lock on" to the darkest pixels of the sunspot umbra.
This ensures the feature remains centered within the extraction window throughout its transit.

The resulting coordinates are returned as NumPy arrays containing the refined heliographic and 
arcsecond positions for each time step. These data are used to generate the extraction crops and 
segmentation maps shown in \autoref{fig:strip}, which displays a time-summed strip of NOAA 12738 as 
it travels accross the solar disk. By calculating these tracks dynamically from the FITS headers, 
the pipeline ensures reproducibility without relying on intermediate external data files.

```python
# Track a solar feature and generate a time-summed visual strip
from spotipy.tracking import track_spots, refine_centering, strip

# 1. Calculate physical tracks from an initial position
tracks = track_spots(hmi_files, hmi_dir, x_start, y_start)

# 2. Refine tracks frame-by-frame to center on the umbra
refined_tracks = []
for i, frame in enumerate(hmi_files):
    # (Internal logic: Load frame, convert arcsec to pixels)
    # new_cx, new_cy = refine_centering(data_crop, rx0, ry0, frame_size)
    # refined_tracks.append(convert_to_arcsec(new_cx, new_cy))
    pass

# 3. Generate a time-summed strip using the refined coordinates
strip_img, _ = strip(
    series=hmi_files,
    directory=hmi_dir,
    tracks=refined_tracks,
    strip_height_arcsec=400.0,
    frame_size=400,
    overlay=True,
    animate=False
)
```

![Time-summed strip of NOAA 12738 travlling accross the solar disk, with tracked
bounding boxes overlaid at each time step. The color gradient from purple to
yellow indicates the progression of
time.\label{fig:strip}](time_summed_strip.pdf){ width=100% }

# Acknowledgements

We acknowledge the use of data from the Solar Dynamics Observatory, provided
by NASA. HMI and AIA data were accessed via the JSOC data export service.

<!-- TODO: Emily needs to add contribution note (rotation rate work) here -->

