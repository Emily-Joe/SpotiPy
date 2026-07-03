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
    corresponding: true
    orcid: 0000-0002-0484-7634
    affiliation: 1
  - name: Carsten Denker
    orcid: 0000-0002-7729-6415
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
spots, and pores in the photosphere, and as filaments and prominences in
the chromosphere. Our ever-improving physical understanding and astronomical
instrumentation has made these features the subject of study for several
centuries [@solanki2003; @Cretignier2024; @palumbo2024].

In recent years, more focus has been placed on the consistent tracking and
observations of active regions with space-based facilities such as the Solar
Dynamics Observatory [SDO, @pesnell2012], the Hinode (Solar-B) Mission [@kosugi2007], the Chinese H$\alpha$ Solar
Explorer [CHASE, @li2022], and the Interface Region Imaging Spectrograph
[IRIS, @depontieu2014], but also in high resolution from the ground with
facilities such as the Swedish Solar Telescope [SST, @scharmer2003] and the
GREGOR solar telescope [@schmidt2012], resulting in strong 
  progression of our understanding of these features.

`SpotiPy` is a Python package developed for feature identification and tracking
in full-disk solar observations. The framework has been
primarily designed for use with data from the Helioseismic and Magnetic Imager
[HMI, @scherrer2012] and the Atmospheric Imaging Assembly [AIA, @lemen2012]
aboard SDO, although the methods are not instrument-specific. The pipeline
enables the identification and tracking of arbitrary fields of view over time.
While initial applications have focused on sunspot studies, the implemented
procedures are equally valid for other types of solar features such as plage,
network, and filaments.

*Corresponding author: Alexander G.M. Pietrow ([apietrow@aip.de](mailto:apietrow@aip.de))*

# Statement of Need

While general-purpose solar physics libraries such as `SunPy` [@sunpy_community2020]
provide functionality for data access, coordinate system handling, and file I/O,
they do not aim to offer an integrated pipeline for the tracking of solar active
regions. In particular, tasks such as time-series download, cross-instrument
alignment between observations, and consistent image preprocessing are typically
handled through custom scripts by individual researchers, which can be difficult
to reproduce and maintain across projects. While foundational libraries such as `astropy`
[@astropy2022] handle the underlying coordinate transformations and instrument-specific 
packages such as `aiapy` [@aiapy2020] provide excellent calibration routines for SDO/AIA data, 
they lack the specific, multi-instrument feature-tracking integration that `SpotiPy` provides.

`SpotiPy` addresses this practical gap by providing a modular framework that
supports coordinated retrieval of large datasets from HMI and AIA, sub-pixel
alignment of multi-instrument observations, and preprocessing steps such as limb
darkening correction. While the package was designed with feature tracking and
active region analysis in mind, its modular components can be applied to other
full-disk studies or as standalone functions within existing pipelines, enabling
reuse and sustained development. `SpotiPy` was developed to facilitate high-throughput, 
reproducible analysis of solar features and has already been successfully employed to 
investigate the center-to-limb variations (CLV) of sunspots, faculae, and the solar network [@pietrow2026].

# Functionality

`SpotiPy` is organized into five independent modules — (1) Data download,
(2) Alignment, (3) Limb Darkening Correction, (4) Segmentation, and (5) Tracking —
that can be used together as a complete pipeline or individually within existing
workflows.

## 1. Data downloading

`SpotiPy` provides a wrapper around `Fido` that handles
JSOC email authentication and batch retrieval of HMI and AIA time series.
Because the module queries JSOC directly, it can retrieve any available data
series. The most common SDO data products used for feature tracking are
listed in Table 1.

Table 1: HMI and AIA full-disk data products and passbands accessible via the downloading module.

| Instrument | Observable | Supported JSOC Series |
| :--- | :--- | :--- |
| HMI | Continuum Intensity | `hmi.Ic_45s`, `hmi.Ic_720s` |
| HMI | Continuum Intensity (Flattened) | `hmi.Ic_nolimbDark_720s`, `hmi.Ic_nolimbDark_720s_nrt` |
| HMI | Line width | `hmi.Lw_45s`, `hmi.Lw_720s` |
| HMI | Line depth | `hmi.Ld_45s`, `hmi.Ld_720s` |
| HMI | LOS Magnetograms | `hmi.M_45s`, `hmi.M_720s` |
| HMI | Dopplergrams | `hmi.V_45s`, `hmi.V_720s` |
| AIA | Continuum Intensity (4500\ Å) | `aia.lev1_vis_1h` |
| AIA | Ultraviolet (1600 Å, 1700\ Å) | `aia.lev1_uv_24s` |
| AIA | Extreme UV (94\ Å, 131\ Å, 193\ Å, 211\ Å, 304\ Å, 335\ Å) | `aia.lev1_euv_12s` |


To retrieve a time series, users must configure the following core parameters:

* **`series`**: the HMI or AIA data product identifier (e.g., `"hmi.Ic_720s"`)
* **`start_date`**: start of the time series in ISO 8601 format (e.g., `"2019-04-07T18:00:00"`)
* **`days`**: duration of the time series in days
* **`cadence_h`**: time step between downloaded frames in hours
* **`email`**: a JSOC-registered email address required for data export

While `download_series` is a general-purpose wrapper capable of retrieving any
JSOC data series, the following example demonstrates fetching the specific HMI
continuum and AIA 1700\ Å data required for the active region pipeline:

```python
# Download HMI continuum intensity and AIA 1700 Å time series
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

## 2. Data Alignment

`SpotiPy` automates the co-alignment of multi-instrument images using World
Coordinate Systems (WCS) and fast interpolation via the `reproject` package 
[@reproject]. The alignment module automatically handles the different
FITS file structures (e.g., data stored in primary versus secondary header data units (HDUs)) across 
different SDO data products.This ensures that HMI and AIA observations share a common
pixel scale of 0.6 arcseconds per pixel, as shown in \autoref{fig:alignment}.

```python
# Reproject AIA onto the HMI pixel grid and orientation
from spotipy.aligning import align_images
align_images(aia_path, hmi_reference_path, output_path)
```

![Reprojection of AIA 1700\ Å intensity onto the HMI pixel grid. Left: original AIA grid.
Center: HMI reference. Right: AIA reprojected to HMI
grid.\label{fig:alignment}](fig2_alignment.pdf){ width=100% }

## 3. Limb Darkening Correction

This function implements radial-profile subtraction to remove the center-to-limb
intensity profile for HMI continuum and AIA 1600\ Å and 1700\ Å observations, 
which is essential for reliable intensity thresholding across
the full solar disk. Limb darkening is a radiative transfer effect induced by viewing geometry
that allows the observer to look deeper into the atmosphere at the center of the solar disc 
than at the limb, where the line of sight passes through the atmosphere at an angle. The decreasing 
temperature with height in the solar photosphere results in a darker limb and a 
brighter center [@Foukal2004; @Gray2005]. Removing this gradient is therefore necessary
before applying uniform intensity thresholds across the disk.

The correction is illustrated in \autoref{fig:limbdark} and \autoref{fig:radialprofile},
where a radial profile is generated and the limb radius is extracted directly from the FITS
header parameters. After removal of this profile a mostly flat disk is left, occasionally
with artifacts due to large active regions (see \autoref{fig:radialprofile}).
However, these have not been shown to affect thresholding methods.

```python
# Remove limb darkening using radial profile subtraction
from spotipy.limbdarkening_removal import remove_limb_darkening, get_header_geometry
(cx, cy), r_pix = get_header_geometry(header)
aia_corrected = remove_limb_darkening(aia_raw, center=(cx, cy), radius_pix=r_pix)
```

![AIA 1700\ Å intensity before (left) and after (right) limb darkening
correction.\label{fig:limbdark}](fig1_limb_darkening_aia.pdf){ width=100% }

![Radial intensity profile of the solar disk (left) and its derivative (right). 
The left panel displays the raw intensity profile overlaid with a 5th-order polynomial 
fit describing the CLV. The dashed vertical 
line marks the limb radius extracted from the FITS header parameters (1564 pixels). 
The derivative profile on the right illustrates the sharp intensity drop at the limb 
boundary.\label{fig:radialprofile}](fig1b_radial_profile.pdf){ width=100% }

## 4. Segmentation {#sec:segmentation}

`SpotiPy` uses intensity thresholding combined with morphological operations to
generate binary masks for umbra, penumbra, and extended spot regions from HMI
continuum intensity, and for plage, network, and quiet Sun regions from AIA
1700\ Å [@Shen2018; @Verma2018]. The segmentation is performed using intensity
thresholding combined with size-based  classification. These thresholds were 
empirically calibrated to $8–43\%$ of the normalized quiet-Sun intensity for the 
umbra, $59–94\%$ for the penumbra, $>15\%$ for the network, and $>20\%$ for plage. These 
values are broadly consistent with those reported in previous studies [@Chapman1989; @gyori1998; @preminger2001; @solanki2003; @kiess2014; @hoeksema2014]. 

To distinguish between network and plage, contiguous regions exceeding the network intensity 
threshold are further classified based on their area: structures larger than $450$ pixels are
identified as plage, while smaller structures are classified as network. We note that this introduces a center-to-limb inconsistency, as a fixed pixel threshold corresponds to a larger physical area at smaller $\mu$. However, the chosen value has been found to provide a robust performance across the disk.

The segmentation masks can be applied across all five HMI observables
simultaneously — continuum intensity (`hmi.Ic_720s`), Doppler velocity
(`hmi.V_720s`), magnetic flux density (`hmi.M_720s`), line depth (`hmi.Ld_720s`),
and line width (`hmi.Lw_720s`) — enabling studies of the CLV within individual solar features.

\autoref{fig:fulldisk} shows the full-disk segmentation map, while
\autoref{fig:spotcrop} shows a ROI of an individual sunspot with
umbral and penumbral masks.

```python
# Generate segmentation masks using intensity thresholding
from spotipy.segmentation import get_masks, get_aia_masks

hmi_masks = get_masks(
    hmi_norm,                 # normalized HMI continuum intensity map
    disk_mask=ondisk,         # boolean mask selecting the solar disk
    umbra_range=(10, 55),     # umbra: 8-bit scale (~$8-43\%$ of quiet Sun)
    penumbra_range=(75, 120), # penumbra: 8-bit scale (~$59-94\%$ of quiet Sun)
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

![Full-disk segmentation map derived from HMI continuum intensity and AIA 1700\ Å intensity. Left:
HMI continuum intensity. Center: AIA 1700\ Å. Right: segmentation map showing
quiet Sun (gray), network (blue), plage (gold), and sunspot (red). Excluded
regions such as pores and regions beyond the limb are shown in
black.\label{fig:fulldisk}](fig3_full_disc_segmentation.pdf){ width=100% }

![Sunspot segmentation of active region NOAA 12738 at central meridian passage. Left: HMI continuum ROI. 
Center: AIA 1700\ Å intensity ROI. Right: Combined segmentation map showing umbra (red), penumbra 
(orange), plage (gold), network (blue), and quiet Sun (gray). Small, isolated dark pores 
excluded (black) to ensure the analysis remains focused on the primary active region.
\label{fig:spotcrop}](fig4_spot_segmentation.pdf){ width=100% }



## 5. Tracking

The tracking module allows for consistent feature tracking at an arbitrary cadence, allowing for the study of a region and its parameters as it makes its way across the solar disk, and allows for the calculation of differential rotation rates. 

This module is based on the processes laid out in [@loessnitz2025]. A core functionality is based on a [@sunpy_community2020] model to compute the expected position of a solar feature over a time series for a given start position, which can be selected via an interactive clicker or by passing starting coordinates. The resulting coordinates are returned and allow for the isolation of an ROI by cropping out a window around it for each time step. To account for local motions and deviations from the rotation model, the pipeline provides an optional recentering function that refines the feature position centering the window around the darkest feature in the cropped image. This is done by smoothening the ROI with a Gaussian and then finding the center of mass of the largest feature. The extracted ROIs can be visualized in time-summed strips, as shown in Figure \autoref{fig:strip}. This series of crops can then be used for segmentation maps (#sec:segmentation) or further analysis.
  

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
![Time-summed strip of active region NOAA 12738 travelling accross the solar disk, with tracked
bounding boxes overlaid at each time step. The color of the bounding boxes, from purple to
yellow, indicates the progression of
time.\label{fig:strip}](time_summed_strip.pdf){ width=100% }


The determination of feature rotation rates can be performed using either the masks produced by the [Segmentation](#sec:segmentation) module or an internal masking process specifically designed for sunspots. Both approaches identify features such as sunspots at each time step. The centers of the sunspot masks are then determined either from image moments or, for approximately circular features, by fitting a minimal-area ellipse. The resulting positions are tracked over time and finally transformed into the heliographic coordinate system using observational parameters such as the solar radius, tilt, and distance obtained from the FITS header. The heliographic coordinate system expresses points on the solar surface by latitude and longitude [@Thompson2006]. This transformation is necessary to correct for projection effects caused by the spherical geometry, as illustrated in Figure \autoref{fig:rotation_rate}.  After deprojection, the tracked positions are used to calculate rotation rates, which are expressed in the Carrington reference frame [@Carrington1863]. This conversion follows the method described in [@Balthasar1979]. The Carrington System rotates as a rigid sphere along with the surface of the Sun, ensuring comparability between rotation studies based on SpotiPy and others. 

![Rotation rate of a single feature (a sunspot in active region NOAA 12738) during its passage across the solar disk. The left panel shows the calculated rotation rate without correcting for projection effects, while the right panel shows the deprojected rotation rate in the Carrington reference frame using the two available methods (image moments and ellipse fitting). Near the solar limb, the apparent motion of the feature corresponds to only a few pixels between consecutive frames. After deprojection, these small positional uncertainties translate into larger uncertainties in the calculated rotation rate. The breaks in the plotted curves occur where the Carrington longitude wraps from 360° to 0°. \label{fig:rotation_rate}](fig7_rotation_plot_comparison.pdf){ width=100% }


# Acknowledgements

We acknowledge the use of data from the Solar Dynamics Observatory, provided
by NASA. HMI and AIA data were accessed via the JSOC data export service. 
AP was supported by grant PI 2102/1-1 from the Deutsche Forschungsgemeinschaft (DFG). 
DeepL Write was used for assistance with proofreading and copy-editing.

# References
