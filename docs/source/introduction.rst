Introduction
============

**SpotiPy** is a modular Python framework designed for the analysis of solar active region feature identification and tracking in full-disk observations. While foundational libraries like ``astropy`` and ``sunpy`` provide core coordinate transformations and data access, **SpotiPy** addresses the practical gap in providing an integrated, multi-instrument feature-tracking pipeline.

The framework is optimized for data from the Helioseismic and Magnetic Imager (HMI) and the Atmospheric Imaging Assembly (AIA) aboard the Solar Dynamics Observatory (SDO). It supports coordinated analysis across five primary observables:

1. **Intensity ($I_c$)** — SDO/HMI Continuum
2. **Magnetogram ($M$)** — SDO/HMI LOS Magnetic Field
3. **Dopplergram ($V$)** — SDO/HMI Velocity
4. **Line Depth ($L_d$)** — SDO/HMI
5. **Line Width ($L_w$)** — SDO/HMI

Core Functionality
------------------
SpotiPy is built to be imported as standalone components or run as a unified pipeline. The core modules include:

* ``spotipy.downloading``: Batch data retrieval using ``sunpy.net.Fido`` with automatic JSOC authentication.
* ``spotipy.aligning``: Sub-pixel co-alignment of AIA images to the HMI grid using WCS reprojection.
* ``spotipy.limbdarkening_removal``: Radial profile-based limb darkening correction for AIA data.
* ``spotipy.segmentation``: Intensity-based masking of umbra, penumbra, plage, and network features.
* ``spotipy.tracking``: Differential rotation modeling, active feedback trajectory correction, and tracking.
* ``spotipy.gui_tools``: Interactive GUI for the manual selection of solar feature starting coordinates.

Available Data Series
---------------------
The downloading module enables seamless batch retrieval for a variety of JSOC data products.

**HMI Observables**
Measurements of the Fe I absorption line at 6173 Å.

.. list-table::
   :widths: 25 20 20 20 15
   :header-rows: 1

   * - Continuum Intensity (Ic)
     - Linewidth (Lw)
     - Linedepth (Ld)
     - Magnetograms (M)
     - Dopplergrams (V)
   * - hmi.Ic_45s
     - hmi.Lw_45s
     - hmi.Ld_45s
     - hmi.M_45s
     - hmi.V_45s
   * - hmi.Ic_720s
     - hmi.Lw_720s
     - hmi.Ld_720s
     - hmi.M_720s
     - hmi.V_720s
   * - hmi.Ic_noLimbDark_720s
     -
     -
     -
     -
   * - hmi.Ic_noLimbDark_720s_nrt
     -
     -
     -
     -

**AIA Observables**
Level 1 observables covering extreme UV, UV, and visible channels.

.. list-table::
   :widths: 33 33 34
   :header-rows: 1

   * - Extreme UV (euv), 12s
     - UV (uv), 24s
     - Visible (vis), 1h
   * - aia.lev1_euv_12s[94]
     - aia.lev1_uv_24s[1600]
     - aia.lev1_vis_1h[4500]
   * - aia.lev1_euv_12s[131]
     - aia.lev1_uv_24s[1700]
     -
   * - aia.lev1_euv_12s[193]
     -
     -
   * - aia.lev1_euv_12s[211]
     -
     -
   * - aia.lev1_euv_12s[304]
     -
     -
   * - aia.lev1_euv_12s[335]
     -
     -
