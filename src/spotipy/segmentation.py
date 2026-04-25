"""
Segmentation Module
-------------------
Tools for masking solar features.

Uses standard intensity thresholding to isolate sunspots (umbra/penumbra)
in HMI continuum images, and a polar-slice Quiet Sun median to find
plage and network in AIA 1700A.

Changes applied vs. Original version:

1. AIA QS reference uses a polar slice at disk centre (r=60-95%, angle
   70-110 deg) rather than the full crop median, preventing plage/network
   contamination of the threshold.
2. Full Spot mask uses MORPH_CLOSE before keep_largest to fill internal
   holes before component filtering.
3. QS mask is defined as pixels that are finite, above zero, not spot,
   not network, and not plage -- matching the original script definition
   and excluding dark calibration artefacts near limb.
4. get_masks() now returns 'all_dark': the raw pre-cleanup dark mask that
   includes small pores discarded by _keep_largest. Passed into
   get_aia_masks() so pores are never misclassified as Quiet Sun or
   Network in the AIA step -- matching dirty_mask logic in original script.

"""

import numpy as np
import cv2
from scipy.ndimage import label
from typing import Dict, Tuple, Optional

# =====================================================================
# HMI SEGMENTATION
# =====================================================================

def get_masks(
    image: np.ndarray,
    disk_mask: Optional[np.ndarray] = None,
    umbra_range: Tuple[float, float] = (10, 55),
    penumbra_range: Tuple[float, float] = (75, 120),
    cleanup: bool = True,
) -> Dict[str, np.ndarray]:
    """
    Generates segmentation masks for Umbra, Penumbra, and the full Active Region.

    Utilizes standard MDI/HMI continuum intensity pipelines. The input array is
    normalized to an 8-bit scale for optimized OpenCV morphological thresholding.

    Parameters
    ----------
    image : np.ndarray
        2D array of solar disk intensity (typically normalized HMI continuum,
        median ~ 1.0).
    disk_mask : np.ndarray, optional
        Boolean mask of valid solar disk pixels to exclude off-limb space.
    umbra_range : tuple
        Intensity bounds (min, max) for the umbra on a 0-255 scale.
    penumbra_range : tuple
        Intensity bounds (min, max) for the penumbra on a 0-255 scale.
    cleanup : bool
        If True, applies a connected-component filter to keep only the largest
        contiguous region in each mask, removing isolated noise pixels.

    Returns
    -------
    dict
        Boolean arrays for keys: 'umbra', 'penumbra', 'spot', 'all_dark'.
        'spot'     - largest contiguous active region (umbra + penumbra).
        'all_dark' - raw dark mask BEFORE _keep_largest, i.e. includes all small pores. Pass this to get_aia_masks() so pores are excluded from Quiet Sun and Network classification.
    """
    # Normalize to 8-bit (0-255) for OpenCV compatibility
    vmin, vmax = 0.0, 2.0
    arr = np.clip(image.astype(float), vmin, vmax)
    norm = (arr - vmin) / max(vmax - vmin, 1e-6)
    img_u8 = (norm * 255.0).astype(np.uint8)

    # Gaussian blur suppresses granular noise and stabilises thresholding
    g_blur = cv2.GaussianBlur(img_u8, (7, 7), 0)

    # 1. Umbra — threshold then erode once to remove single-pixel hits
    u_raw = cv2.inRange(g_blur, int(umbra_range[0]), int(umbra_range[1]))
    u = cv2.erode(u_raw, np.ones((3, 3), np.uint8), iterations=1)

    # 2. Penumbra — bandpass excluding the dilated umbra core
    p_band = cv2.inRange(g_blur, int(penumbra_range[0]), int(penumbra_range[1]))
    u_dil = cv2.dilate(u, np.ones((7, 7), np.uint8))
    p = cv2.bitwise_and(p_band, cv2.bitwise_not(u_dil))

    # 3. Full Spot — loose threshold closed morphologically to fill internal holes,
    #    then OR with the refined umbra to guarantee the core is always included.
    #    MORPH_CLOSE fills gaps before keep_largest runs, matching the original.
    f_loose = cv2.morphologyEx(
        cv2.inRange(g_blur, 0, int(penumbra_range[1])),
        cv2.MORPH_CLOSE,
        np.ones((5, 5), np.uint8),
        iterations=1,
    )
    f_spot = cv2.bitwise_or(f_loose, u)

    # Apply disk constraint
    if disk_mask is not None:
        dm = disk_mask.astype(bool)
        u_bool = (u > 0) & dm
        p_bool = (p > 0) & dm
        f_bool = (f_spot > 0) & dm
    else:
        u_bool = (u > 0)
        p_bool = (p > 0)
        f_bool = (f_spot > 0)

    # Save all_dark BEFORE _keep_largest strips the pores.
    # This is the equivalent of dirty_mask in the original reference script.
    all_dark = f_bool.copy()

    # keep_largest on umbra and penumbra individually (noise removal),
    # but apply it to Full Spot AFTER the morphological close so holes are
    # already filled — identical behaviour to the original clean_spot_mask().
    if cleanup:
        u_bool = _keep_largest(u_bool)
        p_bool = _keep_largest(p_bool)
        f_bool = _keep_largest(f_bool)

    return {"umbra": u_bool, "penumbra": p_bool, "spot": f_bool, "all_dark": all_dark}


# =====================================================================
# AIA SEGMENTATION
# =====================================================================

def _polar_slice_mask(
    shape: Tuple[int, int],
    center: Tuple[float, float],
    r_pix: float,
    r_min_frac: float = 0.60,
    r_max_frac: float = 0.95,
    th_start: float = 70.0,
    th_end: float = 110.0,
) -> np.ndarray:
    """
    Returns a boolean mask selecting a polar annular slice of the solar disk.

    Used to sample a clean Quiet Sun region away from the active region and
    the limb, matching the polar_slice_mask() approach in the original script.

    Parameters
    ----------
    shape : tuple
        (H, W) of the full-disk image.
    center : tuple
        (cx, cy) solar disk centre in pixels.
    r_pix : float
        Solar radius in pixels.
    r_min_frac, r_max_frac : float
        Radial range as fraction of r_pix (default 0.60 - 0.95).
    th_start, th_end : float
        Angular range in degrees, measured from the positive x-axis
        (default 70 - 110 deg corresponds roughly to the solar north polar cap).

    Returns
    -------
    np.ndarray
        Boolean mask of the same shape as the input image.
    """
    H, W = shape
    yy, xx = np.indices((H, W))
    X = xx - center[0]
    Y = yy - center[1]
    r_norm = np.sqrt(X**2 + Y**2) / float(r_pix)
    theta = (np.degrees(np.arctan2(Y, X)) + 360.0) % 360.0
    return (r_norm >= r_min_frac) & (r_norm <= r_max_frac) & (theta >= th_start) & (theta <= th_end)


def get_aia_masks(
    crop_aia: np.ndarray,
    spot_mask: np.ndarray,
    all_dark_mask: Optional[np.ndarray] = None,
    plage_excess_pct: float = 20.0,
    qs_tol_pct: float = 15.0,
    min_area: int = 450,
    blur_sigma: int = 3,
    # Compute QS median from polar slice
    solar_center: Optional[Tuple[float, float]] = None,
    r_pix: Optional[float] = None,
    full_disk_aia: Optional[np.ndarray] = None,
    ps_r_min: float = 0.60,
    ps_r_max: float = 0.95,
    ps_th_start: float = 70.0,
    ps_th_end: float = 110.0,
    # Crop origin within the full disk (needed for polar slice)
    x0: int = 0,
    y0: int = 0,
) -> Dict[str, np.ndarray]:
    """
    Segments AIA 1700 A UV data into Plage, Network, and Quiet Sun.

    Polar-slice Quiet Sun:
        When solar_center, r_pix, and full_disk_aia are provided the Quiet Sun
        median is calculated from a polar slice of the full-disk image
        (r = 60-95% of R_sun, angle 70-110 deg) that lies entirely outside the
        active region. This prevents plage/network pixels in the crop from
        contaminating the median, which would push thresholds up and cause
        under-detection of plage.
        If those arguments are not provided (e.g. in test mode) the function
        falls back to the crop median outside all dark pixels.

    Pore exclusion via all_dark_mask:
        Small pores are discarded by _keep_largest in get_masks() and therefore
        absent from spot_mask. Without an explicit exclusion they fall into the
        AIA intensity range of Quiet Sun, contaminating QS statistics.
        Passing all_dark_mask (the pre-cleanup HMI dark mask) ensures every
        HMI-dark pixel -- main spot AND pores -- is excluded from QS and Network,
        matching the dirty_mask logic of the original reference script.

    QS mask definition:
        QS is defined as pixels that are finite, above zero, not dark in HMI,
        not plage, and not network.

    Parameters
    ----------
    crop_aia : np.ndarray
        Aligned, limb-darkening-corrected, cropped AIA 1700 A image.
    spot_mask : np.ndarray
        Boolean mask of the main sunspot from HMI (same shape as crop_aia).
    all_dark_mask : np.ndarray, optional
        Raw HMI dark mask before _keep_largest filtering, i.e. includes pores.
        Returned as 'all_dark' by get_masks(). When provided, pores are excluded
        from Network and Quiet Sun classification.
    plage_excess_pct : float
        % excess above QS median to classify as Plage (default 20%).
    qs_tol_pct : float
        % excess above QS median to classify as Network (default 15%).
    min_area : int
        Minimum contiguous pixel area for a plage region (default 450 px).
    blur_sigma : int
        Gaussian blur kernel size for plage candidate smoothing (default 3).
    solar_center : tuple, optional
        (cx, cy) solar centre in the FULL-DISK image, from FITS header.
    r_pix : float, optional
        Solar radius in pixels, from FITS header.
    full_disk_aia : np.ndarray, optional
        Full-disk aligned AIA image before cropping, used for polar-slice median.
    ps_r_min, ps_r_max : float
        Radial range of the polar slice as fraction of R_sun.
    ps_th_start, ps_th_end : float
        Angular range of the polar slice in degrees.
    x0, y0 : int
        Top-left corner of the crop in the full-disk frame.
        Used to translate solar_center into the crop coordinate system.

    Returns
    -------
    dict
        Boolean arrays for keys: 'plage', 'network', 'qs'.

    """
    # ----------------------------------------------------------------
    # Build combined dark exclusion mask (main spot + pores)
    # ----------------------------------------------------------------
    if all_dark_mask is not None:
        dark_excl = spot_mask | all_dark_mask.astype(bool)
    else:
        dark_excl = spot_mask.astype(bool)

    # ----------------------------------------------------------------
    # Compute QS median from polar slice on the full-disk image
    # ----------------------------------------------------------------
    use_polar_slice = (
        solar_center is not None
        and r_pix is not None
        and full_disk_aia is not None
    )

    if use_polar_slice:
        ps_mask = _polar_slice_mask(
            full_disk_aia.shape,
            center    = solar_center,
            r_pix     = r_pix,
            r_min_frac= ps_r_min,
            r_max_frac= ps_r_max,
            th_start  = ps_th_start,
            th_end    = ps_th_end,
        )
        qs_vals = full_disk_aia[ps_mask]
        qs_vals = qs_vals[np.isfinite(qs_vals) & (qs_vals > 0)]
        qs_median = np.nanmedian(qs_vals) if qs_vals.size > 0 else np.nanmedian(crop_aia[~dark_excl])
    else:
        # use crop median outside all dark pixels (test mode / no header info)
        qs_median = np.nanmedian(crop_aia[~dark_excl])

    # ----------------------------------------------------------------
    # Emission thresholds
    # ----------------------------------------------------------------
    plage_thresh = qs_median * (1.0 + plage_excess_pct / 100.0)
    net_thresh   = qs_median * (1.0 + qs_tol_pct / 100.0)

    # ----------------------------------------------------------------
    # 1. Plage — blurred to reduce UV transient noise, then area filter
    # ----------------------------------------------------------------
    # Exclude all dark pixels (spot + pores) before blurring
    aia_nodark = np.where(dark_excl, np.nan, crop_aia.astype(float))
    blurred = cv2.GaussianBlur(
        np.nan_to_num(aia_nodark).astype(np.float32),
        (blur_sigma, blur_sigma), 0
    )

    plage_candidates = (blurred > plage_thresh) & (~dark_excl)

    plage = np.zeros_like(plage_candidates, dtype=bool)
    if np.any(plage_candidates):
        lbl, n = label(plage_candidates)
        for i in range(1, n + 1):
            region = (lbl == i)
            if np.sum(region) >= min_area:
                plage[region] = True

    # ----------------------------------------------------------------
    # 2. Network — emission between QS tolerance and plage threshold
    # ----------------------------------------------------------------
    network = (crop_aia > net_thresh) & (~plage) & (~dark_excl)

    # ----------------------------------------------------------------
    # QS — finite, positive, not dark (spot or pore), not plage, not network
    # ----------------------------------------------------------------
    valid = np.isfinite(crop_aia) & (crop_aia > 0)
    qs = valid & (~dark_excl) & (~plage) & (~network)

    return {"plage": plage, "network": network, "qs": qs}


# =====================================================================
# UTILITIES
# =====================================================================

def _keep_largest(mask: np.ndarray) -> np.ndarray:
    """
    Isolates the largest contiguous region in a binary mask.
    Suppresses isolated noise pixels outside the primary active region.
    """
    if mask is None or not np.any(mask):
        return mask
    lbl, n = label(mask)
    if n <= 1:
        return mask
    sizes = np.bincount(lbl.ravel())
    sizes[0] = 0  # ignore background label
    return (lbl == np.argmax(sizes))
