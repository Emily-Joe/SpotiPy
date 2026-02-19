"""
Segmentation Module
-------------------
Tools for masking solar features.
Uses standard intensity thresholding to isolate sunspots (umbra/penumbra) in HMI continuum images,
and Quiet Sun mdian thresholding to find plage and network in AIA 1700Å.
"""

import numpy as np
import cv2
from scipy.ndimage import label
from typing import Dict, Tuple, Optional

def get_masks(
    image: np.ndarray,
    disk_mask: Optional[np.ndarray] = None,
    umbra_range: Tuple[float, float] = (10, 55),
    penumbra_range: Tuple[float, float] = (75, 120),
    cleanup: bool = True
) -> Dict[str, np.ndarray]:
    """
    Generates segmentation masks for Umbra, Penumbra, and the full Active Region.

    Utilizes standard MDI/HMI continuum intensity pipelines. The input array is
    normalized to an 8-bit scale for optimized OpenCV morphological thresholding.

    Parameters
    ----------
    image : np.ndarray
        2D array of solar disk intensity (typically normalized HMI continuum).
    disk_mask : np.ndarray, optional
        Boolean mask of valid solar disk pixels to exclude off-limb space.
    umbra_range : tuple
        Intensity bounds (min, max) for the umbra on a 0-255 scale.
    penumbra_range : tuple
        Intensity bounds (min, max) for the penumbra on a 0-255 scale.
    cleanup : bool
        If True, applies a connected-component filter to remove isolated noise.

    Returns
    -------
    dict
        Boolean arrays corresponding to 'umbra', 'penumbra', 'spot', and 'full_spot'.
    """
    # Normalize to 8-bit (0-255) for OpenCV compatibility
    vmin, vmax = 0.0, 2.0
    arr = np.clip(image.astype(float), vmin, vmax)
    norm = (arr - vmin) / (max(vmax - vmin, 1e-6))
    img_u8 = (norm * 255.0).astype(np.uint8)

    # Apply Gaussian blur to suppress granular noise and stabilize thresholding
    g_blur = cv2.GaussianBlur(img_u8, (7, 7), 0)

    # 1. Isolate Umbra
    u_raw = cv2.inRange(g_blur, int(umbra_range[0]), int(umbra_range[1]))
    u = cv2.erode(u_raw, np.ones((3, 3), np.uint8), iterations=1)

    # 2. Isolate Penumbra (Bandpass excluding the dilated umbra to ensure separation)
    p_band = cv2.inRange(g_blur, int(penumbra_range[0]), int(penumbra_range[1]))
    p = cv2.bitwise_and(p_band, cv2.bitwise_not(cv2.dilate(u, np.ones((7, 7), np.uint8))))

    # 3. Define the broader Active Region (Loose Spot Mask)
    f_loose = cv2.morphologyEx(
        cv2.inRange(g_blur, 0, int(penumbra_range[1])),
        cv2.MORPH_CLOSE,
        np.ones((5, 5), np.uint8)
    )
    f_spot = cv2.bitwise_or(f_loose, u)

    # Apply global disk constraints if provided
    if disk_mask is not None:
        dm = disk_mask.astype(bool)
        u_bool, p_bool, f_bool = (u > 0) & dm, (p > 0) & dm, (f_spot > 0) & dm
    else:
        u_bool, p_bool, f_bool = (u > 0), (p > 0), (f_spot > 0)

    # Filter out spurious magnetic artifacts leaving only the primary feature
    if cleanup:
        u_bool = _keep_largest(u_bool)
        p_bool = _keep_largest(p_bool)
        f_bool = _keep_largest(f_bool)

    return {"umbra": u_bool, "penumbra": p_bool, "spot": f_bool, "full_spot": f_bool}

def get_aia_masks(
    crop_aia: np.ndarray,
    spot_mask: np.ndarray,
    plage_excess_pct: float = 20.0,
    qs_tol_pct: float = 15.0,
    min_area: int = 450
) -> Dict[str, np.ndarray]:
    """
    Segments AIA 1700 Å UV data into Plage, Network, and Quiet Sun.

    Establishes a statistical Quiet Sun baseline by explicitly excluding
    sunspot pixels, preventing the dark umbra/penumbra from skewing the median.
    Applies configurable percentage-based emission thresholds.

    Parameters
    ----------
    crop_aia : np.ndarray
        Aligned and cropped 2D array of the AIA 1700 Å observation.
    spot_mask : np.ndarray
        Boolean mask of the active sunspot (from HMI) to be explicitly ignored.
    plage_excess_pct : float
        Percentage emission excess above the Quiet Sun median to classify as Plage.
    qs_tol_pct : float
        Percentage emission excess above the Quiet Sun median to classify as Network.
    min_area : int
        Minimum contiguous pixel area required to classify a feature as Plage.

    Returns
    -------
    dict
        Boolean arrays corresponding to 'plage', 'network', and 'qs' (Quiet Sun).
    """
    # Calculate Quiet Sun median strictly outside the sunspot boundary
    qs_median = np.nanmedian(crop_aia[~spot_mask])

    # Calculate absolute emission thresholds
    plage_thresh = qs_median * (1.0 + plage_excess_pct / 100.0)
    net_thresh = qs_median * (1.0 + qs_tol_pct / 100.0)

    # 1. Isolate Plage candidates (highest emission, excluding sunspot)
    plage_candidates = (crop_aia > plage_thresh) & (~spot_mask)

    # Apply spatial filtering to eliminate small-scale UV transient noise
    plage = np.zeros_like(plage_candidates, dtype=bool)
    if np.any(plage_candidates):
        lbl, n = label(plage_candidates)
        for i in range(1, n + 1):
            region = (lbl == i)
            if np.sum(region) >= min_area:
                plage[region] = True

    # 2. Isolate Magnetic Network (emission between QS tolerance and Plage)
    network = (crop_aia > net_thresh) & (~plage) & (~spot_mask)

    # 3. Define the pure Quiet Sun (baseline emission, non-magnetic)
    qs = (crop_aia <= net_thresh) & (~spot_mask)

    return {'plage': plage, 'network': network, 'qs': qs}

def _keep_largest(mask: np.ndarray) -> np.ndarray:
    """
    Isolates the largest contiguous region in a binary mask.
    Used to suppress minor magnetic/intensity artifacts outside the primary active region.
    """
    if mask is None or not np.any(mask):
        return mask
    lbl, n = label(mask)
    if n <= 1:
        return mask
    sizes = np.bincount(lbl.ravel())
    sizes[0] = 0  # Ignore background
    return (lbl == np.argmax(sizes))
