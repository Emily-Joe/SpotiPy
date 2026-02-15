"""
Segmentation Module
-------------------
This module provides tools for detecting and masking solar features
(spots, umbra, penumbra) using intensity thresholding and morphological operations.
"""

import numpy as np
import cv2
from scipy.ndimage import label

def get_masks(image, disk_mask=None, umbra_range=(10, 55), penumbra_range=(75, 120), cleanup=True):
    """
    Generates segmentation masks for Umbra, Penumbra, and the full Spot.

    This function mimics standard intensity thresholding workflows (like standard
    MDI/HMI pipelines) by converting the image to an 8-bit range and applying
    thresholds followed by morphological cleanup.

    Parameters
    ----------
    image : np.ndarray
        The 2D array containing the solar disk intensity data.
    disk_mask : np.ndarray, optional
        A boolean mask where True indicates valid solar disk pixels.
        Pixels outside this mask will be ignored.
    umbra_range : tuple, optional
        The (min, max) intensity values (0-255 scaled) for the Umbra.
        Default is (10, 55).
    penumbra_range : tuple, optional
        The (min, max) intensity values (0-255 scaled) for the Penumbra.
        Default is (75, 120).
    cleanup : bool, optional
        If True, keeps only the largest connected component (removes noise).
        Default is True.

    Returns
    -------
    masks : dict
        A dictionary containing the boolean masks:
        - 'umbra': Pixels identified as umbra.
        - 'penumbra': Pixels identified as penumbra.
        - 'spot': The combined feature mask (loose definition).
    """

    # 1. Normalize to 8-bit (0-255) for OpenCV
    # We clip values between 0.0 and 2.0 (standard solar norms)
    vmin, vmax = 0.0, 2.0
    arr = np.clip(image.astype(float), vmin, vmax)
    norm = (arr - vmin) / (max(vmax - vmin, 1e-6))
    img_u8 = (norm * 255.0).astype(np.uint8)

    # 2. Pre-processing (Blur to reduce noise)
    g_blur = cv2.GaussianBlur(img_u8, (7, 7), 0)

    # 3. Apply Thresholds
    # Umbra Mask
    u_raw = cv2.inRange(g_blur, int(umbra_range[0]), int(umbra_range[1]))
    u = cv2.erode(u_raw, np.ones((3, 3), np.uint8), iterations=1)

    # Penumbra Mask logic: Bandpass AND NOT Umbra
    p_band = cv2.inRange(g_blur, int(penumbra_range[0]), int(penumbra_range[1]))
    # Exclude the umbra from the penumbra (dilated slightly to ensure separation)
    p = cv2.bitwise_and(p_band, cv2.bitwise_not(cv2.dilate(u, np.ones((7, 7), np.uint8))))

    # Loose "Spot" Mask (Broad detection)
    f_loose = cv2.morphologyEx(
        cv2.inRange(g_blur, 0, int(penumbra_range[1])),
        cv2.MORPH_CLOSE,
        np.ones((5, 5), np.uint8)
    )
    f_spot = cv2.bitwise_or(f_loose, u)

    # 4. Apply Disk Mask (if provided)
    if disk_mask is not None:
        dm = disk_mask.astype(bool)
        u_bool = (u > 0) & dm
        p_bool = (p > 0) & dm
        f_bool = (f_spot > 0) & dm
    else:
        u_bool = (u > 0)
        p_bool = (p > 0)
        f_bool = (f_spot > 0)

    # 5. Cleanup (Keep only largest feature)
    if cleanup:
        u_bool = _keep_largest(u_bool)
        p_bool = _keep_largest(p_bool)
        f_bool = _keep_largest(f_bool)

    return {
        "umbra": u_bool,
        "penumbra": p_bool,
        "spot": f_bool
    }

def _keep_largest(mask):
    """Internal helper to keep only the largest connected component."""
    if mask is None or not np.any(mask):
        return mask
    lbl, n = label(mask)
    if n <= 1:
        return mask
    sizes = np.bincount(lbl.ravel())
    sizes[0] = 0  # Ignore background
    return (lbl == np.argmax(sizes))
