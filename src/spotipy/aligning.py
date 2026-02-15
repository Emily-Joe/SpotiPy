"""
Alignment Module
----------------
This module handles the co-alignment of solar images (e.g., AIA to HMI)
using WCS reprojection.
"""

import warnings
from astropy.wcs import WCS
from astropy.io import fits
from reproject import reproject_interp
import numpy as np

# RENAMED: align_aia_to_hmi -> align_images (Matches __init__.py)
def align_images(aia_path, hmi_path, output_path=None):
    """
    Reprojects an AIA image to match the geometry of an HMI image.

    Parameters
    ----------
    aia_path : str
        Path to the source AIA FITS file.
    hmi_path : str
        Path to the reference HMI FITS file.
    output_path : str, optional
        Path to save the aligned FITS file.

    Returns
    -------
    aligned_data : np.ndarray
        The AIA data reprojected to the HMI grid.
    header : astropy.io.fits.Header
        The new header for the aligned image.
    """

    with fits.open(aia_path) as hdul_a, fits.open(hmi_path) as hdul_h:
        data_a = hdul_a[1].data if len(hdul_a) > 1 and hdul_a[1].data is not None else hdul_a[0].data
        header_a = hdul_a[1].header if len(hdul_a) > 1 and hdul_a[1].data is not None else hdul_a[0].header

        data_h = hdul_h[1].data if len(hdul_h) > 1 and hdul_h[1].data is not None else hdul_h[0].data
        header_h = hdul_h[1].header if len(hdul_h) > 1 and hdul_h[1].data is not None else hdul_h[0].header

    wcs_a = WCS(header_a)
    wcs_h = WCS(header_h)

    print(f"Reprojecting AIA to match HMI grid ({data_h.shape})...")

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        aligned_data, _ = reproject_interp((data_a, wcs_a), wcs_h, shape_out=data_h.shape)

    new_header = header_h.copy()
    new_header['HISTORY'] = 'Aligned to HMI using SpotiPy'

    if output_path:
        fits.writeto(output_path, aligned_data, new_header, overwrite=True)
        print(f"Saved aligned image to: {output_path}")

    return aligned_data, new_header
