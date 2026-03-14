"""
Alignment Module
----------------
This module handles the co-alignment of solar images (e.g., AIA to HMI)
using WCS reprojection. It takes the file paths and saving the aligned
output directly.
"""

import os
import warnings
from astropy.wcs import WCS
from astropy.io import fits
from reproject import reproject_interp

def read_fits(path):
    """To read SDO FITS files which may have data in index 0 or 1."""
    try:
        with fits.open(path) as hdul:
            if len(hdul) > 1 and hdul[1].data is not None:
                return hdul[1].data, hdul[1].header
            return hdul[0].data, hdul[0].header
    except Exception as e:
        print(f"Error reading {path}: {e}")
        return None, None

def align_images(aia_path, hmi_path, out_path):
    """
    Reprojects an AIA FITS file to match the geometry of an HMI FITS file,
    saving the result to disk.

    Parameters
    ----------
    aia_path : str
        Path to the source AIA FITS file.
    hmi_path : str
        Path to the reference HMI FITS file.
    out_path : str
        Path to save the aligned FITS file.

    Returns
    -------
    bool
        True if successful or file already exists, False if it failed.
    """

    # 1. Skip if already done
    if os.path.exists(out_path) and os.path.getsize(out_path) > 0:
        return True

    # 2. Read the files
    a_dat, a_hdr = read_fits(aia_path)
    h_dat, h_hdr = read_fits(hmi_path)

    if a_dat is None or h_dat is None:
        print(f"Failed to read input files for alignment.")
        return False

    # 3. Create WCS objects
    wcs_a = WCS(a_hdr)
    wcs_h = WCS(h_hdr)

    print(f"Reprojecting {os.path.basename(aia_path)} to HMI grid...")

    # 4. Perform the reprojection
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        reproj, _ = reproject_interp((a_dat, wcs_a), wcs_h, shape_out=h_dat.shape)

    # 5. Create the new header based on HMI
    new_hdr = h_hdr.copy()
    new_hdr['HISTORY'] = 'AIA reprojected to HMI grid using SpotiPy'

    # 6. Save the new FITS file
    fits.writeto(out_path, reproj, new_hdr, overwrite=True)

    return True
