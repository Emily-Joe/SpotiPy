"""
Limb Dark Removal Module
-----------------
This module contains tools for correcting solar data artifacts,
specifically removing limb darkening using radial profiles.
"""

import numpy as np
from scipy.interpolate import interp1d
from astropy.io import fits

def remove_limb_darkening(image, center, radius_pix, mask=None):
    """
    Removes limb darkening by calculating and subtracting the radial profile.

    Parameters
    ----------
    image : np.ndarray
        The 2D solar image.
    center : tuple
        The (x, y) coordinates of the solar center in pixels.
    radius_pix : float
        The radius of the sun in pixels.
    mask : np.ndarray, optional
        A boolean mask of regions to EXCLUDE from the profile calculation
        (e.g., active regions/spots). True = Exclude.

    Returns
    -------
    corrected_image : np.ndarray
        The flattened image.
    """

    # 1. Calculate the Radial Profile
    # Get distance of every pixel from the center
    yy, xx = np.indices(image.shape)
    r = np.sqrt((xx - center[0])**2 + (yy - center[1])**2)

    # Identify the solar disk
    disk = r < radius_pix

    # Exclude masked regions (spots) from the background calculation
    if mask is not None:
        disk &= (~mask)

    # Bin the data by radius to find the average intensity at each distance
    r_int = r.astype(int)

    # Use weights to compute mean (sum of intensities / count of pixels)
    # We use a masked array to ignore non-disk pixels
    arr = np.ma.masked_array(image, mask=~disk)

    # Bincount is fast!
    # Sum of intensity per radial bin
    tbin = np.bincount(r_int[disk], weights=arr[disk])
    # Count of pixels per radial bin
    nr = np.bincount(r_int[disk])

    # Avoid division by zero
    radial_prof = tbin / np.maximum(nr, 1)

    # Normalize profile (optional, but good for stability)
    if radial_prof.size > 20:
        med = np.median(radial_prof[:20]) # Center brightness
        if med > 0:
            radial_prof = radial_prof / med

    # 2. Create the 2D Correction Plane
    # Interpolate the 1D profile back into a 2D image
    interp_func = interp1d(
        np.arange(len(radial_prof)),
        radial_prof,
        bounds_error=False,
        fill_value='extrapolate'
    )

    background_plane = interp_func(r.ravel()).reshape(image.shape)

    # 3. Correct the Image
    # Divide the original image by the background plane
    # Handle edges/nans gracefully
    with np.errstate(divide='ignore', invalid='ignore'):
        corrected = image / np.where((background_plane > 0), background_plane, 1.0)

    return corrected

def get_header_geometry(header):
    """Helper to extract center and radius from a FITS header."""
    cx = header.get('CRPIX1')
    cy = header.get('CRPIX2')
    dx = header.get('CDELT1')
    rsun = header.get('RSUN_OBS')

    if rsun and dx:
        r_pix = rsun / dx
    else:
        r_pix = None

    return (cx, cy), r_pix
