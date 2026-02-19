"""
Tracking Module
---------------
This module calculates the expected position of solar features over time
using differential rotation models.
"""

import numpy as np
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import SkyCoord
from sunpy.coordinates import RotatedSunFrame, Helioprojective
from scipy.ndimage import label, center_of_mass

def track_spots(start_time, start_pos_arcsec, duration_days, time_steps):
    """
    Calculates the theoretical track of a feature using differential rotation.

    Parameters
    ----------
    start_time : str
        ISO format start time.
    start_pos_arcsec : tuple
        The (x, y) coordinates in arcseconds.
    duration_days : float
        Total duration to track.
    time_steps : list
        A list of Time objects.

    Returns
    -------
    track_coords : list of tuples
        The (x, y) coordinates in arcseconds for each time step.
    """

    t0 = Time(start_time)
    c0 = SkyCoord(
        start_pos_arcsec[0]*u.arcsec,
        start_pos_arcsec[1]*u.arcsec,
        frame=Helioprojective,
        obstime=t0,
        observer="earth"
    )

    print(f"Tracking feature starting at {start_pos_arcsec} arcsec...")

    tracks = []

    for t in time_steps:
        current_time = Time(t)
        delta_t = (current_time - t0).to(u.day)

        rotated_coord = RotatedSunFrame(base=c0, duration=delta_t)

        current_pos = rotated_coord.transform_to(
            Helioprojective(obstime=current_time, observer="earth")
        )

        x_arc = current_pos.Tx.to(u.arcsec).value
        y_arc = current_pos.Ty.to(u.arcsec).value

        tracks.append((x_arc, y_arc))

    return tracks

from scipy.ndimage import label, center_of_mass

def refine_centering(image_crop, current_x0, current_y0, frame_size):
    """
    Refines the crop window by centering on the darkest feature (the sunspot).

    This prevents the spot from drifting out of the cropped windown by calculating
    the center of mass of the umbra and shifting the coordinates to match.
    """
    # Find the average brightness to establish a baseline
    median_val = np.nanmedian(image_crop)
    if median_val <= 0:
        return current_x0, current_y0

    # Identify pixels that are significantly darker (the sunspot)
    # 0.90 is a standard threshold to isolate the spot from the quiet sun
    mask = image_crop <= (0.90 * median_val)
    labeled_pixels, num_features = label(mask)

    if num_features > 0:
        # Find the largest dark feature in the crop
        feature_sizes = np.bincount(labeled_pixels.ravel())
        feature_sizes[0] = 0  # Ignore the background
        main_feature_id = int(np.argmax(feature_sizes))

        # Calculate the center of the sunspot
        y_center, x_center = center_of_mass(labeled_pixels == main_feature_id)

        if np.isfinite(x_center) and np.isfinite(y_center):
            # Shift the original x0, y0 so the spot is in the middle of the frame
            refined_x0 = int(np.clip(current_x0 + x_center - frame_size / 2, 0, image_crop.shape[1]))
            refined_y0 = int(np.clip(current_y0 + y_center - frame_size / 2, 0, image_crop.shape[0]))
            return refined_x0, refined_y0

    return current_x0, current_y0
