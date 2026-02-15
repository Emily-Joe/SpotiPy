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

# RENAMED: calculate_tracks -> track_spots (Matches __init__.py)
def track_spots(start_time, start_pos_arcsec, duration_days, time_steps):
    """
    Calculates the track of a solar feature across the disk.

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
