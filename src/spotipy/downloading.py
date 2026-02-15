"""
Downloading Module
------------------
This module handles searching and retrieving solar data from public archives.
It wraps SunPy's Fido interface to provide a streamlined experience for
HMI and AIA data retrieval.
"""

import os
from datetime import timedelta
import astropy.units as u
from sunpy.net import Fido, attrs as a
from astropy.time import Time

# Default Series mapping
SERIES_MAP = {
    "hmi": "hmi.Ic_720s",
    "hmi_m": "hmi.M_720s",
    "hmi_ld": "hmi.Ld_720s",
    "aia": "aia.lev1_uv_24s"
}

def download_data(start_time, duration_days, instrument='hmi', cadence='12m', email=None, output_dir='data/'):
    """
    Downloads solar data for a specific time range and instrument.

    This function abstracts the complexity of Fido searches, automatically
    handling JSOC notification requirements and series selection.

    Parameters
    ----------
    start_time : str or astropy.time.Time
        The start date/time (e.g., "2014-05-20T00:00:00").
    duration_days : float
        The duration of the observation in days.
    instrument : str, optional
        The instrument key. Options: 'hmi', 'hmi_m', 'hmi_ld', 'aia'.
        Default is 'hmi' (Intensity Continuum).
    cadence : str, optional
        The sampling cadence (e.g., '12m', '1h'). Default is '12m'.
    email : str, optional
        Your email address. REQUIRED for downloading HMI data from JSOC.
    output_dir : str, optional
        The directory where files will be saved. Default is 'data/'.

    Returns
    -------
    files : list
        A list of absolute paths to the downloaded FITS files.

    Raises
    ------
    ValueError
        If email is missing for JSOC queries or if no data is found.
    """

    # 1. Input Validation
    if instrument.startswith('hmi') and email is None:
        raise ValueError("JSOC (HMI) data requires a registered email address.")

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 2. Time Handling
    t0 = Time(start_time)
    t1 = t0 + timedelta(days=duration_days)

    print(f"[Query] Searching for {instrument} from {t0.iso} to {t1.iso}...")

    # 3. Build the Search Attributes
    series_name = SERIES_MAP.get(instrument, instrument)

    search_attrs = [
        a.Time(t0, t1),
        a.Sample(u.Quantity(cadence))
    ]

    # Add Instrument-specific attributes
    if 'hmi' in instrument:
        search_attrs.append(a.jsoc.Series(series_name))
        search_attrs.append(a.jsoc.Notify(email))
    elif 'aia' in instrument:
        search_attrs.append(a.Instrument('AIA'))
        search_attrs.append(a.Wavelength(1700 * u.angstrom))
        if email:
            search_attrs.append(a.jsoc.Notify(email)) # Sometimes AIA is on JSOC too

    # 4. Search
    result = Fido.search(*search_attrs)

    if len(result) == 0:
        raise RuntimeError(f"No data found for {instrument} in this range.")

    print(f"[Found] {len(result)} records. Downloading to {output_dir}...")

    # 5. Download
    # We use {file} to keep original filenames
    files = Fido.fetch(result, path=os.path.join(output_dir, "{file}"))

    return files
