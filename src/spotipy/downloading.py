"""
Downloading Module
------------------
This module handles searching and retrieving solar data from public archives.
It wraps SunPy's Fido interface and automatically maps user requests
to the correct JSOC series (HMI) or AIA wavelength channels.
"""

import os
from datetime import timedelta
import astropy.units as u
from sunpy.net import Fido, attrs as a
from astropy.time import Time

def get_aia_series(wavelength_angstrom):
    """
    Returns the correct AIA series name based on the wavelength.

    Mappings:
    - 4500 Å       -> Visible (1 hour cadence)
    - 1600, 1700 Å -> UV (24 sec cadence)
    - Others       -> EUV (12 sec cadence)
    """
    w = int(wavelength_angstrom)

    if w == 4500:
        return "aia.lev1_vis_1h"
    elif w in [1600, 1700]:
        return "aia.lev1_uv_24s"
    else:
        # Default to EUV for 94, 131, 193, 211, 304, 335
        return "aia.lev1_euv_12s"

def download_data(start_time, duration_days, instrument='hmi', wavelength=None, cadence='12m', email=None, output_dir='data/'):
    """
    Downloads solar data for a specific time range, instrument, and optional wavelength.

    Parameters
    ----------
    start_time : str or astropy.time.Time
        The start date/time (e.g., "2017-09-06T12:00:00").
    duration_days : float
        The duration of the observation in days.
    instrument : str
        The instrument or physical observable.
        - For HMI: 'Ic', 'M', 'V', 'Ld', 'Lw', 'Ic_noLD' (or full series like 'hmi.M_45s')
        - For AIA: 'aia' (must provide wavelength)
    wavelength : int, optional
        Required for AIA. The wavelength in Angstroms (e.g., 1700, 4500, 1600).
    cadence : str, optional
        The sampling cadence for Fido (e.g., '12m', '1h', '45s').
        For HMI, this also helps select the series if not provided explicitly.
    email : str
        REQUIRED. Your registered JSOC email address.
    output_dir : str
        The directory where files will be saved.

    Returns
    -------
    files : list
        A list of absolute paths to the downloaded FITS files.
    """

    # 1. Input Validation
    if email is None:
        raise ValueError("JSOC data download requires a registered email address.")

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 2. Time Handling
    t0 = Time(start_time)
    t1 = t0 + timedelta(days=duration_days)

    print(f"[Query] Searching {instrument} ({wavelength if wavelength else ''}) from {t0.iso} to {t1.iso}...")

    # 3. Construct Search Attributes
    search_attrs = [
        a.Time(t0, t1),
        a.jsoc.Notify(email) # Always required for HMI/AIA via JSOC
    ]

    # --- HMI LOGIC ---
    if 'hmi' in instrument.lower() or instrument in ['Ic', 'M', 'V', 'Ld', 'Lw', 'Ic_noLD']:

        # If user passed a full series name (e.g., "hmi.M_45s"), use it directly
        if "hmi." in instrument:
            series_name = instrument

        # Otherwise, construct it smartly
        else:
            # Clean up the input (e.g., "M" -> "M")
            obs = instrument.replace("hmi_", "")

            # Special Case: Limb Darkening Removed (Only available in 720s)
            if obs == "Ic_noLD":
                series_name = "hmi.Ic_nolimbDark_720s"

            # Standard Observables (Ic, M, V, Ld, Lw)
            else:
                # Decide cadence suffix based on user input string
                if '45s' in cadence:
                    suffix = "45s"
                else:
                    suffix = "720s"

                series_name = f"hmi.{obs}_{suffix}"

        print(f"   -> Identified HMI Series: {series_name}")
        search_attrs.append(a.jsoc.Series(series_name))

        # Add sampling cadence to prevent getting too many files if Fido finds duplicates
        # (Note: For HMI 45s series, asking for 12m sample is valid to reduce data volume)
        if cadence:
             search_attrs.append(a.Sample(u.Quantity(cadence)))

    # --- AIA LOGIC ---
    elif 'aia' in instrument.lower():
        if wavelength is None:
            raise ValueError("For AIA, you MUST specify a 'wavelength' (e.g., 1600, 4500).")

        # 1. Pick the correct series based on wavelength (e.g., 4500 -> vis_1h)
        series_name = get_aia_series(wavelength)
        print(f"   -> Identified AIA Series: {series_name} for {wavelength}A")

        # 2. Add Attributes
        search_attrs.append(a.jsoc.Series(series_name))
        search_attrs.append(a.jsoc.Segment('image')) # Get the image segment
        search_attrs.append(a.Wavelength(int(wavelength) * u.angstrom))

        if cadence:
            search_attrs.append(a.Sample(u.Quantity(cadence)))

    else:
        raise ValueError(f"Unknown instrument: {instrument}")

    # 4. Search
    result = Fido.search(*search_attrs)

    if len(result) == 0:
        print(f"Error: No data found for {series_name}. Check your dates or cadence.")
        return []

    print(f"[Found] {len(result)} records. Downloading to {output_dir}...")

    # 5. Download
    # We use {file} to keep original filenames from JSOC
    files = Fido.fetch(result, path=os.path.join(output_dir, "{file}"))

    return files
