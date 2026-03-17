"""
Downloading Module
------------------
A general-purpose tool for downloading any JSOC/SDO data series.

This module provides a single clean function `download_series()` that can be
used standalone for any JSOC series, or called in a loop from a pipeline
script for multiple series.

Example (standalone use):
    from spotipy.downloading import download_series

    download_series(
        series    = "hmi.M_720s",
        start_date= "2019-04-07T00:00:00",
        days      = 11,
        cadence_h = 6,
        email     = "you@institute.de",
        out_dir   = "data/M/"
    )
"""

import os
from sunpy.net import Fido, attrs as a
import astropy.units as u
from astropy.time import Time, TimeDelta


def ask_yn(msg: str) -> bool:
    """Prompts the user for a yes/no answer."""
    while True:
        ans = input(f"{msg} (y/n): ").strip().lower()
        if ans in {"y", "n"}:
            return ans == "y"
        print("Please answer 'y' or 'n'.")


def save_file_list(out_dir, list_path):
    """
    Saves a sorted list of all .fits files in out_dir to a text file.
    Allows the pipeline to reload file paths without re-querying JSOC.
    """
    files = sorted([f for f in os.listdir(out_dir) if f.endswith(".fits")])
    with open(list_path, "w") as fh:
        for f in files:
            fh.write(f + "\n")
    print(f"  File list saved: {list_path} ({len(files)} files)")


def load_file_list(list_path, out_dir):
    """
    Loads file paths from a saved list if it exists, otherwise falls back
    to scanning the directory directly for .fits files.
    This ensures manually placed FITS files are always found.
    """
    if os.path.exists(list_path):
        with open(list_path, "r") as fh:
            files = [line.strip() for line in fh if line.strip()]
        if files:
            return [os.path.join(out_dir, f) for f in files]

    # No list file found — scan the folder directly
    if os.path.isdir(out_dir):
        files = sorted([
            os.path.join(out_dir, f)
            for f in os.listdir(out_dir)
            if f.endswith(".fits")
        ])
        if files:
            print(f"  [INFO] No file list found for {out_dir} — scanned folder directly ({len(files)} files)")
        return files

    return []


def download_series(
    series,
    start_date,
    days,
    cadence_h,
    email,
    out_dir,
    wavelength_angstrom=None,
    preview=True,
):
    """
    Downloads a single JSOC data series.

    Can be used standalone for any series, or called in a loop from a
    pipeline script for multiple series.

    Parameters
    ----------
    series : str
        JSOC series name, e.g. 'hmi.Ic_720s', 'aia.lev1_uv_24s'.
    start_date : str
        ISO format start time, e.g. '2019-04-07T18:17:00'.
    days : int or float
        Number of days to cover.
    cadence_h : int or float
        Cadence in hours between downloaded frames.
    email : str
        JSOC-registered email address for data export notifications.
    out_dir : str
        Directory to save downloaded FITS files.
    wavelength_angstrom : float, optional
        If provided, adds a wavelength filter to the query.
        Required for AIA to avoid downloading all 10 wavelengths
        (e.g. wavelength_angstrom=1700).
    preview : bool
        If True, shows a preview of search results and asks for
        confirmation before downloading. Default is True.

    Returns
    -------
    list
        Sorted list of downloaded file paths. Empty list if failed or
        no records found.
    """
    os.makedirs(out_dir, exist_ok=True)

    safe_series = series.replace(".", "_")
    list_path = os.path.join(out_dir, f"{safe_series}_files.txt")

    existing = load_file_list(list_path, out_dir)
    if existing:
        print(f"\n[{series}] Found {len(existing)} existing files.")
        if not ask_yn(f"Re-download {series}?"):
            print(f"  Skipping — using existing files.")
            return existing

    # Build query
    t_start = Time(start_date)
    t_end   = t_start + TimeDelta(days * 24 * 3600, format="sec")

    query = [
        a.Time(t_start.iso, t_end.iso),
        a.jsoc.Series(series),
        a.Sample(cadence_h * u.hour),
        a.jsoc.Notify(email),
    ]
    if wavelength_angstrom is not None:
        query.append(a.Wavelength(wavelength_angstrom * u.angstrom))

    print(f"\n[{series}] Searching JSOC: {t_start.iso}  →  {t_end.iso}  (cadence: {cadence_h}h)")
    res = Fido.search(*query)

    if len(res) == 0:
        print(f"  No records found for {series}.")
        save_file_list(out_dir, list_path)
        return []

    print(f"  Records found: {len(res)}")

    if preview:
        print(res.show("TELESCOP", "INSTRUME", "T_OBS"))
        if not ask_yn(f"Download {len(res)} files for {series}?"):
            print("  Download skipped by user.")
            return []

    print(f"  Fetching {series}...")
    downloaded = Fido.fetch(res, path=os.path.join(out_dir, "{file}"))
    files = sorted(list(downloaded))
    print(f"  Downloaded {len(files)} files to: {out_dir}")
    save_file_list(out_dir, list_path)
    return files
