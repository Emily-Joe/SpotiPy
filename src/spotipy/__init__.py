"""
SpotiPy: Solar Spot Analysis Tool
----------------------------------
A library for downloading, processing, and analysing solar features
from SDO/HMI and SDO/AIA data.
"""

from .downloading import download_series, load_file_list, save_file_list, ask_yn
from .aligning import align_images
from .limbdarkening_removal import remove_limb_darkening, get_header_geometry
from .segmentation import get_masks, get_aia_masks
from .tracking import track_spots, refine_centering
from .gui_tools import get_manual_coordinates

__all__ = [
    # Downloading
    "download_series",
    "load_file_list",
    "save_file_list",
    "ask_yn",
    # Aligning
    "align_images",
    # Limb Darkening
    "remove_limb_darkening",
    "get_header_geometry",
    # Segmentation
    "get_masks",
    "get_aia_masks",
    # Tracking
    "track_spots",
    "refine_centering",
    # GUI
    "get_manual_coordinates",
]
