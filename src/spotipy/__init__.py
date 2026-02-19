"""
SpotiPy: Solar Spot Analysis Tool
---------------------------------
A library for downloading, processing, and tracking solar features.
"""

# Import the modules so users can access them directly
# e.g., spotipy.download_data()

from .downloading import download_data
from .aligning import align_images
from .limbdarkening_removal import remove_limb_darkening
from .segmentation import get_masks, get_aia_masks
from .tracking import track_spots, refine_centering
from .gui_tools import get_manual_coordinates

__all__ = [
    "download_data",
    "aligning",
    "remove_limb_darkening",
    "get_masks",
    "get_aia_masks",
    "track_spots",
    "refine_centering",
    "get_manual_coordinates"
]
