"""
GUI Tools Module
----------------
Provides interactive graphical interfaces for (solar) feature selection 
and coordinate calibration in solar physics datasets.
"""

import matplotlib
matplotlib.use('TkAgg')  # Required for interactive window popups
import matplotlib.pyplot as plt
import astropy.io.fits as f
import numpy as np

def get_manual_coordinates(fits_path):
    """
    Interactive tool to select a solar feature and calculate its arcsec coordinates.

    This function opens a FITS image, allows the user to click a target (e.g., 
    a sunspot center), and translates the pixel position into Helioprojective 
    arcseconds using the FITS header metadata.

    Parameters
    ----------
    fits_path : str
        Path to the FITS file used for coordinate calibration.

    Returns
    -------
    tuple or None
        A tuple (x_arcsec, y_arcsec) representing the clicked position.
        Returns None if the window is closed without a selection.
    """
    coords_arcsec = []
    
    # Open FITS and handle multi-extension formats (common in HMI/AIA)
    with f.open(fits_path) as hdul:
        # Check for data in primary or secondary HDU 
        data = hdul[1].data if len(hdul) > 1 and hdul[1].data is not None else hdul[0].data
        header = hdul[1].header if len(hdul) > 1 and hdul[1].data is not None else hdul[0].header

    # Extract WCS (World Coordinate System) metadata 
    # CRPIX: Reference pixel (Solar Center)
    # CDELT: Plate scale (Arcseconds per pixel)
    mid_x, mid_y = header['CRPIX1'], header['CRPIX2']
    delta_x, delta_y = header['CDELT1'], header['CDELT2']

    # Setup the interactive plot
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Normalize display for better feature visibility (clipping outliers)
    ax.imshow(data, origin='lower', cmap='gray', 
              vmin=np.nanpercentile(data, 1), vmax=np.nanpercentile(data, 99))
    
    ax.set_title("Manual Selection: Click the center of the Active Region")
    ax.set_xlabel("Pixel X")
    ax.set_ylabel("Pixel Y")

    def onclick(event):
        """Internal callback to handle mouse click events."""
        if event.xdata is not None and event.ydata is not None:
            # Transformation logic based on your Clicker tool [cite: 2]
            # (Pixel_Coordinate - Reference_Pixel) * Plate_Scale * Orientation_Multiplier
            x_arc = (event.xdata - mid_x) * delta_x * (-1)
            y_arc = (event.ydata - mid_y) * delta_y * (-1)
            
            coords_arcsec.append((x_arc, y_arc))
            print(f"[Input] Target Locked: x={x_arc:.2f}\", y={y_arc:.2f}\"")
            plt.close()

    # Connect the matplotlib event manager and wait for user input
    fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()
    
    return coords_arcsec[0] if coords_arcsec else None
