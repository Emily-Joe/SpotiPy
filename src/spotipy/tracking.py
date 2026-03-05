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


def get_time_steps(file_array, directory):
    """
    calculates the time beween two consecutive FITS-files (in the order of thier position in file_name_array). To be corrosponsive to file_name_array, the first entry is set to 0, because the duration to the previous file is zero, since its the first.

    Parameters
    ----------
    length: int
            length of the file-name-list

    Returns
    -------
    durations
            array of the amount of time between the corresponding two images in the name-list. This array corresponds to the file_name_array
    :Authors:
        Emily Joe Loessnitz (2024)
    """
    print('Calculating Time Steps...')

    time_steps = [0]  #set 0 as first entry, since the is no previous entry

    length = len(file_array)
    # for every file after the first:
    for i in range(1, length):
        hdr = f.getheader(os.path.join(directory, file_array[i]))
        hdr = f.open(os.path.join(directory, file_array[i]))
        time_i = hdr[1].header['DATE-OBS']
        # get the time_stamp from the previous image
        hdr2 = f.getheader(os.path.join(directory, file_array[i-1]))
        hdr2 = f.open(os.path.join(directory, file_array[i-1]))
        time_i0 = hdr2[1].header['DATE-OBS']
        #divide time from image with the one from the previous file
        time_step = (Time(time_i)-Time(time_i0)).value
        time_steps.append(time_step) # add the caculated time difference to the array

    time_steps = (time_steps)*u.hour*24 #make sure durations are in the correct units

    print('Done!')

    return time_steps ;


def calculate_tracks(file_array, directory, x0, y0, time_steps):
    """
    Takes a starting position (in arcsec) and calculates the next position of this point after a desired time step in accordance with differential rotation rates at the given latitude of the starting point.

    Parameters
    ----------
    start_time : str
        ISO format start time.
    start_pos_arcsec : tuple
        The (x, y) coordinates in arcseconds.
    duration_days : float
        Total duration to track. A full crossing of the Solar Disc will take up to 11 days.
    time_steps : list
        A list of Time objects.

    Returns
    -------
    track_coords : list of tuples
        The (x, y) coordinates in arcseconds for each time step.
    """

    print('Calculating frame locations...')
    ##############################################################################
    # First, load an observation and define a coordinate in its coordinate
    # frame (here, helioprojective Cartesian).  The appropriate rate of rotation
    # is determined from the heliographic latitude of the coordinate.

    #make sunpy Map from the information of the first FITS file of the series
    name  = file_array[0]
    data  = f.getdata(os.path.join(directory, name))
    hdr   = f.getheader(os.path.join(directory, name))
    hdr['cunit1'] = 'arcsec'
    hdr['cunit2'] = 'arcsec'
    aiamap = sunpy.map.Map(data, hdr)

    #defines start point which we can rotate
    start_point = SkyCoord(x0*u.arcsec, y0*u.arcsec, frame=aiamap.coordinate_frame)

    ##############################################################################
    # We can differentially rotate this coordinate by using
    # `~sunpy.coordinates.metaframes.RotatedSunFrame` with an array of observation
    # times

    diffrot_point = SkyCoord(RotatedSunFrame(base=start_point, duration=np.cumsum(time_steps)))

    ##############################################################################
    # To see what this coordinate looks like in "real" helioprojective
    # Cartesian coordinates, we can transform it back to the original frame.
    # Since these coordinates are represented in the original frame, they will not
    # account for the changing position of the observer over this same time range.

    transformed = diffrot_point.transform_to(aiamap.coordinate_frame)

    x_arc = transformed.Tx.to_value(u.arcsec)
    y_arc = transformed.Ty.to_value(u.arcsec)

    tracks = np.column_stack((x_arc, y_arc))

    print('Done!')

    return tracks


def track_region(image, header, tracks, frame_size, i):
    """
    Crops out a frame of specified size around the active region or feature
    """
    print(f'Tracking Region...step {i}')
    # get image height H and width W of full image
    H, W = image.shape[:2]

    #set up header of the FITS file form 'image'
    #hdr = f.getheader(image.path)
    #hdr = f.open(image.path)


    # get center pixels that align with center of the sun and arcsec-per-pixel for x- and y-direction:
    cx = header['CRPIX1']
    dx = header['CDELT1']
    cy = header['CRPIX2']
    dy = header['CDELT2']

    # the rotated frame locations (x_arc and y_arc) must be converted into pixels using the header information
    x_arc = tracks[i, 0]  #from previous step
    y_arc = tracks[i, 1]  #from previous step

    # converst to pixels:
    x_pix = int(round(cx - (x_arc / dx)))
    y_pix = int(round(cy - (y_arc / dy)))

    #y_pix = int(round(cy - (y_start_arcsec / dy))) if (cy and dy) else H//2
    #x_pix = int(round(x_hint_pix))


    # check if frame is outside of original FITS files range. This can happen with a big frame size and features close to either side of the solar edge
    # FAIL-SAVE:

    half_frame = int(frame_size//2) #rename for convinience

    # OUTDATED:

    if x_pix<(half_frame):
        #set minimum position to 0 + half_frame:
        x_f = int(half_frame)
        pic_crop = image[y_pix-half_frame:y_pix+half_frame, 0:int(frame_size)]

        #frame_location.append(x_failsave) #overwrite old location, add new one

    elif x_pix>(W-(half_frame)):
        x_f = int(W-half_frame)
        #set maximum position to full width W - half frame
        pic_crop = image[y_pix-half_frame:y_pix+half_frame, (int(W)-int(frame_size)):int(W)]

    else:
        x_f = int(x_pix)
        pic_crop = image[y_pix-half_frame:y_pix+half_frame, x_pix-half_frame:x_pix+half_frame]

        #frame_location.append(x_location)

    y_f = y_pix
    # as a fail_save
    #x_f = int(np.clip(x_pix - half_frame, 0+half_frame, W-half_frame))
    #y_f = int(np.clip(y_pix - half_frame, 0+half_frame, H-half_frame))

    pic_crop = image[y_f-half_frame:y_f+half_frame, x_f-half_frame:x_f+half_frame]

    pic_rot = np.rot90(pic_crop, 2)              #rotate image the right way
    crop = np.nan_to_num(pic_rot, nan=1) #fill nan's outside the solar disc with 1's
    #(this last step was neccessary, because otherwise the nans would appear dark on png images made from the data, which would mess with the center-finding process)

    #update frame locations in case fail-save was used:
    x_arc_f =  (-1)* x_f + cx * dx


    #frame_location = np.column_stack((x_arc_f, y_arc))
    frame_pos = np.array([x_arc_f, y_arc], dtype=float)

    #print('Done!')

    #return (x_arc_f, y_f), crop

    return frame_pos, crop


def strip(series, directory, tracks, strip_height_arcsec=None, frame_size=None, animate=False, overlay=True, save_path=None):
    """
    Unified strip generator.

    Creates a time-summed strip cut-out from a series of full-disk images.
    Optionally overlays tracking frames or creates animation frames.

    Parameters
    ----------
    series : list of str
        FITS filenames (ordered in time)
    directory : str
        Directory containing FITS files
    tracks : ndarray (N, 2)
        x_arcsec, y_arcsec positions
    strip_height_arcsec : float or None
        Total strip height in arcsec (default: 200)
    frame_size : int or None
        Frame size in pixels (square). If None → no frames.
    animate : bool
        If True, return moving-frame sequence (for animation or GIF-making)
    overlay : bool
        If True, the strip will consist of the averaged or overlayed images of all files in "series"
    save_path : str or None
        Path to save static strip plot, should contain the desired name
    """

    # double check input
    if animate and frame_size is None:
        raise ValueError("animate=True requires frame_size to be specified!")
    frames = None   # default


    # Load reference image + header
    first_file = series[0]
    first_img = f.getdata(os.path.join(directory, first_file))
    first_img = np.flipud(np.fliplr(first_img))  # FITS orientation fix, because FITS are rotated 180deg compared to observation

    H, W = first_img.shape# Height, Width

    hdr = f.open(os.path.join(directory, first_file))[1].header
    cx = hdr['CRPIX1']
    cy = hdr['CRPIX2']
    dx = hdr['CDELT1']
    dy = hdr['CDELT2']

    # setting up strip geometry
    # Define strip thickness: (default = 400 arcsec, unless stated otherwise)
    if strip_height_arcsec is None:
        strip_height_arcsec = 400.0

    strip_half_arcsec = strip_height_arcsec / 2
    strip_half_pix = strip_half_arcsec / abs(dy)

    y_arc_center = np.mean(tracks[:, 1])
    y_pix_center = cy + y_arc_center / dy   #conversion to pixels using header info

    y_min = int(np.clip(y_pix_center - strip_half_pix, 0, H))   # strip window lower limit
    y_max = int(np.clip(y_pix_center + strip_half_pix, 0, H))   # strip window upper limit
    strip_height = y_max - y_min

    # set up extend in (Solar Disc-) arcsec coordinates:
    x_left  = (0 - cx) * dx
    x_right = (W - cx) * dx
    y_min_arc = (y_min - cy) * dy
    y_max_arc = (y_max - cy) * dy

    # Frame geometry (if requested)
    if frame_size is not None:
        half_x = (frame_size / 2) * abs(dx)
        half_y = (frame_size / 2) * abs(dy)

    # ============================================================
    # CASE 1: TIME-SUMMED STRIP (overlay=True)
    # ============================================================
    if overlay:

        print('Overlaying images...')
        cube = np.zeros((len(series), strip_height, W))

        for i, fname in enumerate(series):
            img = f.getdata(os.path.join(directory, fname))
            img = np.flipud(np.fliplr(img))
            cube[i] = img[y_min:y_max, :]

        strip_img = np.sum(cube, axis=0)
        print('Done!')

        # STATIC STRIP (no animation)
        if not animate:
            fig, ax = plt.subplots(figsize=(12, 3))

            ax.imshow(
                strip_img,
                origin='lower',
                cmap='hot',
                extent=[x_left, x_right, y_min_arc, y_max_arc]
            )

            ax.set_xlabel('x [arcsec]')
            ax.set_ylabel('y [arcsec]')
            #ax.set_title('time-summed strip')

            # optional frames
            if frame_size is not None:
                print('Drawing frames...')
                colors = plt.cm.viridis(np.linspace(0, 1, len(tracks)))

                for i, (x_arc, y_arc) in enumerate(tracks):
                    rect = plt.Rectangle(
                        (x_arc - half_x, y_arc - half_y),
                        2 * half_x,
                        2 * half_y,
                        linewidth=2,
                        edgecolor=colors[i],
                        facecolor='none',
                        alpha=0.7
                    )
                    ax.add_patch(rect)
                print('Done!')

            plt.tight_layout()
            plt.show()

            if save_path is not None:
                yn = input("Save strip? (y/n): ").strip().lower()
                if yn == 'y':
                    plt.savefig(
                        os.path.join(save_path, 'time_summed_strip.png'),
                        dpi=300
                    )
                    print("Saved.")

            plt.close(fig)
            return strip_img, None


        # ANIMATION MODE (moving frame)
        frames = []
        print('Drawing animation frames...')

        for idx in range(len(tracks)):
            fig, ax = plt.subplots(figsize=(12, 3))
            canvas = FigureCanvas(fig)

            ax.imshow(
                strip_img,
                origin='lower',
                cmap='hot',
                extent=[x_left, x_right, y_min_arc, y_max_arc]
            )


            # --- single moving frame ---
            x_arc, y_arc = tracks[idx]
            rect = plt.Rectangle(
                (x_arc - half_x, y_arc - half_y),
                2 * half_x,
                2 * half_y,
                linewidth=2.5,
                edgecolor='dodgerblue',
                facecolor='none',
                alpha=0.95
            )
            ax.add_patch(rect)

            fig.tight_layout()
            canvas.draw()

            w, h = canvas.get_width_height()
            img_rgb = np.frombuffer(
                canvas.buffer_rgba(), dtype='uint8'
            ).reshape(h, w, 4)[:, :, :3]

            frames.append(img_rgb)
            plt.close(fig)

        print('Done!')
        return strip_img, frames

    # ============================================================
    # CASE 2: TIME-RESOLVED STRIPS (overlay=False)
    # ============================================================
    else:
        if not animate:
            raise ValueError("overlay=False requires animate=True")

        frames = []
        print('Creating time-resolved strip frames...')

        for i, fname in enumerate(series):
            img = f.getdata(os.path.join(directory, fname))
            img = np.flipud(np.fliplr(img))
            strip_img = img[y_min:y_max, :]

            fig, ax = plt.subplots(figsize=(12, 3))
            canvas = FigureCanvas(fig)

            #ax.imshow(
            #    strip_img,
            #    origin='lower',
            #    cmap='hot',
            #    extent=[x_left, x_right, y_min_arc, y_max_arc]
            #)

            ax.imshow(
                strip_img,
                origin='lower',
                cmap = 'gray',
                extent=[x_left, x_right, y_min_arc, y_max_arc]
            )

            if frame_size is not None:
                x_arc, y_arc = tracks[i]
                rect = plt.Rectangle(
                    (x_arc - half_x, y_arc - half_y),
                    2 * half_x,
                    2 * half_y,
                    linewidth=2.5,
                    edgecolor='dodgerblue',
                    facecolor='none',
                    alpha=0.95
                )
                ax.add_patch(rect)

            fig.tight_layout()
            canvas.draw()

            w, h = canvas.get_width_height()
            img_rgb = np.frombuffer(
                canvas.buffer_rgba(), dtype='uint8'
            ).reshape(h, w, 4)[:, :, :3]

            frames.append(img_rgb)

            if save_path is not None:
                frame_dir = 'strip_frames'
                folder_path= os.path.join(save_path, frame_dir)
                os.makedirs(folder_path, exist_ok=True)
                plt.imsave(
                    os.path.join(folder_path, f"frame_{i:04d}.png"),
                    img_rgb
                )

            plt.close(fig)
            print('Done!')

        return None, frames


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
