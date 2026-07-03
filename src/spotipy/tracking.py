"""
Tracking Module
---------------
This module calculates the expected position of solar features over time
using differential rotation models.
"""
import os
import numpy as np
import numpy.ma as ma
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cv2
import astropy.io.fits as f
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import SkyCoord
import sunpy.map
from sunpy.coordinates import RotatedSunFrame
from scipy.ndimage import label, center_of_mass



def get_time_steps(file_array, directory):
    """
    Calculate time differences between consecutive FITS files.
    The function reads the observation time ('DATE-OBS') from each FITS file and computes the time difference to the previous file. The first entry    is set to 0 since there is no previous file.

    Parameters
    ----------
    file_array : list of str
        List of FITS file names in chronological order.
    directory : str
        Path to the directory containing the FITS files.

    Returns
    -------
    time_steps : numpy.ndarray
        Array of time differences in days between consecutive files. The first element is 0. The array length matches `file_array`.
    """
    print('Calculating Time Steps...')

    time_steps = [0]
    # for every file after the first:
    for i in range(1, len(file_array)):
        hdr1 = f.open(os.path.join(directory, file_array[i]))
        hdr2 = f.open(os.path.join(directory, file_array[i-1]))

        t1 = Time(hdr1[1].header['DATE-OBS'])
        t0 = Time(hdr2[1].header['DATE-OBS'])

        dt = (t1 - t0).to(u.day).value # time steps in unit [days]
        time_steps.append(dt)

    time_steps = np.array(time_steps)

    print('Done!')

    return time_steps ;


def calculate_tracks(file_array, directory, x0, y0, time_steps):
    '''
    Compute the tracked position of a solar feature using differential rotation.

    A starting coordinate (in helioprojective arcseconds) is propagated over time using SunPy's differential rotation model. The time evolution is defined by the cumulative sum of the provided time steps.

    Parameters
    ----------
    file_array : list of str
        List of FITS file names. The first file is used to define the reference coordinate frame.
    directory : str
        Path to the directory containing the FITS files.
    x0 : float
        Initial x-coordinate in arcseconds (helioprojective).
    y0 : float
        Initial y-coordinate in arcseconds (helioprojective).
    time_steps : array-like of float
        Time differences between consecutive frames in days. The first entry should be 0.

    Returns
    -------
    tracks : numpy.ndarray
        Array of shape (N, 2) containing the tracked (x, y) positions in arcseconds for each time step.
    '''

    print('Calculating frame locations...')
    ##############################################################################
    # First, load an observation and define a coordinate in its coordinate frame (here, helioprojective Cartesian).  The appropriate rate of rotation is determined from the heliographic latitude of the coordinate.

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
    # We can differentially rotate this coordinate by using `~sunpy.coordinates.metaframes.RotatedSunFrame` with an array of observation times

    diffrot_point = SkyCoord(RotatedSunFrame(base=start_point, duration=np.cumsum(time_steps)))

    ##############################################################################
    # To see what this coordinate looks like in "real" helioprojective Cartesian coordinates, we can transform it back to the original frame. Since these coordinates are represented in the original frame, they will not account for the changing position of the observer over this same time range.

    transformed = diffrot_point.transform_to(aiamap.coordinate_frame)

    x_arc = transformed.Tx.to_value(u.arcsec)
    y_arc = transformed.Ty.to_value(u.arcsec)

    tracks = np.column_stack((x_arc, y_arc))

    print('Done!')

    return tracks


def track_region(image, hdr, tracks, frame_size, i=None):
    """
    Extract a cropped region around a tracked solar feature.

    The function converts tracked helioprojective coordinates (arcsec)into pixel coordinates using FITS header information, crops a square region around that position, and applies basic post-processing.

    Parameters
    ----------
    image : numpy.ndarray
        Two-dimensional image array containing the FITS data.
    hdr : astropy.io.fits.HDUList
        FITS file object containing the image header information
        (extension ``[1]`` is used).
    tracks : numpy.ndarray
        Either an array of shape ``(N, 2)`` containing tracked
        helioprojective coordinates (x, y) in arcseconds for multiple
        frames, or a single coordinate pair ``(x, y)`` if ``i`` is ``None``.
    frame_size : int
        Size of the square crop in pixels.
    i : int, optional
        Index of the current frame in ``tracks``. If ``None``, ``tracks`` is
        interpreted as a single coordinate pair.

    Returns
    -------
   frame_pos : numpy.ndarray
        Corrected helioprojective coordinates ``(x, y)`` in arcseconds.
        The x-coordinate may differ from the input if it was adjusted to keep
        the crop within the image boundaries.
    crop : numpy.ndarray
        Square image crop of shape ``(frame_size, frame_size)`` with NaN
        values replaced by ``1``.
    """

    # get image height H and width W of full image
    H, W = image.shape[:2]

    # get center pixels that align with center of the sun and arcsec-per-pixel for x- and y-direction:
    cx = hdr[1].header['CRPIX1']
    dx = hdr[1].header['CDELT1']
    cy = hdr[1].header['CRPIX2']
    dy = hdr[1].header['CDELT2']

    if i is not None:
        # the rotated frame locations (x_arc and y_arc) must be converted into pixels using the header information
        x_arc = tracks[i, 0]  #from previous step
        y_arc = tracks[i, 1]  #from previous step
    else:
        x_arc = tracks[0]
        y_arc = tracks[1]

    # converst to pixels:
    x_pix = int(round(cx + (x_arc / dx)))
    y_pix = int(round(cy + (y_arc / dy)))

    # FAIL-SAVE for when the frame would extend over the initial full disk image size:
    half_frame = int(frame_size//2) #rename for convinience

    if x_pix<(half_frame):
        #set minimum position to 0 + half_frame:
        x_f = int(half_frame)
        pic_crop = image[y_pix-half_frame:y_pix+half_frame, 0:int(frame_size)]

    elif x_pix>(W-half_frame):
        #frame_location.append(x_failsave) #overwrite old location, add new one
        x_f = int(W-half_frame)
        #set maximum position to full width W - half frame
        pic_crop = image[y_pix-half_frame:y_pix+half_frame, (int(W)-int(frame_size)):int(W)]

    else:
        x_f = int(x_pix)
        pic_crop = image[y_pix-half_frame:y_pix+half_frame, x_pix-half_frame:x_pix+half_frame]


    y_f = y_pix # no fail save in y-direction needed, since most active regions are not near the poles
    pic_crop = image[y_f-half_frame:y_f+half_frame, x_f-half_frame:x_f+half_frame]
    crop = np.nan_to_num(pic_crop, nan=1)        # fill nan's outside the solar disc with 1's

    #update frame locations in case fail-save was used:
    x_arc_f = (x_f - cx) * dx

    frame_pos = np.array([x_arc_f, y_arc], dtype=float)

    return frame_pos, crop


def strip(series, directory, tracks, strip_height_arcsec=None, frame_size=None, animate=False, overlay=True, save_path=None):
    """
    Generate time-summed or time-resolved strip images from a series of FITS files.

    The function extracts a horizontal strip (in helioprojective coordinates) centered on the tracked feature and either:
    - sums all strips into a single image (overlay=True), or
    - creates a time-resolved sequence of strip images (overlay=False).

    Optionally, tracking frames can be overlaid and animation frames generated.

    Parameters
    ----------
    series : list of str
        FITS filenames in chronological order.
    directory : str
        Path to the directory containing the FITS files.
    tracks : numpy.ndarray
        Array of shape (N, 2) containing (x, y) positions in arcseconds.
    strip_height_arcsec : float, optional
        Total height of the strip in arcseconds. Default is 400.
    frame_size : int, optional
        Size of the square tracking frame in pixels. If None, no frames are drawn.
    animate : bool, optional
        If True, returns a list of RGB frames for animation.
    overlay : bool, optional
        If True, produces a time-summed strip. If False, produces time-resolved strips (requires animate=True).
    save_path : str, optional
        Directory path for saving output images or frames.

    Returns
    -------
    strip_img : numpy.ndarray or None
        Time-summed strip image if overlay=True, otherwise None.
    frames : list of numpy.ndarray or None
        List of RGB frames for animation if animate=True, otherwise None.
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

def refine_centering(image_crop, current_x0, current_y0, frame_size, max_shift=None, tracking_mode="largest"):
    """
    Refine the center position of a cropped image by locating the darkest
    feature within the frame.

    The image is first smoothed with a Gaussian filter before dark regions
    are identified using an intensity threshold. Depending on the selected
    tracking mode, the function either centers on the largest contiguous
    dark feature or on the center of mass of all detected dark pixels. If
    the calculated shift exceeds ``max_shift`` (if specified), the refinement
    is rejected and the original center is returned.

    Parameters
    ----------
    image_crop : numpy.ndarray
        Two-dimensional cropped image array.
    current_x0 : int or float
        Current x-coordinate (pixel) of the crop center in the original image.
    current_y0 : int or float
        Current y-coordinate (pixel) of the crop center in the original image.
    frame_size : int
        Size of the square crop in pixels.
    max_shift : float, optional
        Maximum allowed correction distance (in pixels) between the current
        crop center and the newly determined center. If the required shift
        exceeds this value, the refinement is rejected. If ``None``, no
        maximum shift is enforced.
    tracking_mode : {"largest", "evolving"}, optional
        Method used to determine the feature center.

        - ``"largest"``: Track the largest contiguous dark feature.
        - ``"evolving"``: Track the center of mass of all detected dark pixels,
          which is useful for extended or evolving active regions.

    Returns
    -------
    refined_x0 : int
        Refined x-coordinate (pixel) of the crop center. If the refinement is
        rejected, the original coordinate is returned.
    refined_y0 : int
        Refined y-coordinate (pixel) of the crop center. If the refinement is
        rejected, the original coordinate is returned.
    accepted : bool
        Whether the refinement was accepted. Returns ``False`` if no valid
        feature was detected or if the required correction exceeded
        ``max_shift``.
    """


    # Smooth the image via a gaussian blur to not be reliant on individual dark pixels:
    image_blur = cv2.GaussianBlur(image_crop,(7,7),0)

    # Find the average brightness to establish a baseline
    median_val = np.nanmedian(image_blur)
    if median_val <= 0:
        return current_x0, current_y0, False


    # Identify pixels that are significantly darker (the sunspot)
    mask = (image_blur <= (0.90 * median_val)) & (image_blur > 0) & np.isfinite(image_blur)

    labeled_pixels, num_features = label(mask)

    if num_features > 0:
        if tracking_mode == "largest":
            # Find the largest dark feature in the crop
            feature_sizes = np.bincount(labeled_pixels.ravel())
            feature_sizes[0] = 0  # Ignore the background
            main_feature_id = int(np.argmax(feature_sizes))

            # Calculate the center of the largest spot
            y_center, x_center = center_of_mass(labeled_pixels == main_feature_id)

        elif tracking_mode == "evolving":
            # Calculate the center of mass using ALL detected dark pixels
            y_center, x_center = center_of_mass(mask)

        else:
            print(f"[WARN] Unknown tracking_mode '{tracking_mode}', defaulting to largest.")
            feature_sizes = np.bincount(labeled_pixels.ravel())
            feature_sizes[0] = 0
            y_center, x_center = center_of_mass(labeled_pixels == int(np.argmax(feature_sizes)))

        if np.isfinite(x_center) and np.isfinite(y_center):
            # Shift the original x0, y0 so the spot is in the middle of the frame
            refined_x0 = int(current_x0  - (frame_size / 2) + x_center)
            refined_y0 = int(current_y0  - (frame_size / 2) + y_center)

            # Check distance from center
            accepted = True
            if max_shift is not None:
                correction_dist = np.sqrt((x_center - frame_size / 2)**2 + (y_center - frame_size / 2)**2)
                if correction_dist > max_shift:
                    accepted = False
                    return current_x0, current_y0, accepted

            return refined_x0, refined_y0, accepted

    # If no valid features are found, reject the refinement
    return current_x0, current_y0, False


def center(image, mask=None, threshold=None, ellipse=False, save_path=None):
    """
    Determine the center of a solar feature in an image.

    The function identifies the largest detected feature and computes its
    center either from its image moments (centroid) or by fitting an ellipse.
    If no binary mask is provided, one is generated from the input image using
    thresholding and morphological opening.

    Parameters
    ----------
    image : numpy.ndarray
        Two-dimensional grayscale image.
    mask : numpy.ndarray, optional
        Binary mask of the feature of interest. If ``None``, a mask is
        generated from ``image`` using intensity thresholding and
        morphological operations.
    threshold : int, optional
        Intensity threshold (0–255) used to generate the binary mask when
        ``mask`` is ``None``. If ``None``, a default threshold of ``175`` is
        used.
    ellipse : bool, optional
        If ``True``, determine the center from an ellipse fitted to the
        largest contour. Otherwise, compute the centroid from the image
        moments of the feature mask.
    save_path : str, optional
        File path prefix for saving a diagnostic image showing the detected
        contour and center. The suffix ``"_elli.png"`` or ``"_cent.png"``
        is appended automatically depending on the selected method.

    Returns
    -------
    cx : int
        x-coordinate of the detected feature center in pixel coordinates.
    cy : int
        y-coordinate of the detected feature center in pixel coordinates.
    """

    # Fail save if mask was not created with segmentation module:
    if mask is None:

        if threshold is None:
            print ('unspecified threshold value! assumed value: 175 ')
            thresh = 175
        else :
            thresh = int(threshold)

        # if masks werent created before it will now create a spot mask:
        # normalize + convert
        grayscaleImage = cv2.normalize(image, None, 0, 255, cv2.NORM_MINMAX)
        grayscaleImage = grayscaleImage.astype(np.uint8)

        #blur the image to make granulation appear lighter in contrast to sunspots
        blur = cv2.GaussianBlur(grayscaleImage,(7,7),0)
        threshValue, binaryImage = cv2.threshold(blur,thresh,255,cv2.THRESH_BINARY_INV)
        kernel = np.ones((3, 3), np.uint8)  #size of the kernel
        opIterations = 3                    #number of iterations
        # Perform opening:
        openingImage = cv2.morphologyEx(binaryImage, cv2.MORPH_OPEN, kernel, None, None, opIterations, cv2.BORDER_REFLECT101)

        mask = np.zeros_like(openingImage) # make a mask from the binary image
        # find the countours of the central figure
        contours,_ = cv2.findContours(openingImage, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        contours    = [c for c in contours if cv2.contourArea(c) > 100]
        # only consider contours that have an area [pixels] over 100 to exclude noise or fine features in granulation
        # we assume the sunspot is the biggest object in frame, so we only consider the largest found contour
        if not contours:
            raise ValueError("No valid feature contour found.")
        biggest_contour = max(contours, key = cv2.contourArea)
        cv2.drawContours(mask, [biggest_contour], -1, 255, thickness=cv2.FILLED) #draw in contours to check

    else :
        # analog to above but with pre-selected mask of feature (sunspot, plage, network,...)
        contours,_ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        contours    = [c for c in contours if cv2.contourArea(c) > 100]
        biggest_contour = max(contours, key = cv2.contourArea)
        cv2.drawContours(mask, [biggest_contour], -1, 255, thickness=cv2.FILLED)

    # finding the center:
    if ellipse :
        ellipse = cv2.fitEllipse(biggest_contour)
        (x, y), (w, h), ang = cv2.fitEllipse(biggest_contour) # fit ellipse and save the parameters
        cx = int(x)    #center x-coordinate
        cy = int(y)    #center x-coordinate

        if save_path is not None:
            #draw  and save image with ellipse + center
            output = cv2.cvtColor(image, cv2.COLOR_GRAY2BGR)
            output = cv2.drawContours(output, contours, -1, (255,0,0), 1)
            output = cv2.ellipse(output, ellipse, (0, 0, 255), 2)
            output = cv2.line(output, (cx,cy), (cx,cy), (0,255,0), 10)
            plt.imshow(output, cmap='gray')
            plt.gca().invert_yaxis()    # invert axis for plotting
            plt.savefig(save_path+'_elli.png')
        plt.close()

    else:
        # Calculate the moments & compute centroid
        imageMoments = cv2.moments(mask)
        cx = int(imageMoments['m10']/imageMoments['m00'])
        cy = int(imageMoments['m01']/imageMoments['m00'])

        if save_path is not None:
            # Draw centroid onto the overlay image:
            bgrImage = cv2.cvtColor(image, cv2.COLOR_GRAY2BGR)
            bgrImage = cv2.drawContours(bgrImage, contours, -1, (255,0,0), 1)
            bgrImage = cv2.line(bgrImage, (cx,cy), (cx,cy), (0,255,0), 10)
            plt.imshow(bgrImage, cmap = 'gray')
            plt.gca().invert_yaxis()     # invert axis for plotting
            plt.savefig(save_path+'_cent.png')
        plt.close()

    return cx, cy


def position_conversion(center, frame, frame_size, solar_center, pix_to_arc, flip_x=1, flip_y=1):
    '''
    Convert local frame coordinates into global pixel coordinates and helioprojective arcseconds relative to the solar disk center.

    Parameters
    ----------
    center : array-like, shape (N, 2)
        Pixel coordinates of the feature(s) within each frame (e.g. [x, y]).
    frame : array-like, shape (N, 2)
        Pixel coordinates of the frame origin (e.g. top-left corner) in the global image for each timestep. Must have the same shape as center.
    frame_size : float
        Size of the frame in pixels (assumed square).
    solar_center : array-like, shape (2,)
        Pixel coordinates of the solar disk center in the full image (constant reference).
    pix_to_arc : float
        Conversion factor from pixels to arcseconds.
    flip_x : int, optional
        Factor to control x-axis orientation (+1 or -1). Default is +1.
    flip_y : int, optional
        Factor to control y-axis orientation (+1 or -1). Default is +1.

    Returns
    -------
    position_pix : ndarray
        Global pixel coordinates relative to the full image.
    position_arc : ndarray
        Coordinates in arcseconds relative to the solar disk center.

    Notes
    -----
    The sign of the coordinate axes may depend on the FITS header conventions (e.g. CDELT values). Axis flips should be handled via flip_x and flip_y instead of hard-coded sign changes.
    '''

    center = np.asarray(center)
    frame = np.asarray(frame)

    frame_pix = solar_center + (frame / pix_to_arc)     # frame pixel position
    position_pix = frame_pix - frame_size/2 + center    # global center pixel position
    position_arc = (position_pix - solar_center) * pix_to_arc   # convert to solar heliographic arcsec

    # fail save if orientation of input data is flipped (then set flip to -1, default=1)
    position_arc[:, 0] *= flip_x
    position_arc[:, 1] *= flip_y

    return position_pix, position_arc


def extract_header_params(file_array, directory):
    """
    Extract key observational and geometric parameters from FITS headers.

    Parameters
    ----------
    file_array : list of str
        FITS filenames.
    directory : str
        Path to the directory containing the FITS files.

    Returns
    -------
    x_center_pix : numpy.ndarray
        Reference pixel x-coordinates (CRPIX1).
    y_center_pix : numpy.ndarray
        Reference pixel y-coordinates (CRPIX2).
    dx : numpy.ndarray
        Pixel scale in x-direction (arcsec/pixel).
    dy : numpy.ndarray
        Pixel scale in y-direction (arcsec/pixel).
    R_sun : numpy.ndarray
        Observed solar radius (RSUN_OBS) in arcseconds.
    L : numpy.ndarray
        Carrington longitude (CRLN_OBS) in degrees.
    B : numpy.ndarray
        Carrington latitude (CRLT_OBS) in degrees.
    time_stamps : numpy.ndarray
        Observation timestamps (DATE-OBS) as strings.
    """

    x_center_pix = []   # center pixel x-coordinate
    y_center_pix = []   # center pixel y-coordinate
    dx = []             # arcsec per pixel conversion in x-direction
    dy = []             # arcsec per pixel conversion in y-direction
    R_sun_list = []     # observed solar radius [in arcsec]
    L_list = []         # Carrington Longitude [deg]
    B_list = []         # Carrington Latitude [deg]
    time_stamps = []    # exact timestamp of observation in UT

    for name in file_array:
        path = os.path.join(directory, name)
        header = f.getheader(path, ext=1)

        x_center_pix.append(header['CRPIX1'])
        y_center_pix.append(header['CRPIX2'])
        dx.append(header['CDELT1'])
        dy.append(header['CDELT2'])
        R_sun_list.append(header['RSUN_OBS'])
        L_list.append(header['CRLN_OBS'])
        B_list.append(header['CRLT_OBS'])
        time_stamps.append(header['DATE-OBS'])


    return np.array(x_center_pix), np.array(y_center_pix), np.array(dx), np.array(dy), np.array(R_sun_list),  np.array(L_list), np.array(B_list), np.array(time_stamps)


def differential_rotation_rate(positions, durations, R_sun, L, B):
    """
    Compute the solar rotation rate from tracked feature positions.The function converts helioprojective coordinates (arcsec) into Carrington longitudes and derives the rotation rate from their temporal evolution.

    Parameters
    ----------
    positions : numpy.ndarray, shape (N, 2)
        Feature positions (x, y) in arcseconds.
    durations : numpy.ndarray, shape (N,)
        Time differences between frames in days (first entry typically 0).
    R_sun : numpy.ndarray, shape (N,)
        Observed solar radius in arcseconds.
    L : numpy.ndarray, shape (N,)
        Carrington longitude of disk center in degrees.
    B : numpy.ndarray, shape (N,)
        Carrington latitude (solar tilt) in degrees.

    Returns
    -------
    omega : numpy.ma.MaskedArray
        Estimated rotation rate in deg/day (including Carrington base rate).
    time_axis : numpy.ma.MaskedArray
        Time axis in hours corresponding to rotation rate values.
    """

    positions = np.asarray(positions)
    x_all = positions[:, 0]
    y_all = positions[:, 1]

    L = np.deg2rad(L)   #convert to radiants
    B = np.deg2rad(B)   #convert to radiants

    l_deg_array = []    # set up array for carrington latitude

    for i, (x, y) in enumerate(zip(x_all, y_all)):
        # we need position of sunspot-center in heliocentric coordinates to converted to carrington system:
        SC = np.sqrt(x**2 + y**2)       # distance between spot (S) and sun-center (C)[arcsec]
        SC_rad = np.deg2rad(SC / 3600)  # converting to [rad]
        SC_vector = np.array([x, y])    # vector to the spot

        # for orientation: asign a sign to both sides of the solar disc center
        N = np.array([0, R_sun[i]])         # vector to north
        cross = np.cross(N, SC_vector)      # cross product
        alpha = np.arccos(y / SC) if SC != 0 else 0.0      # angle between SC and north
        alpha *= -1 if cross < 0 else 1  # if crossproduct negative then the angle is negative as well per definition (positive from north eastwards)

        # conversion to heliocentric coordinate system:
        sigma = np.arcsin(SC / R_sun[i]) - SC_rad       # heliocentric angle sigma [rad]
        # heliographic latitude:
        b = np.arcsin(np.cos(sigma) * np.sin(B[i]) + np.sin(sigma) * np.cos(B[i]) * np.cos(alpha))
        # heliographic longitude:
        longitude_central = np.arcsin(np.sin(alpha) * np.sin(sigma) / np.cos(b))
        longitude = L[i] - longitude_central
        l_deg_array.append(np.rad2deg(longitude))

    car_longitude = np.asarray(l_deg_array)

    # carrington rotation rate:
    longitude_diff = np.diff(car_longitude)
    # mask (=ignore) all points where the difference is smaller than 1. this case can only happen if the start of a new carringtion rotation started during the observation (in that case the difference would be negative!)
    mask = longitude_diff < 1
    masked_longitude_diff = ma.masked_where(~mask, longitude_diff)

    # timesteps/durations masking:
    durations = durations[1:] #this is new, but I think this is wrong.
    masked_durations = ma.masked_where(~mask, durations)

    # time axis:
    cumulative_time = np.cumsum(durations)                  # cumulative sum to set up time axis
    masked_time = ma.masked_where(~mask, cumulative_time)   # mask jumps as before
    time_axis_days = masked_time + durations / 2
    time_axis = time_axis_days * 24

    # Final Rotation Rate (deg/day) in Carrington System:
    omega = masked_longitude_diff/ masked_durations
    omega_carrington = omega + 14.184

    return omega_carrington, time_axis


def track_spots(file_array, directory, x0, y0):
    """
    wrapper to handle time steps and coordinate tracking only

    Parameters
    ----------
    file_array : list of str
        FITS filenames in chronological order.
    directory : str
        Path to FITS files.
    x0, y0 : float
        Initial feature position in arcseconds.

    Returns
    -------
    tracks : numpy.ndarray
        Array of shape (N, 2) containing the tracked (x, y) positions in arcseconds for each time step.
    """
    time_steps = get_time_steps(file_array, directory)
    tracks = calculate_tracks(file_array, directory, x0, y0, time_steps)
    return tracks


def sunspot_rotation_rate(file_array, directory, x0, y0, frame_size, threshold=None, refine=True, ellipse=False, save_path=None):
    """
    Wrapper to compute the solar rotation rate from a sequence of FITS images, specifically for sunspots.

    This function performs the full analysis pipeline:
    1. Compute time steps between observations
    2. Predict feature positions via differential rotation
    3. Extract header parameters
    4. Track and crop regions around the feature
    5. Determine feature center positions
    6. Convert coordinates to heliographic system
    7. Compute Carrington rotation rate

    Parameters
    ----------
    file_array : list of str
        FITS filenames in chronological order.
    directory : str
        Path to FITS files.
    x0, y0 : float
        Initial feature position in arcseconds.
    frame_size : int
        Size of the cropped region in pixels.
    threshold : int, optional
        Threshold for feature detection (passed to `center`).
    refine : bool, optional
        If True, apply iterative recentering.
    ellipse : bool, optional
        If True, use ellipse fitting for centering.
    save_path : str, optional
        Directory path to save images with contours and center overlay.

    Returns
    -------
    omega : numpy.ndarray
        Rotation rate in deg/day.
    time_axis : numpy.ndarray
        Time axis in hours.
    position_arc : numpy.ndarray
        Feature positions in arcseconds.
    """
    # --- 1. Time ---
    print('Step 1: Calculate Time Steps...')
    durations = get_time_steps(file_array, directory)

    # --- 2. Predict tracks ---
    print('Step 2: Calculate Positions...')
    tracks = calculate_tracks(file_array, directory, x0, y0, durations)

     # --- 3. Header info ---
    print('Step 3: Get Header Info...')
    sun_cx, sun_cy, dx, dy, R_sun, L, B, times = extract_header_params(file_array, directory)

    solar_center = np.array([sun_cx[0], sun_cy[0]])
    pix_to_arc = dx[0]  # assumed constant

    # set up arrays
    center_positions = []
    frame_positions = []
    if refine:
        frame_positions_refined = []

    # --- 4. Loop over images ---
    print('Step 4: Crop Regions & Find Center Coordinates...')

    cropped_cube = np.zeros([len(file_array), frame_size, frame_size])

    for i, fname in enumerate(file_array):

        image = f.getdata(os.path.join(directory, fname))
        header = f.getheader(os.path.join(directory, fname))
        header = f.open(os.path.join(directory, fname))

        # FITS file is rotated 180degrees against observation data, flip to match:
        image_flip = np.rot90(image, 2)

        # --- crop region series  ---
        frame_pos, crop = track_region(image_flip, header, tracks, frame_size, i)

        # --- optional refinement ---
        if refine:
            # get global pixel location of frame for each given step
            hdr = header
            cx = hdr[1].header['CRPIX1']
            cy = hdr[1].header['CRPIX2']
            dx = hdr[1].header['CDELT1']
            dy = hdr[1].header['CDELT2']

            frame_pix_x = int(round(cx + frame_pos[0] / dx))
            frame_pix_y = int(round(cy + frame_pos[1] / dy))

            # recenter the image around a feature:
            new_x, new_y, accept = refine_centering(crop, frame_pix_x, frame_pix_y, frame_size, tracking_mode="largest")

            #convert back to arc_sec
            new_x_arc = (new_x - cx) * dx
            new_y_arc = (new_y - cy) * dy
            frame_pos_refined = np.array([new_x_arc, new_y_arc])

            frame_pos, crop = track_region(image_flip, header, frame_pos_refined, frame_size)

        cropped_cube[i] = crop #fill cube

        # --- find center ---
        if save_path is not None:
            name = 'center_'+str(i)
            new_path = save_path + name

        else:
            new_path = None

        cx, cy = center(crop, threshold=threshold, ellipse=ellipse, save_path=new_path)

        center_positions.append([cx, cy])
        frame_positions.append(frame_pos)

    # write final positions of center coordinates in crop and global frame position
    center_positions = np.array(center_positions)
    frame_positions = np.array(frame_positions)

    f.writeto(os.path.join(directory, 'cropped_cube_TEST.fits'), np.nan_to_num(cropped_cube), overwrite=True) #save cube


    # --- 5. Convert coordinates ---
    print('Step 5: Convert Coordinates...')
    position_pix, position_arc = position_conversion(center_positions, frame_positions, frame_size, solar_center, pix_to_arc, flip_x=1, flip_y=1)


    # --- 6. Rotation Rate ---
    print('Step 6: Calculate Carrington Rotation Rate...')
    omega, time_axis = differential_rotation_rate(position_arc, durations, R_sun, L, B)
    print(f'Rotation Rate: {omega}')

    print('DONE!')

    return omega, time_axis, position_arc








