import cv2
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import astropy.io.fits as f
from astropy.time import Time, TimeDelta
import astropy.units as u
import datetime
from sunpy.net import Fido, attrs as a
import os
import astropy.units as u
from astropy.coordinates import SkyCoord
import sunpy.map
from sunpy.coordinates import RotatedSunFrame
from scipy.optimize import curve_fit

def check_path(path, folder):
    '''
    checks if a specific directionary ('folder') in a path already exists and creates it, if not.

    Parameters
    ----------
    path:   string
            gives the path in which the wanted directionary should be in

    folder: string
            the name of the folder/directionary, that should be checked for its existence or be created, if it does not already exist
    Returns
    -------
    folder_path
        string, a combined path to the folder that was checked/created

    :Authors:
        Emily Joe Loessnitz (2024)
    '''
    folder_path = os.path.join(path, folder)
    existence = os.path.exists(folder_path)
    if not existence:
        os.makedirs(folder_path)
    return folder_path;

def get_header_info(i):
    """
    opens a .FITS file and reads the contant of its header. Gives out a multitude of parameters

    Parameters
    ----------
    i  :    int
            determines the positition of the file in the download folder and therefor its name in the name-list if i is 0, for example, the first .FITS file would be opened and its header would be read

    Returns
    -------
    Center_pix_x
            x-coordinate of the solar-disk center in the image, in pixel
    Center_pix_y
            y-coordinate of the solar-disk center in the image, in pixel
    delta_x
            conversion-rate from pixels to arcsec in x-direction
    delta_y
            conversion-rate from pixels to arcsec in y-direction
    R_sun
            observed solar radius [in arcsec]
    R_sun_pix
            observed solar radius [in pixel] (this value is converted from R_sun with delta_x)
    L_Carrington
            Carrington Longitude [in degree]
    B_Carrington
            Carrington Latitude [in degree]
    time_stamp
            exact timestamp of the .FITS file in UT (Time)

    :Authors:
        Emily Joe Loessnitz (2024)
    """
    #set up header from .FITS file
    hdr = f.getheader(FITS_path+ '/' +str(file_name_array[i]))
    hdr = f.open(FITS_path+ '/' +str(file_name_array[i]))

    # get pixels that align with center of the sun
    Center_pix_x = hdr[1].header['CRPIX1']       #center x-pixel for x=0
    Center_pix_y = hdr[1].header['CRPIX2']       #center y-pixel for y=0

    # get arcsec/pixel for x- and y-direction
    delta_x = hdr[1].header['CDELT1']       #x-direction
    delta_y = hdr[1].header['CDELT2']       #y-direction

    R_sun = hdr[1].header['RSUN_OBS']  # read out the observerd angular radius of the sun [arcsec]
    R_sun_pix = int(R_sun/delta_x) # conversion to pixels

    #get the paramters for Carrington Rotation
    L_Carrington = hdr[1].header['CRLN_OBS']  #in [deg]
    B_Carrington = hdr[1].header['CRLT_OBS']  #in [deg]

    time_stamp = hdr[1].header['DATE-OBS'] # exact time of observation

    return Center_pix_x, Center_pix_y, delta_x, delta_y, R_sun, R_sun_pix, L_Carrington, B_Carrington, time_stamp;

def get_time_steps(length):
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
    durations = [0]  #set 0 as first entry, since the is no previous entry

    # for every file after the first:
    for i in range(1, length):
        hdr = f.getheader(FITS_path+ '/' +str(file_name_array[i]))
        hdr = f.open(FITS_path+ '/' +str(file_name_array[i]))
        time_i = hdr[1].header['DATE-OBS']

        # get the time_stamp from the previous image
        hdr = f.getheader(FITS_path+ '/' +str(file_name_array[i-1]))
        hdr = f.open(FITS_path+ '/' +str(file_name_array[i-1]))
        time_i0 = hdr[1].header['DATE-OBS']

        #divide time from image with the one from the previous file
        time_step = (Time(time_i)-Time(time_i0)).value
        durations.append(time_step) # add the caculated time difference to the array

    durations = (durations)*u.hour*24 #make sure durations are in the correct units

    return durations ;

def get_frame_locations(file_name_array, x_start_position, y_start_position, durations):
    """
    calculates the approximate position the frame need to be in by using a sunpy module to differentially rotate from a giving starting position
    ----------
    file_name_array: array (with str)
            array containing all the FITS-file names in chronological order
    x_start_position: int
            x-coordinate of the starting position [in arcsec] where the spot is first visible on the solar disc
    y_start_position: int
            y-coordinate of the starting position [in arcsec] where the spot is first visible on the solar disc
    durations: array
            array containing the time difference between to images in series
    Returns
    -------
    durations
            array of the amount of time between the corresponding two images in the name-list. This array corresponds to the file_name_array
    :Authors:
        Emily Joe Loessnitz (2024)
    """

    print('Calculating frame locations...')
    ##############################################################################
    # First, load an observation and define a coordinate in its coordinate
    # frame (here, helioprojective Cartesian).  The appropriate rate of rotation
    # is determined from the heliographic latitude of the coordinate.

    #make sunpy Map from the information of the first FITS file of the series
    name  = file_name_array[0]
    data  = f.getdata(FITS_path+'/'+name)
    header= f.getheader(FITS_path+ '/' +name)
    header['cunit1'] = 'arcsec'
    header['cunit2'] = 'arcsec'
    aiamap = sunpy.map.Map(data, header)

    #defines start point which we can rotate
    point = SkyCoord((x_start_position)*u.arcsec, (y_start_position)*u.arcsec, frame=aiamap.coordinate_frame)

    ##############################################################################
    # We can differentially rotate this coordinate by using
    # `~sunpy.coordinates.metaframes.RotatedSunFrame` with an array of observation
    # times

    diffrot_point = SkyCoord(RotatedSunFrame(base=point, duration=np.cumsum(durations)))

    ##############################################################################
    # To see what this coordinate looks like in "real" helioprojective
    # Cartesian coordinates, we can transform it back to the original frame.
    # Since these coordinates are represented in the original frame, they will not
    # account for the changing position of the observer over this same time range.

    transformed_diffrot_point = diffrot_point.transform_to(aiamap.coordinate_frame).to_string(unit='arcsec')
    # we only care about the change in x-direction right now:
    frame_location_arcsec = np.array([float(item.split()[0]) for item in transformed_diffrot_point])
    # relate this to the center coordinate to get the x_coordinate in pixels:
    sunpy_location = Center_pix_x - (frame_location_arcsec / delta_x)

    print('done!')

    return sunpy_location ;

def track_region(series, sunpy_location, frame_size, file_name_array, y_start_position, Center_pix_y, delta_y):
    """
    this function now finally crops out the sunspot and surrounding area and saves the crop-outs as pngs in a spereate folder. additionally it saves them to a data cube, which will also be saved in the output folder. In the process the coordinates (related to the entire solar disc) of the center of each cropped images will be saved. The y-coordinate stays constant and the respective x-coordinate gets saved into an array

    Parameters
    ------------
    series :    string
            name of the SDO data series
    sunpy_location: array
            contains the frame-locations calculated by the sunpy-module
    frame_size: int
            lenghth [in pixels] of the side of a quadratic frame that follows the active region (f.e. if set to 100 the image would be 100x100 pixels)
    file_name_array: array (with str)
            array containing all the FITS-file names in chronological order
    y_start_position: int
            y-coordinate of the starting position [in arcsec] where the spot is first visible on the solar disc
    Center_pix_y: float
            reference central coordinates of the solar disc [in pixels](extracted from the header of the FITS files)
    delta_y: float
            conversion from pixels to arcsec in y-direction (extracted from the header of the FITS files)
    Returns
    -------
    frame_location
            array of the final locations (pixels in x-direction) of the frame-center
    y_location
            float of the final location of the frame-center in y-direction. Since the change in y-direction is assumed to be small against the change in x-direction, this does not need to be an array and is instead constant for the series
    cropped_cube
            besides saving the cropped images as pngs in the output folder, the data gets also saved in a data cube [FITs file]

    :Authors:
        Emily Joe Loessnitz (2024)
    """

    #set up an empty data cube of the right dimensions to store the cropped FITS-files
    cropped_cube = np.zeros([length, frame_size, frame_size])
    #set up empty array for the final coordinates of the frame
    frame_location = []

    #for all images in the series:
    for i in range(length):
        name = file_name_array[i]
        pic = f.getdata(FITS_path+'/'+name) #load data

        print(f'tracking the active region... current image: {i}')

        y_location = int(Center_pix_y - (y_start_position/ delta_y)) #starting position in y (stays constant for full rotation) gets converted into pixel
        x_location = int(sunpy_location[i]) #read out corresponding entry from the locations obtained by the differential rotation (sunpy module)

        #close to the edge of the visible solar disc it can happen that the loctaion of the frame is too close to the borders of the entire image, so that the frame-size would extend over it. For this reason the coordinates will now be checked, so if they would be too close for the choosen frame-size to fit, it will overwrite the frame location with the one closest possible one:
        if x_location<(frame_size/2):           #if its too close to the left border
            x_location_failsave=int(frame_size/2)
            pic_crop = pic[y_location-(frame_size//2):y_location+(frame_size//2), x_location_failsave-(frame_size//2):x_location_failsave+(frame_size//2)]
            frame_location.append(x_location_failsave) #overwrite old location, add new one

        elif x_location>(4096-(frame_size/2)):  #if its too close to the right border
            x_location_failsave=int(4096-(frame_size/2))
            pic_crop = pic[y_location-(frame_size//2):y_location+(frame_size//2), x_location_failsave-(frame_size//2):x_location_failsave+(frame_size//2)]
            frame_location.append(x_location_failsave) #overwrite old location, add new one

        else:  #if there is no problem, the frame-loctation from the sunpy module will be copied
            pic_crop = pic[y_location-(frame_size//2):y_location+(frame_size//2), x_location-(frame_size//2):x_location+(frame_size//2)]
            frame_location.append(x_location)

        pic_rot = np.rot90(pic_crop, 2)         #rotate image the right way
        pic_rim = np.nan_to_num(pic_rot, nan=1) #fill nan's outside the solar disc with 1's
        #(this last step was neccessary, because otherwise the nans would appear dark on png images made from the data, which would mess with the center-finding process)

        #save cropped images as pngs:
        plt.imsave(png_path+'/'+file_name_array[i]+'.png', pic_rim, cmap=plt.cm.gray, origin='lower')

        cropped_cube[i] = pic_rim #fill cube

    f.writeto(new_path+ '/cropped_cube.fits', np.nan_to_num(cropped_cube), overwrite=True) #save cube

    plt.close()
    print('tracking complete!')

    return frame_location, y_location, cropped_cube ;

def make_grayscale(filename):
    """
    this function loads the pngs from a sunspot-rotation series, aligns the orientation properly and converts into a grayscale image. Finally a guassian blur is applied to smooth the image, so that fine features in the granulation appear lighter in contrast to the sunspots

    Parameters
    ------------
    series :    string
            name of the SDO data series
    Returns
    -------
    grayscaleImage
            grayscale version of input image (with axes flipped)
    blur
            grayscaleImage with applied Gaussian Blur

    :Authors:
        Emily Joe Loessnitz (AIP 2023)
    """
    #read pngs of the sunspot-series from the corresponding folder
    inputImage = cv2.imread(png_path+'/'+filename+'.png')
    flip = inputImage[::-1,:,:] # revise height in (height, width, channel)
    #(this flip is because of the way the FITS files were open, so that the axis are reversed)
    # Convert BGR to grayscale:
    grayscaleImage = cv2.cvtColor(flip, cv2.COLOR_BGR2GRAY)
    #blur the image to make granulation appear lighter in contrast to sunspots
    blur = cv2.GaussianBlur(grayscaleImage,(7,7),0)

    return grayscaleImage, blur ;

def make_binary(image, threshold):

    """
    creates an binary image of an (grayscale) image by sorting every pixel depending on if it surpasses the selected threshold

    Parameters
    ----------
    image:  str
            selected image which should be turned into a binary image

    threshold : int
                the threshold value between 0 (black) and 255 (white)

    Returns
    -------
    binaryImage
        Sorted array of data with renumbered 'id' values.

    :Authors:
        Emily Joe Loessnitz (AIP 2024)
    """
    #apply thresholding
    threshValue, binaryImage = cv2.threshold(image,threshold,255,cv2.THRESH_BINARY_INV)

    plt.imshow(binaryImage, cmap ='gray' )
    plt.gca().invert_yaxis() #making sure the axis are converted

    return binaryImage ;

def perform_opening(image):
    """
    applies erosion/structuring element to binary image to smooth out noises. The scale of this erosion depends on the chosen kernel size which results in a rectengular structuring element. Using the cv / cv2 library an opening of the image will be performed afterwards.

    Parameters
    ----------
    image :     str
                name of image-file

    Returns
    -------
    image :
        eroded and opened image (binary, if input was binary)

    :Authors:
        Emily Joe Loessnitz (AIP 2024)
    """
    # Get the structuring element:
    kernel = np.ones((3, 3), np.uint8)  #size of the kernel
    opIterations = 3 #number of iterations
    # Perform opening:
    openingImage = cv2.morphologyEx(image, cv2.MORPH_OPEN, kernel, None, None, opIterations, cv2.BORDER_REFLECT101)

    return openingImage ;

def get_centroid(image, overlay_image, SaveImage):
    """
    finds the central moments of a figure in the (binary) input image by first selecting the biggest contour (which would be the sunspot) and then find its center via the cv2.moments module. The center can be drawn into the original image (overlay-image) and saved, if wanted

    Parameters
    ----------
    image : array (?)
        name of the image file
    overlay_image : array(?)
        image where the center shall be drawn in
    SaveImage :  boolean
        'True' if images with center drawn in should be saved as pngs, 'False' if not
    Returns
    -------
    Cx_centroid
        center x-coordinate of the sunspot/figure [pixel] relative to image
    Cy_centroid
        center y-coordinate of the sunspot/figure [pixel] relative to image
    :Authors:
        Emily Joe Loessnitz (AIP 2024)
    """
    mask = np.zeros_like(image) # make a mask from the binary image
    # find the countours of the central figure (=sunspot/penumbra/umbra)
    contours,_ = cv2.findContours(image, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    contours    = [c for c in contours if cv2.contourArea(c) > 100]
    # only consider contours that have an area [pixels] over 100 to exclude noise or fine features in granulation
    # we assume the sunspot is the biggest object in frame, so we only consider the largest found contour
    biggest_contour = max(contours, key = cv2.contourArea)
    cv2.drawContours(mask, [biggest_contour], -1, 255, thickness=cv2.FILLED) #draw in contours to check

    # Calculate the moments
    imageMoments = cv2.moments(mask)
    # Compute centroid
    cx = int(imageMoments['m10']/imageMoments['m00'])
    cy = int(imageMoments['m01']/imageMoments['m00'])

    Cx_centroid = cx  #center x-coordinate
    Cy_centroid = cy  #center x-coordinate

    # draw the center into the original sunspot image and save those as pngs in the corresponding folder
    if SaveImage:
        # Draw centroid onto the overlay image (BGR):
        bgrImage = cv2.cvtColor(overlay_image, cv2.COLOR_GRAY2BGR)
        bgrImage = cv2.line(bgrImage, (cx,cy), (cx,cy), (0,255,0), 10) #set color and thickness of center
        plt.imshow(bgrImage, cmap = 'gray')
        plt.gca().invert_yaxis()            # flip axis to be the right way around
        plt.savefig(centroid_path+'/'+file_name_array[i]+"_central_moment"+'.png') # save image (with centroid and contours drawn in)
    plt.close()

    return Cx_centroid, Cy_centroid ;

def get_ellipseCenter(image, overlay_image, SaveImage):
    """
    finds the central moments of a figure in the (binary) input image by first dialating the figure slighlty, selecting the biggest contour (which would be the sunspot) and then find its center by fitting an ellipse of minimal area to it.  The center can be drawn into the original image (overlay-image) and saved if wanted

    Parameters
    ----------
    image : array (?)
        name of the image file
    overlay_image : array(?)
        image where the center shall be drawn in
    SaveImage :  boolean
        'True' if images with center drawn in should be saved as pngs, 'False' if not
    Returns
    -------
    Cx_ellipse
        center x-coordinate of the sunspot/figure [pixel] relative to image
    Cy_ellipse
        center y-coordinate of the sunspot/figure [pixel] relative to image
    :Authors:
        Emily Joe Loessnitz (AIP 2024)
    """

    #testing different methode, using countours to lay a circle of minimal area around the sunspot
    thresh = cv2.morphologyEx(image, cv2.MORPH_DILATE, np.ones((1, 1)))
    plt.imshow(thresh, cmap ='gray' )
    plt.gca().invert_yaxis()

    contours, _ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    contours    = [c for c in contours if cv2.contourArea(c) > 100]
    # only consider contours that have an area [pixels] over 100 to exclude noise or fine features in granulation
    # we assume the sunspot is the biggest object in frame, so we only consider the largest found contour
    biggest_contour = max(contours, key = cv2.contourArea)

    #find radius to minimalize circle area
    ellipse = cv2.fitEllipse(biggest_contour)
    (x, y), (w, h), ang = cv2.fitEllipse(biggest_contour) # fit ellipse and save the parameters

    Cx_ellipse = (int(x))    #center x-coordinate
    Cy_ellipse = (int(y))    #center x-coordinate


    # draw the center into the original sunspot image and save those as pngs in the corresponding folder
    if SaveImage:
        output = overlay_image.copy() # the image, where the center and the ellipse should be drawn on
        #draw ellipse
        cv2.ellipse(output, ellipse, (0, 255, 0), 2)
        #draw center
        cv2.circle(output, (int(x), int(y)), 3, (0, 255, 0), -1)
        plt.imshow(cv2.cvtColor(output, cv2.COLOR_BGR2RGB))
        plt.gca().invert_yaxis()
        cv2.drawContours(output, contours, -1, (0,255,0), 3),
        plt.savefig(ellipse_path+'/'+file_name_array[i]+"_ellipse_center"+'.png') #save final pictures into new folder (as .png)

    plt.close()

    return Cx_ellipse, Cy_ellipse ;

def translate_pixel_location(center_coordinate, frame_coordinate, frame_size, SolarDiskCenter, PixToArc):
    """
    this function first calculates the position of the sunspot-center in a given frame in regard to the full image, by adding

    Parameters
    ----------
    center_coordinate: float
        position (in one direction) calculated to be the center of the sunspot in the cropped image (in frame) [pixel]
    frame_coordinate: int
        position of the frame center (on the full image) for this image of the series [pixel]
    frame_size: int
        lenghth [in pixels] of the side of a quadratic frame which was set before
    SolarDiskCenter: float
        position (in same direction as center_coordinate) of the center of the solar disk [pixel] (this information can be obtained from the header)
    PixToArc: float
        conversion from pixel to arcsec (in that direction; can also be obtained from the header)
    Returns
    -------
    fullimage_postition
        position of sunspot center [in pixels] on the full image in chosen direction (x or y) (reference point = bottem left of full solar disk image)
    fullimage_position_arcsec
        position of sunspot center [in arcsec] on the full image in chosen direction (x or y) (reference point = solar disc center)
    :Authors:
        Emily Joe Loessnitz (AIP 2024)
    """
    # right now we only know the center of a sunspot at each point in the series in regards to the cropped image, meaning the coordinates obtained during the centroid or ellipse process describe the location of the center in the cropped image. In order to translate this to the actual solar disc, we use the fact that we know the position of each cropped image (its center aka the frame-location) and saved those in an array (in x-direction) or have a constant position (in y-direction). We can calculate the actual position on the solar disc as follows:
    fullimage_position = frame_coordinate + frame_size/2 - center_coordinate # in [pixel]
    # this position can now be converted into arcsec by first relating to the center of the solar disc and then converting pixels to arcsec, both need information from the FITS header:
    fullimage_position_arcsec = (-1) *(fullimage_position - SolarDiskCenter) * PixToArc

    # it is important to note that this function only converts for one direction (x or y) at a time and must be run twice with different input parameters for the whole set

    return fullimage_position, fullimage_position_arcsec ;

def rotation_speed(col, durations):
    """
    calculate the rotation speed at each step of the series during a rotation

    Parameters
    ----------
    col : int
        number of the collum corresponding to the x-position from the method of choice
    durations: array
        contains the times between two consecutive images in the series

    Returns
    -------
    omega_with_Carrington
        rotation velocity with the Carrington Rotation
    time_axis
        timestemps for the data-points for plot

    :Authors:
        Emily Joe Loessnitz (AIP 2024)
    """

    #set up empty arrays to save data in for later
    l_deg_array = []    # longitude in degree
    time_array = []     # time-steps

    # open and read the cooridante list
    with open(new_path+ '/' + str(NOAA) + '_' + "coordinates_list_full.txt", "r") as coordinates:
        next(coordinates) #skip first line (header)
        for line in coordinates:

            currentline = line.split(",")   # split lines at , to make data accessable
            i = int(float(currentline[0]))  # the number in the series
            x = float(currentline[col])     # x-position [arcsec]
            y = float(currentline[(col+1)]) # y-position [arcsec]

            Center_pix_x, Center_pix_y, delta_x, delta_y, R_sun, R_sun_pix, L_Carrington, B_Carrington, time_stamp = get_header_info(i) # information from header for each fits-file (i)

            # Carrington Longitude (L) and Latitude (B) from the header need to be converted:
            L_C = L_Carrington * ( np.pi/180)       #convert to [rad]
            B_C = B_Carrington * ( np.pi/180)       #convert to [rad]

            # we need position of sunspot-center in heliocentric coordinates to converted to carrington system:
            SC = np.sqrt(x**2 + y**2)               #distance between spot (S) and sun-center (C)[arcsec]
            SC_rad = SC * np.pi /(180 * 3600)       #converting to [rad]
            SC_vector = np.array([x, y])            #vector to the spot

            # we need orientation, to assign a sign to both sides of the solar disc center
            N  = np.array([0, R_sun])               #North Vector
            cross_product = np.cross(N,SC_vector)   #vector product to get sign of alpha (angle)
            alpha = np.arccos(y / SC)               #angle between North and SC
            if cross_product < 0:  # if crossproduct negative then the angle is per definition negative too (positive from north eastwards)
                    alpha2 = alpha * (-1)
            else:
                    alpha2 = alpha

            #now we have the heliocentric coordinate system:
            sigma = np.arcsin(SC / R_sun) - SC_rad  #heliocentric angle sigma [rad]
            # heliographic latitude:
            b = np.arcsin( np.cos(sigma) * np.sin(B_C) + np.sin(sigma) * np.cos(B_C)* np.cos(alpha2))
            # heliographic longitude:
            longitude_central = np.arcsin( np.sin(alpha2) * np.sin(sigma) * (1/np.cos(b)) ) #distance from central meridian
            # in the carrington system the longitude is given in respect to the carrington longitude L_C which rotates with the system:
            longitude = L_C - longitude_central
            longitude_deg = longitude * (180/np.pi) # convert to degree
            l_deg_array.append(longitude_deg)       # append to array

    # the longitude in the carrington system was calculated and saved for every point in the series
    car_longitude = np.asarray(l_deg_array)

    # for the rotation speed we need to know how far the spot moved between to images:
    longitude_diff = np.diff(car_longitude) # array with distance moved between images
    mask = (longitude_diff<1)               # mask (=ignore) all points where the difference is smaller than 1 for this case can only happen if the start of a new carringtion rotation started during the observation (in that case the difference will be negative!)
    masked_longitude_diff = ma.masked_where(mask == False, longitude_diff)  #apply mask

    # to calculate speed we also need to know in which time the spot moved for those distances:
    masked_durations = ma.masked_where(mask== False, durations)         # apply the mask
    masked_time = ma.masked_where(mask== False, np.cumsum(durations))   # cumulative sum to set up time axis
    time_axis = (masked_time) + (durations/2)   # assign the velocity to the time between two images in a series, since the velocity was calcultated from the distance and time between those

    # calculate the rotational velcity by taking the distance moved (delta x or in this case longitude_diff) and divide by the time the movement took (in h) and multiply by 24 to get [deg/day]
    omega = masked_longitude_diff *24 / masked_durations # rotational velocity omega in [deg/day]
    omega_with_Carrington = omega + 14.184               # add rotation speed of the Carrington System (14.184 deg/day)

    return omega_with_Carrington, time_axis ;

##########################################################################################
#       PLEASE SET BEFORE RUNNING:
##########################################################################################
# set email for conformation from the SDO Database
myemail = 'loessnitz@uni-potsdam.de'

# set directionary where the output should appear
mydir = '/work2/loessnitz/'

# NOAA of the to be observed active region
NOAA = 12218

# Latitude [deg]
latitude = 16

# Spefify the cadance between two pictures [hours]
dt = 6

# set the size of the region that should be cropped around the active region. It will define the lenghth [in pixels] of the side of a quadratic frame (f.e. if set to 100 the image would be 100x100 pixels)
frame_size = 400

# Specifiy threshold value (umbra/penumbra)
penumbra_threshold = 175
umbra_threshold = 83

#if different JSOC data series is wanted, please edit:
No_LD_Series = 'hmi.Ic_noLimbDark_720s'

# please select if you want images (png) showing the center and contours of sunspots saved ('True') or not ('False')
SaveImage = True

##########################################################################################
#       PLEASE GIVE SUNSPOT A RATING:
##########################################################################################
#some alpha spots are even more ALPHA than others. If there is f.e. missing data, some sidespots, general weirdness going on consider giving a lower score
#score between 0(unusable) and 3 (the ideal alpha-spot):
weight = 3
####################################################################################################

####################### FOLDER- & OUTPUT-STRUCTURE ########################################
# set up new directionary for output for this spot
main_dir = 'NOAA_' + str(NOAA) + '_dt_' + str(dt) + 'h'
new_path = check_path(mydir, main_dir)

#inside the new directionary, set up sub-dirs for images and different methods (optinal?)
FITS_path    = check_path(new_path, "FITS_files")
png_path     = check_path(new_path, "png_images")
centroid_path= check_path(new_path, "central_moments")
ellipse_path = check_path(new_path, "ellipse_center")

##########################################################################################

# read out the starting coordinates of the active region from the file created with the clicker
checked_coordindates = np.loadtxt(new_path+'/'+"coordinate_checked.txt", dtype=str, comments="#", delimiter=",", unpack=False)
x_start_position = float(checked_coordindates[3])+5
y_start_position = float(checked_coordindates[4])+15

# set up array containing the file-names of all the files from the series
file_name_array = np.loadtxt(new_path+'/'+'NOAA_'+str(NOAA)+"file_names_list.txt", dtype=str, comments="#", delimiter=",", unpack=False)
length = len(file_name_array)   # the length of the filename-array aka how many FITS-files there are

# create array with the how much time passed between two images in series:
durations = get_time_steps(length)

# obtain information stored in the header of the first .FITS file for reference
Center_pix_x, Center_pix_y, delta_x, delta_y, R_sun, R_sun_pix, L_Carrington, B_Carrington, time_stamp  = get_header_info(0)

# get the location (x-position in pixel) of the frame for each image
sunpy_location = get_frame_locations(file_name_array, x_start_position, y_start_position, durations)

#crop images, save them as pngs and as a fits cube
frame_location, y_location, cropped_cube = track_region(No_LD_Series, sunpy_location, frame_size, file_name_array, y_start_position, Center_pix_y, delta_y)

# set up new txt-list to write the coordinates to:
coordinates_list = open(new_path+'/'+ str(NOAA)+'_'+"coordinates_list_full.txt", "w+")
#set up a header for the table: (repeat same layout like original center coordinate table)
coordinates_list.write ("image#, x_centroid[px], y_centroid[px], x_centroid[arcsec], y_centroid[arcsec], x_ellipse[px], y_ellipse[px], x_ellipse[arcsec], y_ellipse[arcsec]\n")

for i in range(length):             # for every file
    filename = file_name_array[i]   # read filename to load image
    number = str(i)                 # number of current step in the series

    print(f'calculating center... current image: {i}')

    grayscaleImage, blur = make_grayscale(filename)     # load image, turn to grayscale and apply gaussian blur
    binaryImage = make_binary(blur, penumbra_threshold) # make image binary according to set threshold
    overlay_image = grayscaleImage                      # set up image to (optinally) draw center coordinates in

    # Find center coordinates with two methods:
    Cx_centroid, Cy_centroid = get_centroid(binaryImage, overlay_image, SaveImage)    # Central Moments Method
    Cx_ellipse, Cy_ellipse = get_ellipseCenter(binaryImage, overlay_image, SaveImage) # Ellispe-Fitting Method
    plt.close()

    # Read out Header-Information for each FITS-file:
    Center_pix_x, Center_pix_y, delta_x, delta_y, R_sun, R_sun_pix, L_Carrington, B_Carrington, time_stamp= get_header_info(i)

    x_location = frame_location[i]  # corresponding coordinate of frame for this step (i)

    # transform AR-center-coordinates from pixels of cropped images to coordinates in arcsec in regards to the solar disc:
    x_centroid_fixed, x_cent_arcsec = translate_pixel_location(Cx_centroid, x_location, frame_size, Center_pix_x, delta_x)  # x-position for centroid method
    y_centroid_fixed, y_cent_arcsec = translate_pixel_location(Cy_centroid, y_location, frame_size, Center_pix_y, delta_y)  # y-position for centroid method
    x_ellipse_fixed, x_elli_arcsec = translate_pixel_location(Cx_ellipse, x_location, frame_size, Center_pix_x, delta_x)                # x-position for ellipse method
    y_ellipse_fixed, y_elli_arcsec = translate_pixel_location(Cy_ellipse, y_location, frame_size, Center_pix_y, delta_y)                # y-position for ellipse method

    # write final coordinates down into a txt-table:
    coordinates_list.write("{0},{1},{2},{3},{4},{5},{6},{7},{8}\n".format(number, round(x_centroid_fixed, 0), round(y_centroid_fixed, 0), x_cent_arcsec, y_cent_arcsec, round(x_ellipse_fixed, 0), round(y_ellipse_fixed, 0), x_elli_arcsec, y_elli_arcsec))

coordinates_list.close()    # close list
print('Done!')

# delete first entry of duration-array, because the first entry has no reference
durations2 = np.delete(durations, 0)

# calculate roational velocity of both methods:
omega_cent, time_axis_cent = rotation_speed(3, durations2)
omega_elli, time_axis_elli = rotation_speed(7, durations2)

############ STATISTICAL ANALYSES ###########################################################
# for centroid-method:
mean_cent = np.mean(omega_cent)     # mean rotation velocity (centroid method)
deviation_cent = np.std(omega_cent) # calculate standart deviation of rot. velocity
cut_cent = omega_cent[3:-3]         # cut of 3 data points from each side
mean_cut_cent = np.mean(cut_cent)   # mean rotational velocity (without the outer edges)
mean_cent_rounded = round(mean_cut_cent, 3) # round the mean for plot
deviation_cut_cent = np.std(cut_cent) # standart deviation without the outer edges
deviation_mean_cent = deviation_cut_cent/np.sqrt(cut_cent.size) # standart deviation of the mean (its statistical uncertainty)

# for ellipse-method:
mean_elli = np.mean(omega_elli)     # mean rotation velocity (ellipse method)
deviation_elli = np.std(omega_elli) # calculate standart deviation of rot. velocity
cut_elli = omega_elli[3:-3]         # cut of 3 data points from each side
mean_cut_elli = np.mean(cut_elli)   # mean rotational velocity (without the outer edges)
mean_elli_rounded = round(mean_cut_elli, 3) # round the mean for plot
deviation_cut_elli = np.std(cut_elli) # standart deviation without the outer edges
deviation_mean_elli = deviation_cut_elli/np.sqrt(cut_elli.size) # standart deviation of the mean (its statistical uncertainty)

################################# PLOT AND FINAL OUTPUT ###############################################
# set up plot:
plt.title('rotational velocity of sunspot NOAA'+str(NOAA))
#plt.xlabel("heliographic longitude")
plt.xlabel("time [h]")
plt.ylabel("rotation omega [deg/day]")
plt.ylim([8, 20])
# plot the rotational velocitys of both methods over time
plt.plot(time_axis_cent, omega_cent, 'b', linestyle="solid", marker='o', label='centroid')
plt.plot(time_axis_elli, omega_elli, 'g', linestyle="solid", marker='o', label='ellipse' )
#plot constant line for mean-value as a reference:
plt.axhline(y = mean_cut_cent, color = 'b', linestyle = 'dashed', label='mean omega (centroid) ='+ str(mean_cent_rounded)+' deg/day')
plt.axhline(y = mean_cut_elli, color = 'g', linestyle = 'dashed', label='mean omega (ellipse) ='+ str(mean_elli_rounded)+' deg/day')
leg = plt.legend()

plt.savefig(new_path+'/'+'rotation_plot_NOAA'+str(NOAA)+'.png') # save plot as png
plt.close()

# write results into txt-file to view in the output-folder
results = open(new_path+'/'+str(NOAA)+'_'+"results.txt", "w+")
results.write('NOAA: {0}\nSCORE: {1}\nCADANCE: {2}\nLATITUDE: {3}\nSTART_X: {4}\t#arcsec\nSTART_Y: {5}\t#arcsec\nTHRESH: {6}\nC_ROT_V: {7}\nC_STD: {8}\nC_MEAN_STD: {9}\nE_ROT_V: {10}\nE_STD: {11}\nE_MEAN_STD: {12}'.format(NOAA, weight, dt, latitude, x_start_position, y_start_position, penumbra_threshold, mean_cut_cent, deviation_cut_cent, deviation_mean_cent, mean_cut_elli, deviation_cut_elli, deviation_mean_elli))
results.close()
