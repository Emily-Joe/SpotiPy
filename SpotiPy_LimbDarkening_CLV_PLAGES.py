import cv2
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as f
from astropy.time import Time
import astropy.units as u
import datetime
from sunpy.net import Fido, attrs as a
from sunpy.time import parse_time
import os
from astropy.coordinates import SkyCoord
import sunpy.map
from sunpy.coordinates import RotatedSunFrame
from scipy.optimize import curve_fit

def check_path(path, folder):
    '''
    checks if a specific directionary ('folder') in a path already exists and creates it, if not./work2/loessnitz/

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

def check_file(path):
    if os.path.exists(path):
        os.remove(path)
        print('existing results deleted')

def downloadSDO(series, use_time_stamps, time_stamps, start_date, dt, days, download_path):

    """
    make a search query to the SDO Database to extract images of the sun in the desired timeframe.
    These are the images at --- and already correct the effects of limbdarkening.
    Checking the database might take some minutes.
    If the search-request is successful, an email will be send to your email account and the download will begin.
    The images will be downloaded as .FITS files into the directory.

    Parameters
    ----------
    t_start :  ?
            stets the start date and time in the format YYYY-MM-DDT00:00:00

    dt  :   float
            specifies the cadance between two images in hours

    days  : float?
            specifies the number of days over wich the observation shall take place

    download_path : str
            sets path and directionary for images to be downloaded to

    Returns
    -------
    numpy.ndarray
        Sorted array of data with renumbered 'id' values.

    :Authors:
        Emily Joe Loessnitz (2023)
    """
    if use_time_stamps==False:

        steps = int(days*24 / dt)

        print(f'Preparing to download {steps} images starting from {start_date}, spanning {days} days, with a {dt} hour cadence.')

        print('\n\n\nChecking Database...')


        # Calculate start and end times for the search
        start_time = start_date #+ datetime.timedelta(hours=i)
        end_time = start_time + datetime.timedelta(hours=dt*steps) #only select the right images

        # Make the search query
        res = Fido.search(a.Time(start_time.iso, end_time.iso), a.jsoc.Series(series), a.Sample(dt*u.hour) , a.jsoc.Notify(myemail))


        print(res.show('TELESCOP', 'INSTRUME', 'T_OBS'))

        #i would like to count the number of images here and double check if these equal the number of calculated steps, because sometimes there is data missing and that messes with the calculations later

        print('Does this look correct?')
        resp = input('Press [y] to download and [n] to exit\n')
        if resp == 'y':
            print('Downloading ...')
            downloaded_files = Fido.fetch(res, path=download_path+'/{file}')
        else:
            print('END')
            pass
        print('END')

    elif use_time_stamps== True:

        array_length = len(time_stamps)
        print(f'Preparing to download {array_length} images according to the provided time-stamp list...')

        # Convert string time stamps to SunPy time objects
        times = [parse_time(ts) for ts in time_stamps]

        # Create an empty list to hold individual search results
        search_results = []

        print('\n\n\nChecking Database...')

        # Perform individual searches for each time stamp
        for time in times:
            result = Fido.search(a.Time(time, time + 1 * u.minute), a.jsoc.Series(series), a.Sample(720 * u.second), a.jsoc.Notify(myemail))
            search_results.append(result)

        # Display the combined results
        for res in search_results:
            print(res)

        #i would like to count the number of images here and double check if these equal the number of calculated steps, because sometimes there is data missing and that messes with the calculations later

        print('Does this look correct?')
        resp = input('Press [y] to download and [n] to exit\n')
        if resp == 'y':
            print('Downloading ...')
            downloaded_files = Fido.fetch(*search_results, path=download_path+'/{file}')
        else:
            print('END')
            pass
        print('END')


    else:
        print('Specify if you want to download in a cadance (dt) or if you want to use an array of specific time stamps. Please set use_time_stamps to either True or False and set a cadance if False')


    # directory/folder path
    dir_path = download_path
    dirs = sorted(os.listdir(dir_path))

    # list to store files
    file_names = []

    if series== 'hmi.Ic_noLimbDark_720s':
         file_names_list = open(new_path + '/' +'NOAA_' + str(NOAA) + "file_names_list.txt", "w+")
    else:
        #file_names_list = open("file_names_list.txt", "w+")
        file_names_list = open(new_path + '/' +'NOAA_' + str(NOAA) + '_LimbDark_' +  "file_names_list.txt", "w+")

    # Iterate directory=
    for file_path in dirs:
        name = str(file_path)
        file_names_list.write(name + '\n')
        file_names.append(name)

    file_names_list.close()

    file_name_array = np.asarray(file_names)


    print('Download complete!')

    return FITS_path, file_name_array ;

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
    sunpy_location

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

    print(f'tracking the active region in dataset: {series}...')
    #set up an empty data cube of the right dimensions to store the cropped FITS-files
    cropped_cube = np.zeros([length, frame_size, frame_size])
    #set up empty array for the final coordinates of the frame
    frame_location = []

    #for all images in the series:
    for i in range(length):
        if series == 'hmi.Ic_noLimbDark_720s':
            name = file_name_array[i]
            pic = f.getdata(FITS_path+'/'+name)
        else:
            name = file_name_array_LimbDark[i]
            pic = f.getdata(LD_FITS_path+'/'+name)

        #print(f'tracking the active region... current image: {i}')

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

        if series == 'hmi.Ic_noLimbDark_720s' :
            #making a .png file for every .fits file in the list, saviing in another folder (named simply as i rn)
            plt.imsave(png_path+'/'+file_name_array[i]+'.png', pic_rim, cmap=plt.cm.gray, origin='lower')

        #make a data cube of all telect file which should be turned into a binary imagehe images
        cropped_cube[i] = pic_rim#fill cube

    if series == 'hmi.Ic_noLimbDark_720s' :
        f.writeto(new_path+ '/cropped_cube.fits', np.nan_to_num(cropped_cube), overwrite=True)
    else:
        f.writeto(new_path+ '/cropped_cube_limbdark.fits', np.nan_to_num(cropped_cube), overwrite=True)

    plt.close()

    print('tracking complete!')

    return frame_location, y_location, cropped_cube ;

def mask_sunspot(filename, threshold):
    """
    This funciton creates a binary mask of the (full) sunspot. It loads the pngs from a sunspot series, aligns the orientation properly and converts into a grayscale image. Finally a guassian blur is applied to smooth the image, so that fine features in the granulation appear lighter in contrast to the sunspots. Via a pre-determined threshold value the image will then be turned into a binary mask of the sunspot

    Parameters
    ------------
    filename :  str
                selected image which should be turned into a binary image
    threshold : int
                the threshold value between 0 (black) and 255 (white)
    Returns
    -------
    fullspot
        binary mask, sorted array of data with renumbered 'id' values.

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

    # apply thresholding
    threshValue, binaryImage = cv2.threshold(blur, threshold,255,cv2.THRESH_BINARY_INV)
    plt.imshow(binaryImage, cmap ='gray' )
    plt.gca().invert_yaxis()

    fullspot = binaryImage

    return fullspot;

def mask_umbra(filename, threshold):
    """
    This funciton creates a binary mask of the isolated umbra of a sunspot. It loads the pngs from a sunspot series, aligns the orientation properly and converts into a grayscale image. Finally a guassian blur is applied to smooth the image, so that fine features in the granulation appear lighter in contrast to the sunspots. Via a pre-determined threshold value the image will then be turned into a binary mask of the umbra, which is then. additionally, a second, slightly dialated mask of the umbra will be created, which can later be subtracted from the fullspot to create a isolated penumbra mask, with no overlap of umbra and penumbra.

    Parameters
    ------------
    filename :  str
                selected image which should be turned into a binary image
    threshold : int
                the threshold value between 0 (black) and 255 (white)
    Returns
    -------
    umbra
        binary mask, sorted array of data with renumbered 'id' values.
    umbra_dilated
        a binary mask dialated with a kernel, sorted array of data with renumbered 'id' values.
    :Authors:
        Emily Joe Loessnitz (AIP 2023)
    """
    inputImage = cv2.imread(png_path+'/'+filename+'.png')

    flip = inputImage[::-1,:,:] # revise height in (height, width, channel)
    #(this flip is because of the way the FITS files were open, so that the axis are reversed)

    # Convert BGR to grayscale:
    grayscaleImage = cv2.cvtColor(flip, cv2.COLOR_BGR2GRAY)

    #blur the image to make granulation appear lighter in contrast to sunspots
    blur = cv2.GaussianBlur(grayscaleImage,(7,7),0)

    # apply thresholding
    threshValue, binaryImage = cv2.threshold(blur,threshold,255,cv2.THRESH_BINARY_INV)
    plt.imshow(binaryImage, cmap ='gray' )
    plt.gca().invert_yaxis()

    # set up kernels for image manipulation
    kernel_erode = np.ones((3, 3), np.uint8)    #slight erosion for normal umbra
    kernel_dialation = np.ones((7, 7), np.uint8)#dialation for later use

    umbra = cv2.erode(binaryImage, kernel_erode ,iterations = 1)            # create a minimal (eroded) mask of the umbra
    umbra_dilated = cv2.dilate(binaryImage, kernel_dialation, iterations=1) # create a maximal (dialated) mask of the umbra, which can then be subtracted from the fullspot, so that there is no overlap

    return umbra, umbra_dilated;

def mask_penumbra(filename, umbra_dilated, threshold):
    """
    This funciton creates a binary mask of the penumbra. It loads the pngs from a sunspot series, aligns the orientation properly and converts into a grayscale image. Finally a guassian blur is applied to smooth the image, so that fine features in the granulation appear lighter in contrast to the sunspots. Via a pre-determined threshold value the image will then be turned into a binary mask of the fullspot and will then be subtracted by the dialated umbra-mask, leaving only the isolated penumbra mask

    Parameters
    ------------
    filename :  str
                selected image which should be turned into a binary image
    umbra_dilated: array
                 a binary mask dialated with a kernel, sorted array of data with renumbered 'id' values.
    threshold : int
                the threshold value between 0 (black) and 255 (white)
    Returns
    -------
    penumbra
        binary mask, sorted array of data with renumbered 'id' values.

    :Authors:
        Emily Joe Loessnitz (AIP 2023)
    """
    inputImage = cv2.imread(png_path+'/'+filename+'.png')
    flip = inputImage[::-1,:,:] # revise height in (height, width, channel)
    #(this flip is because of the way the FITS files were open, so that the axis are reversed)

    # Convert BGR to grayscale:
    grayscaleImage = cv2.cvtColor(flip, cv2.COLOR_BGR2GRAY)

    #blur the image to make granulation appear lighter in contrast to sunspots
    blur = cv2.GaussianBlur(grayscaleImage,(7,7),0)

    # apply thresholding
    threshValue, binaryImage = cv2.threshold(blur, threshold,255,cv2.THRESH_BINARY_INV)
    plt.imshow(binaryImage, cmap ='gray' )
    plt.gca().invert_yaxis()

    mask = np.zeros_like(binaryImage) # set up empty mask

    contours,_ = cv2.findContours(binaryImage, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    contours    = [c for c in contours if cv2.contourArea(c) > 100]
    biggest_contour = max(contours, key = cv2.contourArea)
    cv2.drawContours(mask, [biggest_contour], -1, 255, thickness=cv2.FILLED)

    # create final penumbra mask by subtracting the dialed umbra from the whole spot
    penumbra = mask - umbra_dilated

    return penumbra;

def check_masks(image, filename, penumbra, umbra, wantCheck, wantSave):
   '''
    depening on what was set in the beginning, this will display the individual images with the masks drawn in and save them, if that is desired
   '''
   if wantCheck==True and wantSave==True:
        plt.imshow(image, origin='lower', cmap='gray')
        plt.imshow(penumbra, origin='lower', alpha=0.5)
        plt.imshow(umbra, origin='lower', alpha=0.3)
        plt.savefig(mask_path+'/'+str(filename)+'_masked.png')
        plt.show()
        plt.close()
   elif wantCheck==True and wantSave==False:
        plt.imshow(image, origin='lower', cmap='gray')
        plt.imshow(penumbra, origin='lower', alpha=0.5)
        plt.imshow(umbra, origin='lower', alpha=0.3)
        plt.show()
        plt.close()
   elif wantCheck==False and wantSave==True:
        plt.imshow(image, origin='lower', cmap='gray')
        plt.imshow(penumbra, origin='lower', alpha=0.5)
        plt.imshow(umbra, origin='lower', alpha=0.3)
        plt.savefig(mask_path+'/'+str(filename)+'_masked.png')
        plt.close()
   return

def LimbDarkeningAR(name, location, binaryImage_P, binaryImage_U, binaryImage_F, r_sun):

    #There is a problem with opening png files like this, as they are rescaled to new values between 0 and 255.
    #This means that we lost the intensity calibration of the SDO files. I suggest that you save the image as a png,
    #and also as a numpy array to avoid this problem. Then just remove the line below and replace name with fov.
    #I'll just continue like everything is good.
    fov = name

    x_size, y_size = np.shape(fov) #If you ever have not square arrays, check if x and y are not flipped.

    #make a grid that contains the x coordinate of each pixel (Its basically a bunch of columns stating 0,1,2,3...400)
    pixel_grid_x = np.tile(np.arange(x_size), y_size).reshape([x_size,y_size])
    #Do the same but for y, so we transpose the array.
    pixel_grid_y = pixel_grid_x.T

    #convert to arcsec (is it just 0.5 or some other number?)
    arcsec_grid_x = pixel_grid_x * 0.5 + location[0]
    arcsec_grid_y = pixel_grid_y * 0.5 + location[1]

    #Now we ravel all the arrays into 1d arrays
    fov_ravel = np.ravel(fov)

    maskF = np.ravel(binaryImage_F) # fullspot
    maskU = np.ravel(binaryImage_U) # umbra
    maskP = np.ravel(binaryImage_P) # penumbra

    arcsec_x_ravel = np.ravel(arcsec_grid_x)
    arcsec_y_ravel = np.ravel(arcsec_grid_y)

    #We then make a new coordinate r, which is incorporates x and y.
    arcsec_r_ravel = np.sqrt(arcsec_x_ravel**2 + arcsec_y_ravel**2)

    #We then convert r to mu
    mu_ravel = np.sqrt(1 - (arcsec_r_ravel/r_sun)**2)

    #Now it is simply a matter of taking out all the values that we want
    #In our case we want to know, for each pixel inside the U and P what the mu, intensity (and fullspot)
    mu_F = mu_ravel[np.where(maskF)]
    int_F = fov_ravel[np.where(maskF)] # Full spot

    mu_U = mu_ravel[np.where(maskU)]
    int_U = fov_ravel[np.where(maskU)] # umbra

    mu_P = mu_ravel[np.where(maskP)]
    int_P = fov_ravel[np.where(maskP)] # penumbra

    #Now we output 3 arrays, one for U and one for P.
    F_arr = np.array([mu_F,int_F, int_F]).T
    U_arr = np.array([mu_U,int_U, int_U]).T
    P_arr = np.array([mu_P,int_P, int_P]).T

    return F_arr, U_arr, P_arr

def append_array_to_text_file(filename, numpy_array):
    """
    Appends information from a NumPy array to a text file.

    Parameters:
    - filename: The name of the text file to append to.
    - numpy_array: The NumPy array containing the data to be appended.

    Returns:
    - True if the operation was successful, False otherwise.
    """
    try:
        # Open the file in append mode
        with open(filename, 'a') as file:
            # Convert the NumPy array to a string and write it to the file
            np.savetxt(file, numpy_array, delimiter=',', fmt='%f')

        return True  # Operation was successful
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return False  # Operation failed

def polynomial_fit(mu, int0, n):
    # Define the polynomial function
    def poly_n(x, *coefficients):
        return sum(coefficients[i] * x**i for i in range(n + 1))

    # Fit the polynomial to the data
    popt, pcov = curve_fit(poly_n, mu, int0, p0=[0.0] * (n + 1))  # Initialize coefficients with zeros

    return popt  # Returns the coefficients as a list

def quadratic_law(x,a,b):
    return 1-a*x-b*x**2

def quadratic_fit(mu, int0):
    # Fit the model to the data
    params, pcov = curve_fit(quadratic_law, (1-mu), int0)  # Initialize coefficients with zeros

    return params  # Returns the coefficients as a list

def non_linear_law(x,a,b,c,d):
    return 1-a*(1-x**(1/2))-b*(1-x)-c*(1-x**(3/2))-d*(1-x**2)

def non_linear_fit(mu, int0):
    # Fit the model to the data
    params, pcov = curve_fit(non_linear_law, (mu), int0)  # Initialize coefficients with zeros

    return params  # Returns the coefficients as a list

#######################################################################################
#       PLEASE SET BEFORE RUNNING
#######################################################################################
# set email for conformation from the SDO Database
myemail = '----'

# set directionary where the output should appear/ where you run the script
mydir = '----'

# NOAA of the to-be-observed active region
NOAA = ---
# Latitude [deg]
latitude = ---
# put either the wavelength (e.g. Fe_1673A, or Halpha, ... ) or the oberservable  (e.g. 'linewidth','magnetogram',...)
wavelength = '---'
# the reference SDO HMI series with limb-darkneing removed (needed to create the masks)
No_LD_Series = '--' # DO NOT CHANGE!
# the data series in which CLV shall be studied:
LimbDark_Series = '---'  # set to the series you want to download/use


################## Time frame: ##############
###### OPTION ONE: ##########
# Specify your start date/time
start_date = Time('2014-05-15T21:00:00')
#it is assumed that the whole rotation will be observed, so the number of days is set to be 11
days = 11
# Spefify the cadance between two pictures [hours]
dt = 6

###### OPTION TWO: ##########
# if you dont want tp specify the cadance and rather use a list (array) of time stamps (to have more variety closer to the limb or to work with telescopes that only get data once per day) and input it here:
time_stamps = ['2014-11-23T21:00:00', '2014-11-23T22:00:00', '2014-11-23T23:00:00']  # Example time stamps
#(set time_stamps to be empty, if not used)

use_time_stamps = False #set to TRUE if you want to use the array above for download

############################################

# using the solar monitor, estimate the location of the first appearance of the active region [arcsec]
x_start_position = 750  #-936
y_start_position = 300   #260

# Spefify the desired size of the cropped 5out frame (square) [pixel]
frame_size = 700

# Specifiy threshold value (umbra/penumbra)
penumbra_threshold = 175
umbra_threshold = 83

plage_threshold = 225 

# if data was already downloaded, please set to FALSE
download_needed = False
#if you want to see the individual masked spots please set to 'TRUE'
wantCheck = True
#if you want to save those images to the output folder, please set to 'True':
wantSave = True


############################################################################################
# set up new directionary for output for this spot
#main_dir = 'NOAA_' + str(NOAA) + '_dt_' + str(dt) + 'h'
main_dir = 'NOAA_' + str(NOAA) + '-' + wavelength + '_dt_' + str(dt) + 'h'
new_path = check_path(mydir, main_dir)
#inside the new directionary, set up sub-dirs for images and different methods
FITS_path    = check_path(new_path, "FITS_files")
LD_FITS_path = check_path(new_path, "FITS_files_LimbDarkening")
png_path     = check_path(new_path, "png_images")
LD_results_path = check_path(new_path, "LimbDarkening_results")
if wantSave==True:
    mask_path = check_path(new_path, "masks")

#check if output file with that name already exists in that folder, delete it and make new one:
check_file(LD_results_path+'/lbd_P.txt')
check_file(LD_results_path+'/lbd_U.txt')
check_file(LD_results_path+'/lbd_F.txt')

##########################################################################################
# follow sunspot, get cropped frames for both series
##########################################################################################


# if data is not yet in the correct folder (True) this will initiate the communication with JSOC and download the FITS files
if download_needed == True:
    # make a request to the SDO datloessnitz@uni-potsdam.deabase and download the selected data
    file_name_array = downloadSDO(No_LD_Series, use_time_stamps, time_stamps, start_date, dt, 11, FITS_path)
    # download the same images without the removed limbdarkening
    file_name_array_LimbDark = downloadSDO(LimbDark_Series, use_time_stamps, time_stamps, start_date, dt, 11, LD_FITS_path)


#name array for data without limb darkening
file_name_array = np.loadtxt(new_path+'/'+'NOAA_'+str(NOAA)+"file_names_list.txt", dtype=str, comments="#", delimiter=",", unpack=False)
#name array for data with limb darkening
file_name_array_LimbDark = np.loadtxt(new_path+'/'+'NOAA_'+str(NOAA)+'_LimbDark_'+ "file_names_list.txt", dtype=str, comments="#", delimiter=",", unpack=False)
length = len(file_name_array)

# calculate time steps between images in series:
durations = get_time_steps(length)

# obtain information stored in the header of the first .FITS file for reference
Center_pix_x, Center_pix_y, delta_x, delta_y, R_sun, R_sun_pix, L_Carrington, B_Carrington, time_stamp  = get_header_info(0)

# get the location (x-position in pixel) of the frame for each image
if download_needed == True:
    frame_location = get_frame_locations(file_name_array, x_start_position, y_start_position, durations)
else:
    frame_location = get_frame_locations(file_name_array, x_start_position, y_start_position, durations)

#crop images, save them as pngs and as a fits cube
frame_location, y_location, cropped_cube = track_region(No_LD_Series, frame_location, frame_size, file_name_array, y_start_position, Center_pix_y, delta_y)
#crop images for series with limbdarkening
frame_location, y_location, cropped_cube_limbdark = track_region(LimbDark_Series, frame_location, frame_size, file_name_array, y_start_position, Center_pix_y, delta_y)


######################################################################################
# create masks and get intensity values for fullspot, penumbra and umbra
#####################################################################################

print(f'creating masks...')

for i in range(length):
    filename = file_name_array[i]
    number = str(i)

    #print(f'creating masks... current image: {i}')

    # get header information from the fits files:
    Center_pix_x, Center_pix_y, delta_x, delta_y, R_sun, R_sun_pix, L_Carrington, B_Carrington, time_stamp = get_header_info(i)

    # create the masks for the sunspot:
    fullSpot = mask_sunspot(filename, penumbra_threshold)                # the full sunspot
    umbra, umbra_dilated = mask_umbra(filename, umbra_threshold)         # umbra only
    penumbra = mask_penumbra(filename, umbra_dilated, penumbra_threshold)# penumbra only

    #translate to determine coordinates of each frame in lower left corner
    origin_X = frame_location[i] + frame_size/2
    origin_X_arcsec = (-1) * (origin_X - Center_pix_x) * delta_x
    origin_Y = y_location + frame_size/2
    origin_Y_arcsec = (-1) * (origin_Y - Center_pix_y) * delta_y
    location = np.array([int(origin_X_arcsec),int(origin_Y_arcsec)])

    # now choose the coresponding image from the series with limbdarkening included
    image = cropped_cube_limbdark[i]
    # apply masks on those and retrieve intensity values
    F_LD_arr, U_LD_arr, P_LD_arr = LimbDarkeningAR(image, location, penumbra, umbra, fullSpot, R_sun)

    # save to file:
    append_array_to_text_file(LD_results_path+'/'+'lbd_F.txt', F_LD_arr)
    append_array_to_text_file(LD_results_path+'/'+'lbd_U.txt', U_LD_arr)
    append_array_to_text_file(LD_results_path+'/'+'lbd_P.txt', P_LD_arr)

    medianF = np.median(F_LD_arr[:,1])
    medianP = np.median(P_LD_arr[:,1])
    medianU = np.median(U_LD_arr[:,1])

    check_masks(image, filename, penumbra, umbra, wantCheck, wantSave)

plt.close()

#####################################################################################
# LIMB DARKENING LAW BY NECKEL AND LABS
lbd = np.array([0.31442323,  1.35348459, -1.8265809 ,  2.35713264, -1.68753739, 0.48907787])
#####################################################################################

print(f'calculate CLV...')

# load the penumbra values and create scatter plot
pen = np.loadtxt(LD_results_path+'/'+'lbd_P.txt', delimiter=',')
plt.scatter(pen[:,0], pen[:,1], alpha=0.2)
plt.title('penumbra of '+str(NOAA)+' ('+wavelength+')')
plt.savefig(LD_results_path+'/'+'prenorm_'+'penumbra.png', dpi=300)
plt.show()

# load the numbra values and create scatter plot
umb = np.loadtxt(LD_results_path+'/'+'lbd_U.txt', delimiter=',')
plt.scatter(umb[:,0], umb[:,1], alpha=0.2)
plt.title('umbra of '+str(NOAA)+' ('+wavelength+')')
plt.savefig(LD_results_path+'/'+'prenorm_'+'umbra.png', dpi=300)
plt.show()

# load the fullspot values and create scatter plot
ful = np.loadtxt(LD_results_path+'/'+'lbd_F.txt', delimiter=',')
plt.scatter(ful[:,0], ful[:,1], alpha=0.2)
plt.title('full sunspot '+str(NOAA)+' ('+wavelength+')')
plt.savefig(LD_results_path+'/'+'prenorm_'+'fullspot.png', dpi=300)
plt.show()

# arange the mu and intensity values
mu_p    = np.flip(pen[:,0][np.argsort(pen[:,0] )]) #mu of penumbra pixels
int_p   = np.flip(pen[:,1][np.argsort(pen[:,0] )]) #intensities in penumbra
int_p_s = np.flip(pen[:,2][np.argsort(pen[:,0] )]) #penumbra with quiet sun

mu_u    = np.flip(umb[:,0][np.argsort(umb[:,0] )]) #mu of umbra
int_u   = np.flip(umb[:,1][np.argsort(umb[:,0] )]) #intensities of umbra
int_u_s = np.flip(umb[:,2][np.argsort(umb[:,0] )]) #umbra with quiet sun

mu_f    = np.flip(ful[:,0][np.argsort(ful[:,0] )]) #mu of fullspot
int_f   = np.flip(ful[:,1][np.argsort(ful[:,0] )]) #intensities of fullspot
int_f_s = np.flip(ful[:,2][np.argsort(ful[:,0] )]) #fullspot with quiet sun


# poly-fit of degree n=3
val_p = polynomial_fit(mu_p, int_p, 3) # penumbra
val_u = polynomial_fit(mu_u, int_u, 3) # umbra
val_f = polynomial_fit(mu_f, int_f, 3) # fullspot

# continuos mu for the fits:
mu_fits = np.arange(1000)/1000

#make inverse mu for each set:
mu_p_inv    =1-mu_p     # penumbra mu
mu_u_inv    =1-mu_u     # umbra mu
mu_f_inv    =1-mu_f     # fullspot mu
mu_fits_inv =1-mu_fits  # cont. mu for plot


#now fit poly law to that inverse data, since the firt parameter will be the y-value at x=0
val_p_inv = polynomial_fit(mu_p_inv, int_p, 3)
val_u_inv = polynomial_fit(mu_u_inv, int_u, 3)
val_f_inv = polynomial_fit(mu_f_inv, int_f, 3)

#the first parameter can therefore be used to normalize the curve for this inverse case
val_p_inv_norm = val_p_inv/val_p_inv[0]
val_u_inv_norm = val_u_inv/val_u_inv[0]
val_f_inv_norm = val_f_inv/val_f_inv[0]

#now we repopulate the fitted curve with values that can then be used to go back to the real mus
val_p_norm = polynomial_fit(mu_fits, np.polyval(val_p_inv_norm[::-1], mu_fits_inv), 3)
val_u_norm = polynomial_fit(mu_fits, np.polyval(val_u_inv_norm[::-1], mu_fits_inv), 3)
val_f_norm = polynomial_fit(mu_fits, np.polyval(val_f_inv_norm[::-1], mu_fits_inv), 3)

#now we also normalize the data to this
int_p_norm = int_p/val_p_inv[0]
int_u_norm = int_u/val_u_inv[0]
int_f_norm = int_f/val_f_inv[0]

#now we can fit data to this too:
#get parameters for quadratic fit:
quad_p = quadratic_fit(mu_p, int_p_norm)
quad_u = quadratic_fit(mu_u, int_u_norm)
quad_f = quadratic_fit(mu_u, int_u_norm)
quad_qs = quadratic_fit(mu_fits, np.polyval(lbd[::-1], mu_fits)) # quiet sun

#get parameters for non-linear fit:
non_linear_p = non_linear_fit(mu_p, int_p_norm)
non_linear_u = non_linear_fit(mu_u, int_u_norm)
non_linear_f = non_linear_fit(mu_f, int_f_norm)


#####################################################################################################
# PLOT THE INTENSITY CURVES
####################################################################################################
print(f'preparing CLV plots...')

##### FULL SPOT: ########
# if you want a plot of the intensities for the entire sunspot:
# set up figure:
fig, ax = plt.subplots(figsize=(8, 4))

# plot intensity values of the whole spot (normalized)
plt.scatter(mu_f,int_f/val_f_inv[0], color='violet', s=0.1, alpha=0.1)

plt.plot(mu_fits, np.polyval(lbd[::-1], mu_fits), color='black', linestyle='dotted', linewidth=2, label='Quiet Sun')
#plt.plot(mu_fits, quadratic_law(1-mu_fits, quad_qs[0], quad_qs[1]), color='gray', linestyle='dashdot', linewidth=2, label='Quadratic Quiet Sun')
plt.plot(mu_fits, np.polyval(val_f_norm[::-1], mu_fits), color='indigo', linewidth=2, label='Polynomial Fit Sunspot')
plt.plot(mu_fits, quadratic_law(1-mu_fits, quad_f[0], quad_f[1]), color='m', linestyle='dashed', linewidth=2, label='Quadratic Fit Sunspot ' )

plt.legend()
plt.title('CLV of full sunspot NOAA '+str(NOAA)+' ('+wavelength+')')
plt.xlabel(r'$\mu$')
plt.ylabel('Int (normalized)')
plt.savefig(main_dir+'/'+'CLV_'+str(NOAA)+'_'+wavelength+'_fullspot.pdf', dpi=300)
plt.show()
plt.close()

###### UMBRA AND PENUMBRA: #############
# if you want umbra and penumbra to be fitted seperatly
# set up figure
fig, ax = plt.subplots(figsize=(8, 4))

# plot intensity values
plt.scatter(mu_u,int_u/val_u_inv[0], s=0.1, alpha=0.2)
plt.scatter(mu_p,int_p/val_p_inv[0], color='indianred', s=0.1, alpha=0.1)

# plot polynomial fits to penumbra, umbra and quiet sun
plt.plot(mu_fits, np.polyval(lbd[::-1], mu_fits), color='black', linestyle='dotted', linewidth=2, label='Quiet Sun')
#plt.plot(mu_fits, quadratic_law(1-mu_fits, quad_qs[0], quad_qs[1]), color='gray', linestyle='dotted', linewidth=2, label='Quadratic Quiet Sun')
plt.plot(mu_fits, np.polyval(val_u_norm[::-1], mu_fits), color='blue', linestyle='solid', linewidth=2, label='Polynomial Fit Umbra')
plt.plot(mu_fits, quadratic_law(1-mu_fits, quad_u[0], quad_u[1]), color='dodgerblue', linestyle='dashed', linewidth=2, label='Quadratic Fit Umbra' )
plt.plot(mu_fits, np.polyval(val_p_norm[::-1], mu_fits), color='darkred', linestyle='solid', linewidth=2, label='Polynomial Fit Penumbra')
plt.plot(mu_fits, quadratic_law(1-mu_fits, quad_p[0], quad_p[1]), color='red', linestyle='dashed', linewidth=2, label='Quadratic Fit Penumbra' )
#plt.plot(mu_fits, non_linear_law(mu_fits, non_linear_p[0], non_linear_p[1], non_linear_u[2], non_linear_u[3]), label='Non-Linear Fit Penumbra' )
#plt.plot(mu_fits, non_linear_law(mu-fits, non_linear_u[0], non_linear_u[1], non_linear_u[2], non_linear_u[3]), label='Non_Linear Fit Umbra' )

plt.legend()
plt.title('CLV of umbra & penumbra NOAA '+str(NOAA)+' ('+wavelength+')')
plt.xlabel(r'$\mu$')
plt.ylabel('Int (normalized)')
plt.ylim([0.25, 2])
plt.savefig(main_dir+'/'+'CLV_'+str(NOAA)+'_'+wavelength+'.pdf', dpi=300)
plt.show()
plt.close()

######## SAVE RESULTS AND FIT PARAMETERS: ##########

# write the resulting fit-parameters to a txt in the output folder:
results = open(LD_results_path+'/NOAA'+str(NOAA)+'_'+wavelength+'_'+"fit_parameter.txt", "w+")
results.write('NOAA: {0}\nwavelength/data-set: {1}\nUMBRA:\npolynomial-fit:  {2}\nquadratic-fit: {3}\n\nPENUMBRA:\npolynomial-fit: {4}\nquadratic-fit: {5}\n\nFULL-SPOT:\npolynomial-fit: {6}\nquadratic-fit: {7}'.format(NOAA, wavelength, val_u_norm, quad_u, val_p_norm, quad_p, val_f_norm, quad_f,))
results.close()
