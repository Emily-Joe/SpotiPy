import os
import numpy as np
import astropy.io.fits as f
from astropy.time import Time
import astropy.units as u
import datetime
from sunpy.net import Fido, attrs as a

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

def downloadSDO(series, start_date, dt, days, download_path):
    """
    make a search query to the SDO Database to extract images of the sun in the desired timeframe and from the desired series.Checking the database might take some minutes.
    If the search-request is successful, an email will be send to your email account and the download will begin automatically. The images will be downloaded as .FITS files into the directory stated in download_path. Additionally an array with all the file-names will be available.

    Parameters
    ----------
    series: str
            the name of the data series from which data shall be downloaded, the common series for this work is 'hmi.Ic_720s'
    start_date:  Time
            stets the start date and time in the format YYYY-MM-DDT00:00:00

    dt  :   int
            specifies the cadance between two images in hours

    days  : int
            specifies the number of days over wich the observation shall take place

    download_path : str
            sets path and directionary for .FITS-files to be downloaded to

    Returns
    -------
    file_name_array
            an array containing the names of all the downloaded files in chronological order

    :Authors:
        Emily Joe Loessnitz (2024)
    """

    #calcultate how many images should be downloaded to meet requiremtns set by the search query
    steps = int(days*24 / dt)

    print(f'Preparing to download {steps} images starting from {start_date}, spanning {days} days, with a {dt} hour cadence.')
    print('\n\n\nChecking Database...')

    # Calculate start and end times for the search
    start_time = start_date #+ datetime.timedelta(hours=i)
    end_time = start_time + datetime.timedelta(hours=dt*steps) #only select the right images

    # Make the search query
    res = Fido.search(a.Time(start_time.iso, end_time.iso), a.jsoc.Series(series), a.Sample(dt*u.hour) , a.jsoc.Notify(myemail))

    # give an overview over available data for requested timeframe. If everything looks correct, enter 'y', if not enter 'n' and try again. If the key 'y' is pressed, the download begins
    print(res.show('TELESCOP', 'INSTRUME', 'T_OBS'))
    print('Does this look correct?')
    resp = input('Press [y] to download and [n] to exit\n')
    if resp == 'y':
        print('Downloading ...')
        downloaded_files = Fido.fetch(res, path=download_path+'/{file}')
    else:
        print('END')
        pass
    print('END')


    # to get all the file-names for later use, read out the download-directionary and set up empty array 'file_names' to save all the file names to
    dir_path = download_path
    #to make sure the entries are in chronological order, the directionary gets sorted first
    dirs = sorted(os.listdir(dir_path))
    # empty array to store filenames
    file_names = []

    #depending on the series choosen the correct name for the namelist.txt is choosen
    if series=='hmi.Ic_720s':
        file_names_list = open(new_path + '/' +'NOAA_' + str(NOAA) + '_LimbDark_' +  "file_names_list.txt", "w+")
    else:
         file_names_list = open(new_path + '/' +'NOAA_' + str(NOAA) + "file_names_list.txt", "w+")

    #Iterate directory, add names from download-folder to array and to the .txt-list
    for file_path in dirs:
        name = str(file_path)
        file_names_list.write(name + '\n')
        file_names.append(name)
    file_names_list.close() #close the txt
    file_name_array = np.asarray(file_names) #make sure the name-array is actually an array

    print('Download complete!')

    return FITS_path, file_name_array ;

##########################################################################################
#       PLEASE SET BEFORE RUNNING
##########################################################################################
# set email for conformation from the SDO Database
myemail = 'loessnitz@uni-potsdam.de'

# set directionary where the output should appear
mydir = '/work2/loessnitz/'

# NOAA of the to be observed active region
NOAA = 12218
# Latitude [deg]
latitude = 16

# Specify your start date/time
start_date = Time('2014-11-23T21:00:00')
time_offset = 0 #[hours]
#it is assumed that the whole rotation will be observed, so the number of days is set to be 11
days = 1

# Spefify the cadance between two pictures [hours]
dt = 1

# if data was already downloaded, please set to FALSE
download_needed = True

# select if limbdarkening in spot should be also studied, this would include an additional download of a data set. If limbdarkening is of imporatance set to 'True'
With_LimbDarkening = True

# the names of the data from JSOC, change if there is a change in naming conventions or a different series is wanted:
No_LD_Series = 'hmi.Ic_noLimbDark_720s' #normal series
LD_Series = 'hmi.Ic_720s'         #for limb-darkening studies

###########################################################################################
###########################################################################################
# set up new directionary for output for this spot
main_dir = 'NOAA_' + str(NOAA) + '_dt_' + str(dt) + 'h'
new_path = check_path(mydir, main_dir)

#inside the new directionary, set up sub-dirs for images and different methods
FITS_path    = check_path(new_path, "FITS_files")
png_path     = check_path(new_path, "png_images")
centroid_path= check_path(new_path, "central_moments")
ellipse_path = check_path(new_path, "ellipse_center")

#if LimbDarkening is wanted, this creates the infrastructure for that:
if With_LimbDarkening:
    LD_FITS_path = check_path(new_path, "FITS_files_LimbDarkening")
    LD_results_path = check_path(new_path, "LimbDarkening_results")

# if data is not yet in the correct folder (True) this will initiate the communication with JSOC and download the FITS files
if download_needed == True:
    if not With_LimbDarkening:
        # make a request to the SDO database and download the selected data
        file_name_array = downloadSDO(No_LD_Series, (start_date+time_offset*u.hour), dt, days, FITS_path)
    elif With_LimbDarkening:
        # make a request to the SDO database and download the selected data
        file_name_array = downloadSDO(No_LD_Series, start_date, dt, days, FITS_path)
        # download the same images without the removed limbdarkening
        file_name_array_LimbDark = downloadSDO(LD_Series, start_date, dt, days, LD_FITS_path)



