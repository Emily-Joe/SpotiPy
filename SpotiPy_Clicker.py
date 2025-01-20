import numpy as np
import matplotlib.pyplot as plt
from scipy import stats, optimize
import astropy.io.fits as f
from astropy.time import Time, TimeDelta
import astropy.units as u
import datetime
import os
from sunpy.net import Fido, attrs as a

def download_preview(start_date, FITS_path):

    print(f'Preparing to download first images of the series at {start_date}')
    print('\n\n\nChecking Database...')

    # Calculate start and end times for the search
    start_time = start_date #+ datetime.timedelta(hours=i)
    end_time = start_time + 100*u.second #only select the right images

    # Make the search query
    res = Fido.search(a.Time(start_time.iso, end_time.iso), a.jsoc.Series('hmi.Ic_noLimbDark_720s'), a.Sample(700*u.second), a.jsoc.Notify(myemail))

    print(res.show('TELESCOP', 'INSTRUME', 'T_OBS'))
    print('Does this look correct?')
    resp = input('Press [y] to download and [n] to exit\n')
    if resp == 'y':
        print('Downloading ...')
        downloaded_files = Fido.fetch(res, path=FITS_path+'/{file}')
    else:
        print('END')
        pass
    print('END')

    dir_path = FITS_path
    dirs = sorted(os.listdir(dir_path))

    # list to store files
    file_names = []

    # Iterate directory=
    for file_path in dirs:
        name = str(file_path)
        file_names.append(name)

    file_name_array = np.asarray(file_names)
    filename = str(file_name_array[0])

    print('Download complete!')

    return filename;

def get_pixel_coordinates(data, name, path):

    #Waits for mouse click on image, saves coordinates, closes image.
    def onclick(event):
        ix, iy = event.xdata-midpoint_x, event.ydata-midpoint_y
        print('x = %d, y = %d' % (ix, iy))
        coords.append((ix, iy))
        if len(coords) == 1:
            fig.canvas.mpl_disconnect(cid)
            plt.close()

    #get midpoints and pixel-to-arcsec
    midpoint_x =  data[1].header['CRPIX1']
    midpoint_y =  data[1].header['CRPIX2']
    delta_x    =  data[1].header['CDELT1']
    delta_y    =  data[1].header['CDELT2']


    #get image
    image = data[1].data

    #Set up the imshow
    coords = []
    coords_arcsec = []

    dpi = 100 # Adjust dpi to control the size of the window
    height, width = image.shape
    figsize = 20, 20

    fig = plt.figure(figsize=figsize, dpi=dpi)
    ax = fig.add_subplot(111)
    ax.imshow(image, origin='lower', cmap='gray')
    cid = fig.canvas.mpl_connect('button_press_event', onclick)

    plt.show()

    # translate to arcsec
    x_start_coord = coords[0][0] * delta_x * (-1)
    y_start_coord = coords[0][1] * delta_y * (-1)
    coords_arcsec.append((x_start_coord, y_start_coord))


    #Now open the datafile and see if this calibration has been done already
    filename = path +'/'+ 'coordinate_checked.txt'
    found = False

    # Check if the file exists
    if not os.path.exists(filename):
    # Create the file if it does not exist
        with open(filename, 'w') as file:
            print(f"File '{filename}' not found. Created new file.")
    #check if this file has been checked before
    with open(filename, 'r') as file:
        for line in file:
            if name in line.strip():
                print('Coordinates for ' + name + ' already saved, skipping!')
                found = True
                break
    #If not, write the name plus coordiantes
    if not found:
        with open(filename, 'a') as file:
            file.write(name+','+str(coords[0][0]) + ',' + str(coords[0][1])+',' + str(coords_arcsec[0][0]) + ',' + str(coords_arcsec[0][1]) + '\n')
            print('Coordinates for ' + name + ' saved.')


    return coords

def extract_columns(filename):
    # gets the coordinates that were saved before
    data = []
    with open(filename, 'r') as file:
        for line in file:
            columns = line.strip().split(',')
            # Skip the 0th column and take the 1st and 2nd columns
            if len(columns) >= 3:
                data.append((float(columns[1]), float(columns[2]), float(columns[3]), float(columns[4])))
    return np.array(data)


##########################################################################
# set email for conformation from the SDO Database
myemail = '---'

# set directionary where the output should appear
mydir = '---'

# NOAA of the to be observed active region
NOAA = 13034
# Latitude [deg]
latitude = 1
# Specify the cadance between two pictures [hours]
dt = 6
# specify wavelentgh or data set:
wavelength = 'Fe'


# if you already downloaded the fits files for this spot, set to TRUE, if you want to download only the first image of the series to find the starting coordinates, set to FALSE
FITS_downloaded = False
# if FALSE, please enter the exact timestamp of the first sighting of the sunspot.
start_date =  Time('2022-06-12T23:00:00')

############################################################################


# Check whether the specified path exists or not
directionary = 'NOAA_' + str(NOAA) + '-' + wavelength + '_dt_' + str(dt) + 'h'
new_path = os.path.join(mydir, directionary)
existence = os.path.exists(new_path)
if not existence:
   os.makedirs(new_path)

SDO_data = "FITS_files"
FITS_path = os.path.join(new_path, SDO_data)
existence = os.path.exists(FITS_path)
if not existence:
   os.makedirs(FITS_path)


# if needed, download first image in series to calculate start positon:
if FITS_downloaded == False:
    dataname = download_preview(start_date, FITS_path)
else:
    file_name_array = np.loadtxt(new_path+'/'+'NOAA_'+str(NOAA)+"file_names_list.txt", dtype=str, comments="#", delimiter=",", unpack=False)
    dataname = file_name_array[0]

data_path = FITS_path+'/'+ dataname
data = f.open(data_path)

# click on spot in image to save coordinates in txt file, that can then be read out:
coords = get_pixel_coordinates(data, dataname, new_path)
start_coordinates = extract_columns(new_path+'/'+'coordinate_checked.txt')
