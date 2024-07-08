# Packages used
import numpy as np
from astropy.io import fits
import astropy.wcs as wcs
from astropy.coordinates import SkyCoord
import astropy.units as u
import glob
import os
from copy import deepcopy
import matplotlib.pyplot as plt
from matplotlib import colors
import pyregion
import cupy as cp
from astropy.nddata.utils import Cutout2D

# Paths of the images
def images_paths(field_name, directory_path = "work/red"): 
    # directory_path: in wich folder we have to find the images, default in the work/red folder
    # field_name: name of the sky field of interest
    
    contents = os.listdir(directory_path)
    images_paths_QHY = []
    images_paths_iKon = []
    images_paths_all = []
    for i in contents:
        days_path = directory_path + "/" + i
        day_contents = os.listdir(days_path)

        for j in day_contents:        
            if field_name in j:
                if "QHY" in j:
                    images_paths_QHY.append(days_path + "/" + j)
                else: 
                    images_paths_iKon.append(days_path + "/" + j)
                images_paths_all.append(days_path + "/" + j)
    return images_paths_QHY, images_paths_iKon, images_paths_all

# Get some parameters
def getinfo(hdu):
    image = cp.asarray(hdu[0].data)
    header = hdu[0].header
    fwhm = hdu[0].header["FWHM"]
    filter = hdu[0].header["FILTER"]
    w_image = wcs.WCS(hdu[0].header)
    dateob = hdu[0].header["JD-OBS"]
    sky = hdu[0].header["FLUXSKY"]
    EXPT1 = hdu[0].header["EXPT1"]
    #readnoise = hdu[0].header['RDNOISE']
    return image, w_image, header, dateob, filter, fwhm, sky, EXPT1#, readnoise

# Open the images of interest
def get_images(image_paths):
    open_images = np.array([fits.open(p) for p in image_paths])
    headers = [o[0].header for o in open_images]
    images = np.array([o[0].data for o in open_images])
    Nimage=len(image_paths)
    return open_images, images, Nimage, headers

# Prameters from header
def header_param(headers):
    readnoises = np.array([headers[i]['RDNOISE'] for i in range(len(headers))])
    fwhms = np.array([headers[i]["FWHM"] for i in range(len(headers))])
    w_images = np.array([wcs.WCS(headers[i]) for i in range(len(headers))])
    return readnoises, fwhms, w_images

# Transform the RA DEC coordinates to XY coordinates
def radec_to_xy(target_radec, w_image):
    c = SkyCoord(target_radec[0], target_radec[1], frame='icrs', unit="deg")
    x, y = w_image.all_world2pix(c.ra,c.dec, 1, quiet=True)
    return x, y

# Cut the images
def crop_images(image, x, y, size, w, header):
    cut = Cutout2D(image, (x,y), (size,size), wcs=w)
    cut_image = cut.data
    cut_wcs = cut.wcs
    header.update(cut_wcs.to_header(relax=True))
    return cut_image

# Plot of the images cropped
def plot_crop_images(Nimage, crops, xgls, ygls):
    figure, axis = plt.subplots(Nimage, len(crops), figsize=(10, 10))
    for i in range(Nimage):
        for j in range(len(crops)):
            axis[i][j].imshow(crops[j][i], cmap='gray', vmin=200, vmax=1000, origin= 'lower')
        axis[i][0].set_title("Ephoc "+str(i+1))
    plt.tight_layout()
    
    # Find stars in the images using Pythonphot
def find_stars(Nimage, images, readnoises, fwhms):
    xstars = []
    ystars = []   
    for i in range(Nimage):   
        skymod, skysig, skyskw = pp.mmm.mmm(images[i],readnoise=readnoises[i])
        hmin = skysig*5
        print(hmin)
        xstar, ystar, flux, sharp, roundness = pp.find.find(images[i],hmin,fwhms[i])
        xstars.append(xstar)
        ystars.append(ystar)
    return xstars,ystars

# Plot the images pointing the stars
def plot_find_stars(Nimage, images, xstars, ystars, readnoises):
    figure, axis = plt.subplots(1, Nimage, figsize=(20, 20))
    for i in range(Nimage):
        axis[i].imshow(images[i], cmap='gray', vmin=200, vmax=1000, origin= 'lower')
        axis[i].plot(xstars[i], ystars[i],'o', ms=10, mfc='none', lw=2, mec='r')
        axis[i].set_title("Ephoc " + str(i+1))
    plt.tight_layout()
    
# Find the stars near the target
def stars_near_target(Nimage, image, xstars, ystars, xgls, ygls, sizeexcl):
    stars_position_ords = []
    stars_sky_ords = []
    bright_ords = []
    for i in range(Nimage):
        excl = int(np.shape(image[i])[1]/2)
        # Stars position excluding the target
        stars_position=np.array([(x,y) for x,y in zip(xstars[i],ystars[i])  if not ((excl-sizeexcl<x<excl+sizeexcl) and (excl-sizeexcl<y<excl+sizeexcl))])
        # Cutout of the stars excluding the target
        stars_sky=[image[i][int(np.round(y)-20):int(np.round(y)+20),int(np.round(x)-20):int(np.round(x)+20)] for x,y in zip(xstars[i],ystars[i]) if not ((excl-sizeexcl<x<excl+sizeexcl) and (excl-sizeexcl<y<excl+sizeexcl))]
        
        # Max bright of the cutouts to organize the stars in function its bright
        bright=np.array([stars_sky[i].max() if len(stars_sky[i])!=0 else 0 for i in range(len(stars_sky))])
        
        # Organization
        bright_ord=np.array(np.sort(bright))
        stars_sky_ord=[stars_sky[list(bright).index(i)] for i in np.sort(bright)]
        stars_position_ord=np.array([stars_position[list(bright).index(i)] for i in np.sort(bright)])
        print("There are", len(stars_position_ord), "stars near the target")
        stars_position_ords.append(stars_position_ord)
        stars_sky_ords.append(stars_sky_ord)
        bright_ords.append(bright_ord)
        
    return stars_position_ords, stars_sky_ords, bright_ords