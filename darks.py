#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Filename: darks.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)
from astropy.io import fits
import numpy
import matplotlib.pyplot as plt
from astropy.visualization import ImageNormalize, LinearStretch, ZScaleInterval
from astropy.stats import sigma_clip

def create_median_dark(dark_list, bias_filename, median_dark_filename):
    """This function must:

    - Accept a list of dark file paths to combine as dark_list.
    - Accept a median bias frame filename as bias_filename (the one you created using
      create_median_bias).
    - Read all the images in dark_list and create a list of 2D numpy arrays.
    - Read the bias frame.
    - Subtract the bias frame from each dark image.
    - Divide each dark image by its exposure time so that you get the dark current
      per second. The exposure time can be found in the header of the FITS file.
    - Use a sigma clipping algorithm to combine all the bias-corrected dark frames
      using the median and removing outliers outside 3-sigma for each pixel.
    - Save the resulting dark frame to a FITS file with the name median_dark_filename.
    - Return the median dark frame as a 2D numpy array.

    """

    # Reads the dark image filepaths from dark_list
    dark_images = []
    for i in range(0, numpy.size(dark_list)):
        data = fits.getdata(f'{dark_list[i]}')
        dark_images.append((data).astype('float32'))

    # Subtracts the median bias frame from each dark image
    dark_corrected = []
    for i in range(0, len(dark_images)):
        dark_corrected.append(dark_images[i] - fits.getdata(bias_filename))

    # Calculates the average amount of dark current noise per second of exposure time
    dark_per_sec = []
    for i in range(0, len(dark_corrected)):
        dark_per_sec.append(dark_corrected[i] / 90)

    # Performs sigma clipping on dark_per_sec images
    dark_images_masked = sigma_clip(dark_per_sec, cenfunc='median', sigma=3, axis=0)

    # Returns median dark per second image
    output = numpy.ma.median(dark_images_masked, axis=0)


    # Sets up the fits file with median_dark_filename, updates header information
    header = fits.getheader(dark_list[0])
    primary = fits.PrimaryHDU(data=output.data, header=header)
    primary.header['EXPTIME'] = 1
    primary.header['EXPOSURE'] = 1
    hdul = fits.HDUList([primary])
    hdul.writeto(median_dark_filename, overwrite=True)
    

    return output
