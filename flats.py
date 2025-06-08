#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Filename: flats.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)
from astropy.io import fits
import numpy
import matplotlib.pyplot as plt
from astropy.visualization import ImageNormalize, LinearStretch, ZScaleInterval
from astropy.stats import sigma_clip
from astroscrappy import detect_cosmics

def create_median_flat(
    flat_list,
    bias_filename,
    median_flat_filename,
    dark_filename,
):
    """This function must:

    - Accept a list of flat file paths to combine as flat_list. Make sure all
      the flats are for the same filter.
    - Accept a median bias frame filename as bias_filename (the one you created using
      create_median_bias).
    - Read all the images in flat_list and create a list of 2D numpy arrays.
    - Read the bias frame.
    - Subtract the bias frame from each flat image.
    - Optionally you can pass a dark frame filename as dark_filename and subtract
      the dark frame from each flat image (remember to scale the dark frame by the
      exposure time of the flat frame).
    - Use a sigma clipping algorithm to combine all the bias-corrected flat frames
      using the median and removing outliers outside 3-sigma for each pixel.
    - Create a normalised flat divided by the median flat value.
    - Save the resulting median flat frame to a FITS file with the name
      median_flat_filename.
    - Return the normalised median flat frame as a 2D numpy array.

    """

    # Gathers the exposure time from the header of each image in flat_list
    exptime = []
    for file in flat_list:
        header = fits.getheader(file)
        exptime.append(header['EXPTIME'])

    # Reads the data from the file paths in flat_list
    flat_images = []
    for i in range(0, numpy.size(flat_list)):
        data = fits.getdata(flat_list[i])
        flat_images.append((data).astype('float32'))

    # Reads the provided bias_image
    bias = fits.getdata(bias_filename)

    # Reads median dark image if one is provided
    if dark_filename:
        dark = fits.getdata(dark_filename)
    else:
        dark = None

    # Corrects flat image by subtracting the median bias and dark images
    flat_corrected = []
    for i in range(0, len(flat_images)):
        if dark is not None:
            flat_corrected.append(flat_images[i] - bias -(dark * exptime[i]))
        else:
            flat_corrected.append(flat_images[i] - bias)

    # Performs sigma_clipping on flat image
    flat_images_masked = sigma_clip(flat_corrected, cenfunc='median', sigma=3, axis=0)

    # Creates flattened median image from sigma clipped arrays
    collapsed_flat = numpy.ma.median(flat_images_masked, axis=0)

    # Normalizes the collapsed_flat image by dividing by its median
    ouput = collapsed_flat / numpy.ma.median(collapsed_flat)

    # Replaces values of 0 in the image with 1 to avoid division errors in science reduction
    ouput[collapsed_flat == 0] = 1

    # Gets a standard flat header
    header = fits.getheader(flat_list[0])

    # Sets up fits write to and updates header with relavent information
    primary = fits.PrimaryHDU(data=ouput.data, header=header)
    primary.header['COMMENT'] = 'Normalized flat field image'
    primary.header['BIASFILE'] = 'bias_sc.fits'
    primary.header['DARKFILE'] = 'dark_sc.fits'
    hdul = fits.HDUList([primary])
    hdul.writeto(median_flat_filename, overwrite=True)

    return ouput


def plot_flat(
    median_flat_filename,
    ouput_filename="median_flat.png",
    profile_ouput_filename="median_flat_profile.png",
):
    """This function must:

    - Accept a normalised flat file path as median_flat_filename.
    - Read the flat file.
    - Plot the flat frame using matplotlib.imshow with reasonable vmin and vmax
      limits. Save the plot to the file specified by output_filename.
    - Take the median of the flat frame along the y-axis. You'll end up with a
      1D array.
    - Plot the 1D array using matplotlib.
    - Save the plot to the file specified by profile_output_filename.

    """
    # Gets median flat image produced by create_median_flat
    median_flat = fits.getdata(median_flat_filename)

    # Defines a normalization to apply to the median_flat image
    norm = ImageNormalize(median_flat, interval=ZScaleInterval(), stretch=LinearStretch())
    
    # Plots the flat image as a 2D array
    plt.clf()
    plt.imshow(median_flat, origin='lower', norm=norm, cmap='YlOrBr_r')
    plt.savefig(ouput_filename)

    # Clears plot
    plt.clf()

    # Plots the flat image along its y-axis
    plt.plot(numpy.median(median_flat, axis = 0))
    plt.savefig(profile_ouput_filename)

    return
