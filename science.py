#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Filename: science.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)
from astropy.io import fits
import numpy
import matplotlib.pyplot as plt
from astropy.visualization import ImageNormalize, LinearStretch, ZScaleInterval
from astropy.stats import sigma_clip
from scipy.stats import mode
from astroscrappy import detect_cosmics


def reduce_science_frame(
    science_filename,
    median_bias_filename,
    median_flat_filename,
    median_dark_filename,
    reduced_science_filename="reduced_science.fits",
):
   
    # Gets the science fits data from science_filename
    science_frame = (fits.getdata(science_filename)).astype('float32')

    # Gets the exposure time of the science image
    sci_exp = (fits.getheader(science_filename))['EXPTIME']

    # Reduces the science image by subtracting the provided median bias and dark file
    reduced_science = science_frame - fits.getdata(median_bias_filename) - (fits.getdata(median_dark_filename)*sci_exp)

    # Divides out the normalizes median_flat image
    reduced_science /= fits.getdata(median_flat_filename)

    # Performs cosmic ray detection
    mask, output = detect_cosmics(reduced_science)

    # Sets up fits write to file
    header = fits.getheader(science_filename)
    primary = fits.PrimaryHDU(data=output.data, header=header)
    hdul = fits.HDUList([primary])
    hdul.writeto(reduced_science_filename, overwrite=True)

    return output
