#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Filename: bias.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

from astropy.io import fits
import numpy
import matplotlib.pyplot as plt
from astropy.visualization import ImageNormalize, LinearStretch, ZScaleInterval
from astropy.stats import sigma_clip

def create_median_bias(bias_list, median_bias_filename):

    # Reads the fits filepaths provided by bias_list 
    bias_images = []
    for i in range(0, numpy.size(bias_list)):
        data = fits.getdata(f'{bias_list[i]}')
        bias_images.append((data).astype('float32'))

        
    # Performs sigma clipping to reduce deviations
    bias_images_masked = sigma_clip(bias_images, cenfunc='median', sigma=3, axis=0)

    # Collapses sigma clipped images into single median bias image.
    output = numpy.ma.median(bias_images_masked, axis=0)

    
    # Sets up the fits file write to using the given median_bias_filename
    primary = fits.PrimaryHDU(data=output.data, header=fits.Header())
    hdul = fits.HDUList([primary])
    hdul.writeto(median_bias_filename, overwrite=True)


    return output
