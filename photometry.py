#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Filename: photometry.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)
from photutils.detection import DAOStarFinder
from photutils.aperture import CircularAperture, aperture_photometry, ApertureStats, CircularAnnulus, SkyCircularAperture, SkyCircularAnnulus
from photutils.centroids import centroid_quadratic
from photutils.profiles import RadialProfile
from photutils.centroids import centroid_1dg
from photutils.datasets import load_star_image
from photutils.utils import calc_total_error
from photutils.background import Background2D
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from science import reduce_science_frame

def do_aperture_photometry(
    image,
    positions,
    radii,
    sky_radius_in,
    sky_annulus_width,
    gain
):

    # Gets reduced science image from image parameter
    science_image = fits.getdata(image)


    # Defines a list of Circular apertures for every radius provided
    ap = [CircularAperture(positions, r) for r in radii]
    

    # Defines an annular region for local background calculation
    annulus = CircularAnnulus(positions, r_in = sky_radius_in, r_out = (sky_radius_in + sky_annulus_width))

    # Gathers stats form annular region
    ap_stats = ApertureStats(science_image, annulus)

    # Defines the local background for each star in the image
    bckgrd = []
    for i in range(0, np.size(ap_stats)):
        bckgrd.append(ap_stats[i].median)

    # Defines a blank 2D array to keep track of background flux 
    rows, columns = (np.size(radii), np.size(bckgrd))
    sky_flux = [[0 for _ in range(columns)] for _ in range(rows)]

  
    # Calculates the background flux for star j with aperture radius i
    for j in range(0, np.size(bckgrd)):
        for i in range(0,np.size(radii)):
            area = (ap[i].area_overlap(science_image))
            sky_flux[i][j] = (bckgrd[j] * area[j])

    # Calculates the error from CCD gain, and local background in each image
    error = calc_total_error(science_image, bckgrd[0], gain)

    # Performs aperture photometry for every star and radius, propogates the error through the photometry analysis
    ap_phot = aperture_photometry(science_image, ap, error = error)


    # Subtracts the background sky flux from each star at each radius in the aperture photometry table
    for i in range(0, np.size(radii)):
        ap_phot[f'aperture_sum_{i}'] = ap_phot[f'aperture_sum_{i}'] - sky_flux[i]

    # Saves the local background of each star to the ap_phot table
    ap_phot['background'] = bckgrd

    return ap_phot

