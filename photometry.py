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
    """This function must:

    - Accept a fully reduced science image as a file and read it.
    - Accept a list of positions on the image as a list of tuples (x, y).
    - Accept a list of aperture radii as a list of floats.
    - Accept a the radius at which to measure the sky background as sky_radius_in.
    - Accept a the width of the annulus as sky_annulus_width.
    - For each position and radius, calculate the flux in the aperture, subtracting
      the sky background. You can do this any way that you like but probably you'll
      want to use SkyCircularAnnulus from photutils.
    - The function should return the results from the aperture photometry. Usually
      this will be an astropy table from calling photutils aperture_photometry, but
      it can be something different if you use a different library.

    Note that the automated tests just check that you are returning from this
    function, but they do not check the contents of the returned data.

    """
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

    # background_array = Background2D(science_image, 1)

    error = calc_total_error(science_image, bckgrd[0], gain)

    # Performs aperture photometry for every star and radius
    ap_phot = aperture_photometry(science_image, ap, error = error)


    # Subtracts the background sky flux from each star at each radius in the aperture photometry table
    for i in range(0, np.size(radii)):
        ap_phot[f'aperture_sum_{i}'] = ap_phot[f'aperture_sum_{i}'] - sky_flux[i]

    # Adds the local background of each star to the ap_phot table
    ap_phot['background'] = bckgrd

    return ap_phot



def plot_radial_profile(aperture_photometry_data, image, output_filename="radial_profile.png"):
    """This function must:

    - Accept a table of aperture photometry data as aperture_photometry_data. This
      is probably a photutils table, the result of do_aperture_photometry, but you
      can pass anything you like. The table/input data can contain a single target
      or multiple targets, but it must include multiple apertures.
    - Plot the radial profile of a target in the image using matplotlib. If you
      have multiple targets, label them correctly.
    - Plot a vertical line at the radius of the sky aperture used for the photometry.
    - Save the plot to the file specified in output_filename.

    """
    # Calls the reduced science image used in do_aperture_photometry
    reduced_science = fits.getdata(image)

    # Establishes a list of radii to plot radial profiles
    length = len(aperture_photometry_data[1])
    radii = list(range(1,length-4))

    # Defines a radial profile for a signle star in reduced_science while subtracting that star's local background
    rp = RadialProfile((reduced_science - aperture_photometry_data['background'][1]), (aperture_photometry_data['xcenter'][1], aperture_photometry_data['ycenter'][1]), radii)

    # Records the radius where the most flux is captured
    max_flux = 0
    
    for i in range(0,np.size(radii)):
        if aperture_photometry_data[f'aperture_sum_{i}'][1] > max_flux:
            max_flux = aperture_photometry_data[f'aperture_sum_{i}'][1]
            max_radii = i

    # Plots the radial profile of the star and marks the radius at which the most signal without background is found
    plt.clf()        
    plt.plot(rp.radius, rp.profile, label = f'Position ({aperture_photometry_data['xcenter'][1]:.2f}, {aperture_photometry_data['ycenter'][1]:.2f})')
    
    plt.axvline(x = max_radii + 1, color = 'r', label = f'Aperture = {max_radii + 1} pixels', linestyle = '--')
    
    plt.legend(loc=0, shadow =True)

    plt.savefig(output_filename)

    pass
