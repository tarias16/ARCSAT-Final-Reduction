#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Filename: reduction.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)
from bias import create_median_bias
from darks import create_median_dark
from flats import create_median_flat, plot_flat
from science import reduce_science_frame
from photometry import do_aperture_photometry
from ptc import calculate_gain, calculate_readout_noise
from timeseries import plot_timeseries, find_period
from photutils.detection import DAOStarFinder
from astropy.stats import sigma_clipped_stats
import os
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.timeseries import TimeSeries

def run_reduction(data_dir, reduced_dir, reduce_images = True):
    
    # Establishes list of file names from data_dir
    # Sorts files into biases, darks, flats, and science images
    cwd = os.getcwd()
    path = os.path.abspath(data_dir)
    files = os.listdir(data_dir)
    for i in range(0, len(files)):
        files[i] = os.path.join(path, files[i])

    dark = 'Dark'
    flat = 'flat'
    bias = 'Bias'
    biases = []
    darks = []
    flats = []
    science_list = []

    for i in range(0, np.size(files)):
        if bias in files[i]:
            biases.append(os.path.abspath(files[i]))
        elif dark in files[i]:
            darks.append(os.path.abspath(files[i]))
        elif flat in files[i]:
            flats.append(os.path.abspath(files[i]))
        else:
            science_list.append(os.path.abspath(files[i]))


    biases = sorted(biases)
    darks = sorted(darks)
    flats = sorted(flats)
    science_list = sorted(science_list)
    
    

    # Gets the files for calculating gain and readout noise
    flat_files_gain = []
    for i in range(1,3):
        flat_files_gain.append(flats[i])
    
    bias_files_rn = []
    for i in range(1,3):
        bias_files_rn.append(biases[i])

    
    # Creates the calibration images and saves their paths

    median_bias = create_median_bias(biases, 'median_bias.fits' )

    bias_path = os.path.abspath('median_bias.fits')

    median_dark = create_median_dark(darks, bias_path, 'median_dark.fits')

    dark_path = os.path.abspath('median_dark.fits')

    median_flat = create_median_flat(flats, bias_path, 'median_flat.fits',
                                    dark_path)

    flat_path = os.path.abspath('median_flat.fits')

    # Plots the median_flat file

    plot_flat(flat_path)


    # Creates lists for the reduced science images and their paths
    reduced_science_list = []
    science_paths = []
    # Creates a path to the specified reduced_dir where reduced images will be saved
    path = os.path.abspath(reduced_dir)

    os.chdir(reduced_dir)

    # Checks to see if the function is called to Reduce Images
    if reduce_images is True:

        # Reduced every image provided in science_list
        for ii in range(0,len(science_list)):
            data = reduce_science_frame(science_list[ii], 
                                bias_path, 
                                flat_path, 
                                dark_path,
                                f'reduced_science{ii+1}.fits')
            reduced_science_list.append(data)
            
            science_paths.append(os.path.join(path, f'reduced_science{ii+1}.fits'))
    else:
        # Gathers reduced science data and paths if reduced images are already provided
        for ii in range(0,len(science_list)):
            reduced_science_list.append(fits.getdata(f'reduced_science{ii+1}.fits'))
            science_paths.append(os.path.abspath(f'reduced_science{ii+1}.fits'))

    # Changes back to original directory
    os.chdir(cwd)

    # Calls the gain calculation function
    gain = calculate_gain(flat_files_gain)
    # Calls the readout noise calculation function
    readout_noise = calculate_readout_noise(bias_files_rn, gain)

    # Defines lists for light curve generation
    times = []
    fluxes = []
    error = []
    
    # Start of the photometry loop for every provided science image
    # Gathers the position of target star in each image
    for i in range(0, len(reduced_science_list)):
        radii = [15]
    
        mean, median, std = sigma_clipped_stats(reduced_science_list[i])
        

        dao = DAOStarFinder(fwhm = 3, threshold = 5*std)

        dets = dao(reduced_science_list[i])
        dets = dets[(dets['xcentroid'] > 470) & (dets['xcentroid'] < 530) & (dets['ycentroid'] > 470) & (dets['ycentroid'] < 530)]
        
        x_pos = dets['xcentroid'].tolist()
        y_pos = dets['ycentroid'].tolist()

        # x = x_pos[(dets['xcentroid'] > 480) & (dets['xcentroid'] < 520)]
        # y = y_pos[(dets['ycentroid'] > 480) & (dets['ycentroid'] < 520)]
        positions = list(zip(x_pos, y_pos))
        
    
    
        # Performs aperture photometry on all images
        ap_phot = do_aperture_photometry(science_paths[i], positions, radii, 30, 5, gain)
        
        # Appends the time, flux, and associated uncertainty of each observation
        times.append(fits.getheader(science_list[i])['DATE-OBS']) 
        fluxes.append(ap_phot['aperture_sum_0'][0]) 
        error.append(ap_phot['aperture_sum_err_0'][0]) 
    
    # Gathers the percent uncertainty in flux based on the ap_phot error
    percent_unc = []
    for i in range(0, len(fluxes)):
        percent_unc.append(error[i]/fluxes[i])


    # Normalizes the flux with respect to its maximum value
    norm_flux = [float(i)/max(fluxes) for i in fluxes]

    # Initializes a table of the times, fluxes, and uncertainty for each observation
    table = Table([times, norm_flux, percent_unc], names = ('time', 'flux', 'uncertainty'))

    # Turns the table into a TimeSeries object
    ts = TimeSeries(table)
    # Writes the table as a CSV file
    ts.write('AE Ursa Majoris Time Series.csv', format = 'csv', overwrite = True)

    # Plots the light curve based on the Timeseries table
    plot_timeseries('AE Ursa Majoris Time Series.csv')

    # Calculates the period of the plotted light curve
    find_period('AE Ursa Majoris Time Series.csv')
    
    return 
