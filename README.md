# ARCSAT-Final-Reduction

This repository contains the files used for the final reduction of my Astro 480 ARCSAT Data Reduction Project. 

The script reduction.py is the main script to run and calls all of the other scripts within it.

The bias.py, flat.py, dark.py, and science.py scripts are used for the reduction of raw science images. 

The ptc.py and photometry.py characterize the gain, readout noise, and perform photometry on provided science images.

The script timeseries.py contains functions for plotting the lightcurve of an object with photometry data and for estimating the period of the light curve using fourier transforms.

The image files utilized for the reduction were found in the /images/AS06/20250511 direcrotry within ARCSAT, science images were names AEUma_g_, and the calibration images used were from the first half of the night around 2:30 UTC.
