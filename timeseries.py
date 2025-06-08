import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.timeseries import TimeSeries, LombScargle
import numpy as np
from astropy.io import ascii
from astropy import units as u

def plot_timeseries(timeseries):

    # Reads in the provides timeseries table
    table = ascii.read(timeseries)

    # Converts the table into a TimeSeries object
    ts = TimeSeries(table)

    # Plots the Light Curve
    fig, ax = plt.subplots( 
        constrained_layout = True
        )

    ax.set_title('AE Ursa Majoris Light Curve')
    ax.set_xlabel('Julian Date')
    ax.set_ylabel('Normalized Flux')

    # Sets up associated error bars as uncertainty
    ax.errorbar(ts.time.jd, ts['flux'], (ts['flux']*ts['uncertainty']), ls = 'none', 
                elinewidth = 0.5,
                capsize = 1,
                ecolor = 'b'
                )
    ax.plot(ts.time.jd, ts['flux'], 'k.', markersize=1)

    # ax.legend(loc=0);
    
    fig.show()

    # Saves the figure 
    fig.savefig('AE Ursa Majoris Time Series.png')

    return




def find_period(timeseries):

    # reads in the provided table
    table = ascii.read(timeseries)

    # Turns table into a masked table
    table_masked = Table(table, masked = True, copy = False)

    # Turns table into a TimeSeries object
    ts = TimeSeries(table_masked)

    ts_unmasked = TimeSeries(table)

    # Masks out outlier data points with flux lower than 0.5
    ts['flux'].mask = ts['flux'] < 0.5

    # Produces a periodogram with a 1 term Fourier Series
    frequency1, power1 = LombScargle.from_timeseries(ts, 'flux', 'uncertainty', normalization = 'model', nterms = 1).autopower(minimum_frequency= 1/(3*u.h), maximum_frequency= 1/(0.2*u.h))
    ls1 = LombScargle.from_timeseries(ts, 'flux', 'uncertainty', normalization = 'model', nterms = 1)

    # Produces a periodogram with a 2 ter Fourier Series
    frequency2, power2 = LombScargle.from_timeseries(ts, 'flux', 'uncertainty', normalization = 'model', nterms = 2).autopower(minimum_frequency= 1/(3*u.h), maximum_frequency= 1/(0.2*u.h))
    ls2 = LombScargle.from_timeseries(ts, 'flux', 'uncertainty', normalization = 'model', nterms = 2)

    # Overplots the LombScargle Power vs frequency for the 1 and 2 term solutions
    fig, ax = plt.subplots( 
        constrained_layout = True)
    
    ax.plot(1/frequency1, power1, label = '1 Term Fourier Fit', color = 'b')

    ax.plot(1/frequency2, power2, label = '2 Term Fourier Fit', color = 'r')

    ax.set_title('AE Ursa Majoris Primary Period')

    ax.set_xlabel('Period[Hours]')
    ax.set_ylabel('LombScargle Power')

    ax.legend(loc = 0)
    fig.savefig('Period Constraint.png')
    fig.show()

    # Produces a best fit curve for the 1 term solution
    best_freq1 = frequency1[np.argmax(power1)]
    fit1 = ls1.model(ts['time'], best_freq1)

    # Produces a best fit curve for the 2 term solution
    best_freq2 = frequency2[np.argmax(power2)]
    fit2 = ls2.model(ts['time'], best_freq2)

    # Plots the 1 term solution with the light curve
    fig, ax = plt.subplots( 
        constrained_layout = True
        )

    ax.set_title('AE Ursa Majoris Fitted Light Curve')
    ax.set_xlabel('Julian Date')
    ax.set_ylabel('Normalized Flux')

    # Plots error bars
    ax.errorbar(ts_unmasked.time.jd, ts['flux'], (ts['flux']*ts['uncertainty']), ls = 'none', 
                elinewidth = 0.5,
                capsize = 1,
                ecolor = 'b',
                )
    

    ax.plot(ts_unmasked.time.jd, ts['flux'], 'k.', markersize=2,
            label = 'Stellar Light Curve',
            )
    ax.plot(ts.time.jd, fit1,
            label = '1 Term Fit'
            )
    

    ax.plot(ts.time.jd, fit2,
            label = '2 Term Fit')

    ax.legend(loc = 2)

    fig.savefig('Light Curve Fit.png')
    fig.show()
    



    return


