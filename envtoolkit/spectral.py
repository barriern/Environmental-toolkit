
""" Functions/classes relative to time-series """

from datetime import datetime, timedelta
import numpy as np
import spectrum as sp
import pylab as plt


def multitaper(xdata, deltat=1., nbandw=3, nfft=None):

    """
    Computes the power spectrum using the multi-taper method 
    
    This method uses adaptive weighting. Adapted from the
    `ptmtPH.m <http://www.people.fas.harvard.edu/~phuybers/Mfiles/>`_
    script by Peter Huybers

    :param numpy.array xdata: input data vector.
    :param float deltat: sampling interval
    :param int nbandw: time bandwidth product, acceptable values range from 
     0:.5:length(x)/2-1. 2*nbandw-1 dpss tapers are applied 
     except if nbandw=0 a boxcar window is applied 
     and if nbandw = .5 (or 1) a single dpss taper is applied.
    :param int nfft:  number of frequencies to evaluate P at, default is
     length(x) for the two-sided transform.

    :return: a tuple (P, s, ci) with:

        - P: Power spectrum computed via the multi-taper method.
        - s: Frequency vector.
        - ci: 95% confidence intervals.

    :rtype: tuple

    """

    if xdata.ndim != 1:
        raise IOError("The dataset must be one-dimensional. " + 
                "It currently has %d dimensions." %xdata.ndim)

    if nfft is None:
        nfft = len(xdata)

    k = np.min([np.round(2*nbandw), len(xdata)])
    k = np.max([k-1, 1])

    # use of the dpss function of the spectrum library
    taper, eigval = sp.dpss(len(xdata), NW=nbandw, k=k)

    if len(xdata) <= nfft:
        fftcoef = np.abs(np.fft.fft(taper*np.transpose(np.tile(xdata, (k, 1))),
                                    n=nfft, axis=0))**2
    else:
        raise IOError('The case of len(xdata)>nfft is not implemented yet.')

    # Iteration to determine adaptive weights:
    if k > 1:

        xdata = np.mat(xdata).T
        sig2 = xdata.T*xdata/len(xdata)  # power
        spectrum = (fftcoef[:, 0] + fftcoef[:, 1]) / 2.   # initial spectrum estimate
        spectrum_temp = np.zeros(nfft)
        eigval = np.mat(eigval)

        k = np.mat(np.ones((1, k)))

        while np.sum(np.abs(spectrum-spectrum_temp)/nfft) > .0005*sig2/nfft:
            weights = np.array((np.mat(spectrum).T*k) /(np.mat(spectrum).T*eigval+np.mat(np.ones((nfft, 1)))*np.mat(sig2*(1-eigval))))  # weights
            weights = weights**2*np.array(np.mat(np.ones((nfft, 1)))*eigval)
            spectrum_temp = np.sum(weights*fftcoef, axis=1) / np.sum(weights, axis=1)
            spectrum, spectrum_temp = spectrum_temp, spectrum  # swap spectrum and spectrum_tep

        error = 2*np.sum(weights**2*np.array(np.mat(np.ones((nfft, 1)))*eigval), axis=1)**2
        error = error/np.sum(weights**4*np.array(np.mat(np.ones((nfft, 1)))*np.mat(np.array(eigval)**2)), axis=1)

    select = np.arange(0, (nfft+1)/2+1)

    spectrum = spectrum[select]
    freq = np.arange(0, 1/deltat, 1/(nfft*deltat))[select]
    error = error[select]
    error = np.array([1/(1-2/(9*error) - 1.96*np.sqrt(2./(9*error)))**3,
                      1/(1-2/(9*error) + 1.96*np.sqrt(2/(9*error)))**3])

    return spectrum, freq, error


def plot_spectra(xdata, deltat, ferror, spec_type='variance',
                 nbandw=3, **kwargs):

    """
    Plot the power spectrum of a time-series. 
    
    The power-spectrum is computed by using 
    the :py:func:`envtoolkit.spectra.multitaper` function. 
    Adapted from the JD_spectra script of Julie Deshayes.

    :param numpy.array xdata: time series whose spectrum to plot
    :param float deltat: time step of the time series in seconds
    :param float ferror: the frequency where to do draw the errorbar
    :param str spec_type: 'variance' for variance spectrum 
     (spectrum, in $unit^2 cpy^{-1}$), else 
     energy spectrum (freq*spectrim, in $unit^2$).
    :param int nbandw: time bandwidth product (multitaper used 2*nbandw-1 tapers)
    :param str ylabel: label of the yaxis
    :param **kwargs: plot additional arguments
    :return: A tuple (Px,F,Pxc) with:

        - Px: the spectrum vector
        - F: the frequency vector
        - PxC: the errorbar vector

    :rtype: tuple
    """


    xdata = xdata - np.mean(xdata)
    [spectrum, freq, error] = multitaper(xdata, nbandw=nbandw)

    freq = freq*365*86400/deltat    # to get the result in cpy
    spectrum = spectrum/(365*86400/deltat)    # to get the result in cpy^{-1}
    error = error/(365*86400/deltat)    # to get the result in cpy^{-1}
    barmax = error[0, 0]*1e2
    barmin = error[1, 0]*1e2
    
    freq = np.ma.masked_equal(freq, 0)

    if spec_type == 'variance':
        plt.plot(freq, spectrum, **kwargs)
        plt.plot([ferror, ferror], [barmin, barmax], marker='_', **kwargs)

    else:
        plt.plot(freq, freq*spectrum, **kwargs)

    plt.grid(True, which='both', axis='x')
    plt.grid(True, which='major', axis='y')
    plt.gca().set_yscale('log')
    plt.gca().set_xscale('log')

    return spectrum, freq, error


def plot_slope(spectrum, freq, fmin=None, fmax=None,
               offy=0, **kwargs):

    """ Draws the slope of the spectrum.

    :param numpy.array spectrum: frequency spectrum (output of the :py:func:`plot_spectra` function)
    :param numpy.array freq: frequency vector (output of the :py:func:`plot_spectra` function)
    :param axes ax: axes on which to draw the slope
    :param float fmin: the first frequency on which to compute the slope
    :param float fmax: the last frequency on which to compute the slope
    :param float offy: the y-offset of the slope
    :param dict **kwargs: additional plotting arguments
    :return: the value of the slope
    :rtype: float
    """


    # if fmin or fmax are None, we set to the min
    # and max values of freq
    if fmin is None:
        fmin = freq.min()
    if fmax is None:
        fmax = freq.max()

    islope = np.nonzero((freq >= fmin) & (freq <= fmax) & (freq != 0))[0]
    fout = freq[islope]
    pxout = spectrum[islope]

    yslope = np.log(pxout)  # pylint: disable=no-member
    xslope = np.log(fout)   # pylint: disable=no-member

    # calculation of the polyfit
    pol = np.polyfit(xslope, yslope, 1)

    # calculation of the trend
    trend = np.exp(pol[0]*xslope + pol[1] + offy)

    # plotting the trend
    plt.plot(fout, trend, **kwargs)

    return pol[0]


def plot_ref_slope(fmin, fmax, yinterc, slope=2,
                   **kwargs):

    """ Draws reference slopes on a log spectral plot. 

    :param float fmin: frequency where to start the reference the slope
    :param float fmax: frequency where to end the reference the slope
    :param float yinterc: y intercept of the slopes
    :param float kval: the slope that we
    :param dict **kwargs: additional line plotting arguments
    """

    xdata = np.linspace(fmin, fmax, 5)

    yout = np.log(yinterc) - slope*(np.log(fmin) - np.log(xdata))  # pylint:disable=no-member
    yout = np.exp(yout)
    plt.plot(xdata, yout, **kwargs)
    plt.text(xdata[-1], yout[-1], ' k = -' + str(slope),
             ha='center', va='center', color='k',
             bbox=dict(boxstyle="round", fc="0.9"))
