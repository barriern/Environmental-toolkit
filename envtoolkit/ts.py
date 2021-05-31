
""" Functions/classes relative to time-series """

import datetime 
import numpy as np
import pylab as plt
import scipy.stats as stats


class Lanczos(object):

    """
    Class for Lanczos filtering. Inspired from 
    NCL's `filwgts_lanczos <http://www.ncl.ucar.edu/Document/Functions/Built-in/filwgts_lanczos.shtml>`_ and `wgt_runave <http://www.ncl.ucar.edu/Document/Functions/Built-in/wgt_runave.shtml>`_ functions.

    :param str filt_type: The type of filter ("lp"=Low Pass, "hp"=High Pass,
     "bp"=Band Pass
    :param int nwts: Number of weights (must be an odd number)
    :param float pca: First cut-off period (low-frequency for band-pass, higher period)
    :param float pcb: Second cut-off period (only for band-pass filters, high-frequency, lower period)
    :param float delta_t: Time-step
    :param float nsigma: A scalar indicating the power of the sigma factor (nsigma >= 0). Note: nsigma=1. is common. 

    """

    def __init__(self, filt_type, nwts, pca, pcb=None, delta_t=1.0, nsigma=1):

        """ Initialisation of the filter """

        self.filt_type = filt_type
        self.nwts = nwts
        self.pca = pca
        self.pcb = pcb
        self.delta_t = delta_t
        self.nsigma = nsigma
        
        if self.pcb is not None:
            if(1 / self.pcb < 1 / self.pca):
                # if periods not well defined, swap them
                print('PCA and PCB have been swapped')
                self.pca, self.pcb = self.pcb, self.pca

        if self.nwts % 2 == 0:
            raise ValueError('Number of weigths must be odd')

        # Because w(n)=w(-n)=0, we would have only nwts-2
        # effective weight. So we add to weights so as to get rid off that
        nwts = self.nwts + 2
        weights = np.zeros((nwts))
        nw = (nwts - 1) // 2  # index of the central value of the filter
                
        if self.filt_type not in ['hp', 'bp', 'lp']:
            raise ValueError('Unknowm filter %s must be "lp", "hp" or "bp"' % filt_type)
    
        weights = self._get_lp_weights(self.pca)
        
        if self.filt_type == 'hp':
    
            weights[0] = 1 - weights[0]
            index = np.arange(1, nw + 1, dtype=int)
            weights[index] = -weights[index]

        elif self.filt_type == 'bp':
    
            if pcb is None:
                raise ValueError("pcb is None but filter is band pass")
                 
            # copy the weights for the pca frequency
            index = np.arange(0, nw + 1, dtype=int)
            
            work = weights[index]
                        
            weights = self._get_lp_weights(self.pcb)
        
            weights[index] -= work
            
        # make weights symetric
        index = np.arange(0, nw, dtype=int)
        work = weights[index + 1]  # copy the values of the weights except central value
        
        weights[nw] = weights[0]  #  define center value

        weights[index] = work[::-1]
        weights[nw+1:] = work
        
        self.wgt = weights[1:-1]
        
            
    def _get_lp_weights(self, pca):

        ''' Computes the normalized weights for a low pass filter (DFILWTQ fortran routine in NCL) '''
        
        # conversion from period to frequency
        fca = self.delta_t / pca
        
        nwts = self.nwts + 2
        nw = (nwts - 1) // 2 
        arg = 2 * np.pi * fca
        
        output = np.zeros((nwts))
        
        output[0] = 2.0 * fca
        
        index = np.arange(1, nw + 1, dtype=float)
        sinx = np.sin(arg * index) / (np.pi * index)
        siny = nw * np.sin(index * np.pi / nw)/ (index * np.pi)
        output[index.astype(int)] = sinx * siny**self.nsigma
        
        # normalize weights
        total = output[0] + 2*np.sum(output[index.astype(int)])
        output /= total
        
        return output    
    
    def wgt_runave(self, data):

        """ Compute the running mean of a ND input array using the filter weights.

        :param numpy.array data: Array to filter out (time must be the first dimension)

        """

        # we retrive the wgt array and initialise the output
        wgt = self.wgt
        output = np.zeros(data.shape)
    
        nwt = len(wgt)        
        nwgt2 = (nwt - 1) // 2
    
        index = np.arange(nwt, dtype=int)
        istart = nwgt2
        iend = data.shape[0] - 1 - nwgt2
        
        for i in range(istart, iend + 1):
            output[i] = np.sum(data[index] * wgt, axis=0)
            index += 1
            
        output[output == 0] = np.nan

        return output


def compute_monthly_clim(data, yyyymm):

    """ 
    Compute a monthly seasonal cycle from a monthly dataset. 

    This script is largely inspired from the 
    `NCL clmMonTLL <https://www.ncl.ucar.edu/Document/Functions/Contributed/clmMonTLL.shtml>`_ function

    :param numpy.ndarray data: dataset whose seasonal 
     cycle to extract (time must be the first dimension)

    :param int date: the date vector (format YYYYMM, 
     with YYYY=year, MM=month). For instance obtained with the 
     :py:func:`envtoolkit.ts.make_yymm` function.

    :return: a numpy array that contains the monthly seasonal cycle. 
     Same dimensions as `data`, except for 
     the first dimension which is 12

    :rtype: numpy.array

    """

    data = np.ma.array(data, mask=data != data)

    yyyymm = yyyymm.astype(np.int)
    year = yyyymm/100
    month = yyyymm-100*year

    monthvec = np.unique(month)

    clim = np.empty(tuple(data[:len(monthvec)].shape))

    for indmonth in range(0, len(monthvec)):

        index = np.nonzero(month == monthvec[indmonth])[0]
        clim[indmonth] = np.mean(data[index], axis=0)

    return clim


def compute_monthly_anom(data, yyyymm, clim):

    """ 
    Computes anomalies relative to a monthly climatology.
    This script is largely inspired from NCL's
    `calcMonAnomTLL <https://www.ncl.ucar.edu/Document/Functions/Contributed/calcMonAnomTLL.shtml>`_ function

    :param numpy.array data: data from which to 
     extract monthly anom (data must be first dimension).

    :param int date: the date vector (format YYYYMM, 
     with YYYY=year, MM=month)

    :param numpy.array dataclim: the montly climatology, computed with the 
     :func:`envtoolkit.ts.compute_clim` function.

    :return: a numpy array that contains the monthly anomalies. 
     Same dimensions as `data`

    :rtype: numpy.array

    """

    data = np.ma.array(data, mask=data != data)

    anoms = np.ma.zeros(data.shape)
    yyyymm = yyyymm.astype(np.int)
    year = yyyymm/100
    month = yyyymm - 100*year

    monthvec = np.unique(month)
    reshape_array = np.ones(data.ndim)

    for indmonth in range(0, len(monthvec)):

        ind_ok = np.nonzero(month == monthvec[indmonth])[0]
        reshape_array[0] = len(ind_ok)
        anoms[ind_ok] = data[ind_ok] - np.tile(clim[indmonth], reshape_array)

    return anoms


def day_of_year(yymmdd):

    """ 
    Returns the day of year of 

    a date vector in format YYYYMMDD

    :param int date: Input date (format YYYYMMDD, 
     with YYYY=year, MM=month, DD=day)

    :return: a numpy.array containing the day 
     of year (1st of January=1, etc)

    :rtype: int

    """

    yymmdd = np.array(yymmdd)

    year = yymmdd/10000
    month = (yymmdd - 10000*year)/100
    day = (yymmdd - 10000*year - 100*month)

    # Computation of the 'day_of_year' vector
    dayofyear = np.zeros(len(year))
    for indtime in range(0, len(year)):
        dateint = datetime.datetime(year[indtime], month[indtime], day[indtime])
        dayofyear[indtime] = dateint.strftime('%j')
    dayofyear = dayofyear.astype(np.int)

    return dayofyear


def compute_daily_clim(data, date, smooth=True, nharm=2):

    """ 
    Computes a daily seasonal cycle from a daily dataset. 
    The user has the possibility to smooth the seasonal cycle, 
    for instance by keeping the first two harmonics (annual and
    semi annual cycles). 
    This script is largely inspired from NCL's
    `calcDayAnomTLL <https://www.ncl.ucar.edu/Document/Functions/Contributed/calcMonAnomTLL.shtml>`_ and
    `smthClmDayTLL <https://www.ncl.ucar.edu/Document/Functions/Contributed/smthClmDayTLL.shtml>`_ functions
    
    :param numpy.ndarray data: dataset whose seasonal 
     cycle to extract (time must be the first dimension)

    :param int date: the date vector (format YYYYMMDD, 
     with YYYY=year, MM=month, DD=day). Obtained by using
     the :py:func:`envtoolkit.ts.make_yymmdd` function.

    :param bool smooth: defines whether the seasonal cycle should 
     be smoothed by using Fast Fourier Transform

    :param int nharm: number of harmonics to retain if `smooth` is True. Use 
     2 to keep annual and semi-annual harmonics

    :return: a numpy array that contains the daily seasonal cycle. 
     Same dimensions as `data`, except for 
     the first dimension which is 366.

    :rtype: numpy.array

    """

    data = np.ma.array(data, mask=data != data)
    dayofyear = day_of_year(date)

    clim = np.empty(tuple(data[:366].shape))
    for indday in range(0, 365):
        i = np.nonzero(dayofyear == indday+1)[0]
        clim[indday] = np.mean(data[i], axis=0)
    clim[-1] = 0.5*(clim[0] + clim[364])

    # smoothing using real fft decomposition
    if smooth:
        clim = smooth_data_fft(clim, nharm=nharm)

    return clim

def smooth_data_fft(data, nharm):

    """ 
    Smooth an input data array by using a FFT filter. 
    
    Inspired from NCL's 
    `smthClmDayTLL<https://www.ncl.ucar.edu/Document/Functions/Contributed/smthClmDayTLL.shtml>` _ function

    The FFT coefficients are extracted, and the 
    signal is reconstructed by using only the 
    first `nharm` harmonics. Is used in the 
    :py:mod:`envtoolkit.ts.compute_daily_clim` function to 
    smooth a daily seasonal cycle.

    :param numpy.ndarray data: dataset whose seasonal 
     cycle to extract (time must be the first dimension)
        
    :param int nharm: number of harmonics to retain
    
    """
    
    fft = np.fft.rfft(data, axis=0)
    fft[nharm] = 0.5*fft[nharm]
    fft[nharm+1:] = 0
    dataout = np.fft.irfft(fft, axis=0)

    return dataout


def compute_daily_anom(data, date, clim):

    """ Computes daily anomalies relative to a daily climatology. Inspired from
    NCL's
    `calcDayAnomTLL <https://www.ncl.ucar.edu/Document/Functions/Contributed/calcMonAnomTLL.shtml>`_ function

    :param numpy.array data: data from which to extract daily 
     anom (time must be the first dimension)

    :param int date: the date vector (format YYYYMMDD, 
     with YYYY=year, MM=month, DD=day). 

    :param numpy.array clim: the daily climatology, computed with the 
     :py:func:`envtoolkit.ts.compute_clim` function

    :return: a numpy array that contains the daily anomalies. 
     Same dimensions as `data`

    :rtype: numpy.array

    """

    if np.any(data[0].shape != clim[0].shape):
        raise ValueError('Leftmost dimensions of data and dataclim ' +
                         ' must be the same. Data = %s and Dataclim = %s'
                         % (data[0].shape, clim[0].shape))

    data = np.ma.array(data, mask=data != data)
    dayofyear = day_of_year(date)
    anoms = np.ma.zeros(data.shape)

    reshape_array = np.ones(data.ndim)

    for indday in range(0, 366):

        ind_ok = np.nonzero(dayofyear == indday+1)[0]

        reshape_array[0] = len(ind_ok)
        anoms[ind_ok] = data[ind_ok] - np.tile(clim[indday], reshape_array)

    return anoms


def make_yymmdd(date):

    """
    Converts a list/array of dates into YYYYMMDD integers.
    
    | YYYY=year
    | MM=month
    | DD=day

    :param list date: A list/array of :py:class:`datetime.datetime` objects
    :return: A numpy array containg the dates in format YYYYMMDD
    :rtype: numpy.array
    """

    year = np.array([d.year for d in date])
    month = np.array([d.month for d in date])
    day = np.array([d.day for d in date])

    return 10000*year + 100*month + day


def make_yymm(date):

    """
    Converts a list/array of date objects into YYYYMM integers.

    | YYYY = year
    | MM = month

    :param list date: A list/array of :py:mod:`datetime.datetime`
    :return: A numpy array containing the dates in format YYYYMM
    :rtype: numpy.array
    """

    year = np.array([d.year for d in date])
    month = np.array([d.month for d in date])

    return 100*year + month


def make_monthly_means(data, yymmdd):

    """
    Computes monthy means from daily values

    :param numpy.array yymmdd: Dates in format YYYYMMDD

    :param numpy.array data: Dataset (time must be first dimension)

    :return: a tuple, containing the monthly output
     dates (format YYYYMM), and the monthly data output.

    :rtype: tuple

    """

    yymmdd = np.array(yymmdd)
    dataout = []
    yyyymm = yymmdd/100

    for yyyint in np.unique(yyyymm):
        iok = np.nonzero(yyyymm == yyyint)[0]
        dataout.append(np.mean(data[iok], axis=0))

    dataout = np.array(dataout)
    return np.unique(yyyymm), dataout


def doy_to_date(year, doy):

    """ 
    Converts from day of year into date

    :param int year: Year 
    :param int doy: Day of year
    """

    return datetime.datetime(year, 1, 1) + datetime.timedelta(doy - 1)

def remove_mean(xdata):

    r""" 
    Remove the mean from the dataset.
    
    .. math::
        
        output = \frac{x-\overline{x}}

    Made on the first dimension (usually time).

    :param numpy.array xdata: Input data
    
    :return: A numpy array with the anomalies of the input array
    
    :rtype: numpy.array

    """

    # masking of data where nan
    xdata = np.ma.masked_where(xdata!=xdata, xdata)

    # calculation of mean/std along first dimension (time)
    xmean = np.mean(xdata, axis=0)
    
    return xdata - xmean


def standardize(xdata, ddof=0):

    r""" 
    Standardizes the xdata array
    
    .. math::
        
        output = \frac{x-\overline{x}}{\sigma_x}

    Made on the first dimension (usually time).

    :param numpy.array xdata: Input data
    
    :param ddof int: Number of degrees of freedom to remove
        for the calculation of :math:`\sigma` (see the 
        :py:func:`numpy.std` documentation.

    :return: A numpy array with the standardized values
        of the input array
    
    :rtype: numpy.array

    """

    # masking of data where nan
    xdata = np.ma.masked_where(xdata!=xdata, xdata)

    # calculation of mean/std along first dimension (time)
    xstd = np.std(xdata, axis=0, ddof=ddof)
    
    output = remove_mean(xdata) / (xstd)

    return output


def xcorr_1d(xdata, ydata, maxlag=None, use_covariance=False, ddof=1):

    '''
    Computes the cross-correlation/cross-covariance
    two one-dimensional arrays. 
    xdata leads at positive lags. Inspired from the 
    :py:func:`pyplot.xcorr` function

    :param xdata numpy.array: x-array (leads at positive lags)
    :param ydata numpy.array: y-array (leads at negative lags)
    :param int maxlag: Number of maximum lag to return. If None,
        maxlag=len(xdata)-1
    :param bool use_covariance: Whether covariance (True)
        or correlation (False) should be returned
    :param int ddof: Number of degrees of freedom to remove if covariance
        is computed. If 1, then covariance is divided by (N-1). 

    :return: A tuple with the lags and the cross-correlation array
    :rtype: tuple

    '''

    # check that the array have the right dimensions

    ntime = len(xdata)

    xdata = np.ma.masked_where(xdata!=xdata, xdata)
    ydata = np.ma.masked_where(ydata!=ydata, ydata)

    xdata = (xdata - xdata.mean())
    ydata = (ydata - ydata.mean())

    output = np.correlate(xdata, ydata, mode=2)

    if use_covariance:
        output /= (len(xdata) - ddof)
    else:
        output /= np.sqrt(np.dot(xdata, xdata) * np.dot(ydata, ydata))

    if maxlag is None:
        maxlag = ntime - 1
    if maxlag > ntime-1:
        maxlag = ntime - 1
    if maxlag < 1:
        maxlag = 1

    lags = np.arange(-maxlag, maxlag + 1)
    index = np.arange(ntime + maxlag -1, ntime -1 - maxlag -1, -1)

    output = output[index]

    return lags, output


def corr_ND(xdata, ydata, use_covariance=False):

    """
    Computes the correlation/covariance between the xdata and ydata
    arrays at 0-lag. Calculation is performed on first dimension 
    (usually time). Loop is optimized by using special loops.
    
    :param xdata numpy.array: x-array (time must be the first dim.)

    :param ydata numpy.array: y-array (time must be the first dim.)
    
    :param bool use_covariance: True if covariance should be computed instead of
        correlation
    
    :return: A tuple with the cross-correlation or cross-covariance array, and the lag array. 
        Cross-correlation has dimensions (``xdata.shape[1:]``, ``ydata.shape[1:]``).

    :rtype: tuple
    """

    if xdata.shape[0] != ydata.shape[0]:
        message = 'The first dimension of xdata array(%d) is different from the first dimension of the ydata array (%d). This program will be stopped.'
        raise ValueError(message)

    ntime = xdata.shape[0]

    # removing the time-mean over each point
    xdata = remove_mean(xdata)
    ydata = remove_mean(ydata)

    # chosing the proper function depending on the
    # use_covariance option: np.cov or np.corrcoef
    if use_covariance:
        corrfunc = np.ma.cov
    else:
        corrfunc = np.ma.corrcoef

    # masking of the data where fill value
    xdata = np.ma.masked_where(xdata != xdata, xdata)
    ydata = np.ma.masked_where(ydata != ydata, ydata)

    # If xdata or ydata is 1d (nt), we convert it into
    # array of dimension (nt, 1)
    if xdata.ndim == 1:
        xdata = np.atleast_2d(xdata).T
    if ydata.ndim == 1:
        ydata = np.atleast_2d(ydata).T

    # extracting the shape of the output from 
    # maxlag and the shapes of the input arrays
    outshape = list(xdata.shape[1:]) + list(ydata.shape[1:])

    # extraction of the number of space dimensions in input arrays
    ndimx = np.prod(xdata.shape[1:])
    ndimy = np.prod(ydata.shape[1:])

    # conversion of the array into a 2D array
    # if their dims is greater than 2
    if xdata.ndim > 2:
        xdata = np.reshape(xdata, (ntime, ndimx))
    if ydata.ndim > 2:
        ydata = np.reshape(ydata, (ntime, ndimy))

    # transposition so that time is the last dimension
    xdata = xdata.T
    ydata = ydata.T

    # calculation of correlation/covariance using special loops
    correlation = [corrfunc(xtemp, ytemp)[0][1] for xtemp in xdata for ytemp in ydata]

    # conversion of corr. into array and reshaping
    correlation = np.array(correlation, dtype=float)
    correlation = np.reshape(correlation, outshape)

    return correlation


def xcorr_ND(xdata, ydata, maxlag=None, use_covariance=False):
    
    """
    Computes the cross-correlation/cross-covariance between the xdata and ydata
    arrays. Calculation is performed on first dimension 
    (usually time). xdata leads for positive lags.

    :param xdata numpy.array: x-array (time must be the first dim.)

    :param ydata numpy.array: y-array (time must be the first dim.)
    
    :param bool use_covariance: True if covariance should be computed instead of
        correlation
    
    :param int maxglag: number of lag to consider. If None, maxlag=ntime
    
    :return: A tuple with the cross-correlation or cross-covariance array, and the lag array. 
        Cross-correlation has dimensions (xdata.shape[1:], ydata.shape[1:], 2*maxlag+1).

    :rtype: tuple

    """
    
    if xdata.shape[0] != ydata.shape[0]:
        message = 'The first dimension of xdata array(%d) is different from the first dimension of the ydata array (%d). This program will be stopped.'
        raise ValueError(message)
    
    ntime = xdata.shape[0]

    if maxlag is None:
        maxlag = ntime - 1
    if maxlag > ntime - 1:
        maxlag = ntime - 1
    if maxlag < 1:
        maxlag = 1
    
    # the correlation function used is the 1D xcorr function
    corrfunc = xcorr_1d
    # setting options of the xcorr function
    if use_covariance:
        dict_xcorr = {"use_covariance": True, "maxlag": maxlag} 
    else:
        dict_xcorr = {"use_covariance": False, "maxlag": maxlag} 

    # If xdata or ydata is 1d (nt), we convert it into
    # array of dimension (nt, 1)
    if xdata.ndim == 1:
        xdata = np.atleast_2d(xdata).T
    if ydata.ndim == 1:
        ydata = np.atleast_2d(ydata).T

    # extracting the spatial dimensions associated with
    # the imput arrays
    ndimx = np.prod(xdata.shape[1:])
    ndimy = np.prod(ydata.shape[1:])
    
    # extracting the shape of the output from 
    # maxlag and the shapes of the input arrays
    tempshape = list(xdata.shape[1:]) + list(ydata.shape[1:])
    outshape = tempshape + [2*maxlag+1]
    
    # conversion of the array into a 2D array
    # if their dims is greater than 2
    if xdata.ndim > 2:
        xdata = np.reshape(xdata, (ntime, ndimx))
    if ydata.ndim > 2:
        ydata = np.reshape(ydata, (ntime, ndimy))

    # transposition so that time is the last dimension
    xdata = xdata.T
    ydata = ydata.T

    # initialises the lag vector
    lag_array = np.arange(maxlag+1)
    
    # calculation of the cross-correlation for each spatial point
    output = [corrfunc(xtemp, ytemp, **dict_xcorr) for xtemp in xdata for ytemp in ydata]
    
    # extraction of the lag array
    lag = output[0][0]

    # extraction of the cross-correlation data 
    output = [var[1] for var in output]

    # conversion into a numpy array and reshaping
    output = np.array(output)
    
    output = np.reshape(output, outshape)
    
    return lag, output

def autolag1(TS1):
    x = TS1
    y = np.empty(x.shape)
    n = len(TS1)
    y[0:n-1]  =x[1:n]
    y[n-1]  =x[0]
    test = np.corrcoef(x,y)[0,1]
    return test


def corr_sig(ts1, ts2, df, coeff, use_bretherton=False):

    n = len(ts1)
    Nstep = 1
    maxlag = n-1
    vectemps = np.arange(0, n)
    a = autolag1(ts1)
    b = autolag1(ts2)
    bretherton = (1 - a*b) / (1 + a*b)

    if not use_bretherton:
        bretherton = 1.0

    rsign = np.empty(2 * maxlag  +1)

    cpt = 0
    for lag in range(-maxlag, maxlag+1, 1):

        veclag = vectemps + lag

        I1 = np.nonzero(veclag >= 0)[0]
        I1 = I1[0]

        I2 = np.nonzero(veclag <= n-1)[0]
        I2 = I2[-1]

        nptcom = I2 - I1 + 1

        dl = (nptcom - df) * bretherton

        rsign[cpt] = coeff/(np.sqrt(dl+coeff*coeff))
        cpt += 1

    veclag = np.arange(-maxlag,maxlag+1, 1)
    return rsign,veclag
