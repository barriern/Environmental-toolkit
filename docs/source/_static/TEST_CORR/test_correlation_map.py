import numpy as np


def cross_correlation(xdata, ydata, use_covariance=False, maxlag=None, stand=True, ddof=0):

    if xdata.shape[0] != ydata.shape[0]:
        message = 'The first dimension of xdata array(%d) is different from the first dimension of the ydata array (%d). This program will be stopped.'
        raise ValueError(message)

    ntime = xdata.shape[0]
    ndimx = np.prod(xdata.shape[1:])
    ndimy = np.prod(ydata.shape[1:])
    
    if maxlag is None:
        maxlag = ntime

    if stand:
        xdata = standardize(xdata, ddof=ddof)
        ydata = standardize(ydata, ddof=ddof)

    # mask the input arrays where nan
    xdata = np.ma.masked_where(xdata!=xdata, xdata)
    ydata = np.ma.masked_where(ydata!=ydata, ydata)


    # extracting the shape of the output from 
    # maxlag and the shapes of the input arrays
    outshape = [2*maxlag+1] + list(xdata.shape[1:]) + list(ydata.shape[1:])
    
    xdata = np.reshape(xdata, (ntime, ndimx))
    ydata = np.reshape(ydata, (ntime, ndimy))

    # initialises the output array
    xleady = np.empty((maxlag+1, ndimx*ndimy), dtype=np.float)
    yleadx = np.empty((maxlag+1, ndimx*ndimy), dtype=np.float)

    # initialises the lag vector
    lag_array = np.arange(maxlag+1)

    # initialisation of the correlation function (avoiding dots)
    # see https://wiki.python.org/moin/PythonSpeed/PerformanceTips#Avoiding_dots...
    if use_covariance:
        corrfunc = np.cov
    else:
        corrfunc = np.corrcoef

    # loop over all the lags
    for lag in lag_array:
    
        cpt = 0

        # loop over the x dimensions
        for indx in range(0, ndimx):

            # loop over the y dimensions
            for indy in range(0, ndimy):
                
                # =================================== Calculation of the ylead correlation
                xtemp = xdata[lag:, indx]
                ytemp = ydata[0:ntime-lag, indy]

                # calculates the mask array: True if both xtemp and ytemp are NOT masked
                mask = (np.ma.getmaskarray(xtemp)==False) & (np.ma.getmaskarray(ytemp)==False) 
                
                # index where to do the correlation calculation
                ind_ok = np.nonzero(mask==True)[0]

                if len(ind_ok)>0:
                    yleadx[lag, cpt] = corrfunc(xtemp[ind_ok], ytemp[ind_ok])[0][1]
                else:
                    yleadx[lag, cpt] = np.nan
                
                # =================================== Calculation of the xleady correlation
                xtemp = xdata[0:ntime-lag, indx]
                ytemp = ydata[lag:, indy]
                
                # calculates the mask array: True if both xtemp and ytemp are NOT masked
                mask = (np.ma.getmaskarray(xtemp)==False) & (np.ma.getmaskarray(ytemp)==False) 
                
                # index where to do the correlation calculation
                ind_ok = np.nonzero(mask==True)[0]
                
                if len(ind_ok)>0:
                    xleady[lag, cpt] = corrfunc(xtemp[ind_ok], ytemp[ind_ok])[0][1]
                else:
                    xleady[lag, cpt] = np.nan 
            
                # iteration of the counter
                cpt += 1

    output = np.empty((2*maxlag+1, ndimx*ndimy), dtype=np.float)
    output[:maxlag, :] = yleadx[:0:-1, :]
    output[maxlag:, :] = xleady

    output = np.reshape(output, outshape)
    lag = np.arange(-maxlag, maxlag+1)

    return output, lag


def fast_cross_corr(xdata, ydata, use_covariance=False, stand=True):

    """
    Calculation of 0-lag correlation. Loop is optimized
    by using special loops.
    """

    if xdata.shape[0] != ydata.shape[0]:
        message = 'The first dimension of xdata array(%d) is different from the first dimension of the ydata array (%d). This program will be stopped.'
        raise ValueError(message)

    # if stand=True, we standardize the data
    if stand:
        xdata = standardize(xdata, ddof=ddof)
        ydata = standardize(ydata, ddof=ddof)

    ntime = xdata.shape[0]

    # extracting the shape of the output from 
    # maxlag and the shapes of the input arrays
    outshape = list(xdata.shape[1:]) + list(ydata.shape[1:])

    ndimx = np.prod(xdata.shape[1:])
    ndimy = np.prod(ydata.shape[1:])

    # conversion of the array into a 2D array
    # and transposition (time=last dimension)
    xdata = np.reshape(xdata, (ntime, ndimx)).T
    ydata = np.reshape(ydata, (ntime, ndimy)).T

    # chosing the proper function:
    # np.cov or np.corrcoef
    if use_covariance:
        corrfunc = np.cov
    else:
        corrfunc = np.corrcoef

    # calculation of correlation/covariance using special loops
    correlation = [corrfunc(xtemp, ytemp)[0][1] for xtemp in xdata for ytemp in ydata]

    # conversion of corr. into array and reshaping
    correlation = np.array(correlation, dtype=float)
    correlation = np.reshape(correlation, outshape)

    return correlation

if __name__ == "__main__":

    from netCDF4 import Dataset
    from pylab import *
    from datetime import datetime

    fin = Dataset('data/data_for_correlation.nc', 'r')
    occu = fin.variables["occu"][:]
    psl = fin.variables["psl"][:]
    lat = fin.variables["lat"][:]
    lon = fin.variables["lon"][:]
    fin.close()

    print psl.shape
    print occu.shape

    psl = np.tile(psl, (1, 2, 22))
    occu = np.tile(occu, (1, 8))
    print psl.shape
    print occu.shape
    
    #    now = datetime.now()
    #    corr1, lag =  cross_correlation(occu, psl, use_covariance=False, maxlag=0, stand=False, ddof=0)
    #    print datetime.now() - now
    #    corr1 = np.squeeze(corr1)
    #    
    #    
    now = datetime.now()
    corr2 = fast_cross_corr(occu, psl)
    print datetime.now() - now
    #
    #    #print corr1.shape
    #    print corr2.shape
    #
    #    indcl = 3
    #    figure()
    #    subplot(2, 1, 1)
    #    cs = pcolormesh(lon, lat, corr1[indcl, :, :])
    #    colorbar(cs)
    
#    subplot(2, 1, 2)
#    cs = pcolormesh(lon, lat, corr2[indcl, :, :])
#    colorbar(cs)
#
#    savefig("psl.png")
