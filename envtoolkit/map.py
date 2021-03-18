
"""
Functions/classes relative to maps
"""

import numpy as np
import pylab as plt
import matplotlib as mpl


def lonflip(lonname, data):

    """
    Reorders a longitude and a data array about the central longitude 
            
    Adapted from the NCL `lonFlip
    <https://www.ncl.ucar.edu/Document/Functions/Contributed/lonFlip.shtml>`_
    function.

    Improved from https://stackoverflow.com/questions/53345442/about-changing-longitude-array-from-0-360-to-180-to-180-with-python-xarray

    :param numpy.array lon: 1D longitude array. 
     The number of longitudes must be even

    :param numpy.array data: the data array. 
     Can be 1D, 2D, 3D or 4D. Longitude must be the 
     last dimension

    :return: A tuple with the first element containing the flipped
     longitude array, and the second element containing the
     flipped data array

    :rtype: tuple

    """

    lon = data[lonname]
    if(lon.min() < 0):
        print('Conversion from -180/180 to 0/360')
        # equivalent to lon[lon < 0] += 360
        lon = (lon + 360) % 360
    else:
        print('Conversion from 0/360 to -180/180')
        # equivalent to lon[lon > 180] -= 360
        lon =  (lon + 180) % 360 - 180
    
    data[lonname] = lon
    dataout = data.sortby(data[lonname])

    return(dataout)




def inpolygon(xin_2d, yin_2d, x_pol, y_pol):

    """ 
    Determines whether points of a 2D-grid are within a polygon. 

    Equivalent to the inpolygon function of Matlab. 
    
    .. note:: If the input polygon is not closed, it is automatically closed

    :param numpy.array xin_2d: 2-D array with the x-coords of the domain
    :param numpy.array yin_2d: 2-D array with the y-coords of the domain
    :param numpy.array x_pol: 1-D array with the x-coords of the polygon
    :param numpy.array y_pol: 1-D array with the y-coords of the polygon
    
    :return: A 2D array (same shape as xin_2d and yin_2d) 
     with 1 when the point is within the polygon, else 0.
    :rtype: numpy.array
    
    """

    from matplotlib import path

    if (xin_2d.ndim!=2):
        raise ValueError("The xin_2d argument must be 2D. %d dimensions" %xin_2d.ndim)
    
    if (yin_2d.ndim!=2):
        raise ValueError("The yin_2d argument must be 2D. %d dimensions" %yin_2d.ndim)
    
    if (x_pol.ndim!=1):
        raise ValueError("The x_pol argument must be 1D. %d dimensions" %x_pol.ndim)
    
    if (y_pol.ndim!=1):
        raise ValueError("The y_pol argument must be 1D. %d dimensions" %y_pol.ndim)

    x_pol = np.array(x_pol)
    y_pol = np.array(y_pol)

    # If the polynom is not closed, we close it
    if (x_pol[0] != x_pol[-1]) | (y_pol[0] != y_pol[-1]):
        x_pol = np.append(x_pol, x_pol[0])
        y_pol = np.append(y_pol, y_pol[0])

    nx, ny = xin_2d.shape

    # creation the input of the path.Path command:
    # [(x1, y1), (x2, y2), (x3, y3)]
    path_input = [(xtemp, ytemp) for xtemp, ytemp in zip(x_pol, y_pol)]

    # initialisation of the path object
    temppath = path.Path(path_input)

    # creation of the list of all the points within the domain
    # it must have a N x 2 shape
    list_of_points = np.array([np.ravel(xin_2d), np.ravel(yin_2d)]).T

    # Calculation of the mask (True if within the polygon)
    mask = temppath.contains_points(list_of_points)

    # reconverting the mask into a nx by ny array
    mask = np.reshape(mask, (nx, ny))

    return mask