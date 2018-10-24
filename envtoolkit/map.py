
"""
Functions/classes relative to maps
"""

from mpl_toolkits.basemap import Basemap
import numpy as np
import pylab as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib as mpl


def make_bmap(lon, lat, **kwargs):

    """ 
    Creates a basemap object by using the limites of input lon/lat arrays. 

    | Lower left corner coordinates = min(lon), min(lat)
    | Upper right corner coordinates = max(lon), max(lat)

    :param numpy.array lon: longitude array (any dimension)
    :param numpy.array lat: latitude array (any dimension)

    :param dict \**kwargs: keyword arguments that are
     passed to the Basemap class (proj, resolution, etc).

    :return: A Basemap object
    :rtype: mpl_toolkits.basemap.Basemap
    """

    # conversion into masked arrays 
    lon = np.ma.array(lon)
    lat = np.ma.array(lat)

    bmap = Basemap(llcrnrlon=lon.min(), llcrnrlat=lat.min(),
                   urcrnrlon=lon.max(), urcrnrlat=lat.max(),
                   **kwargs)

    return bmap


def lonflip(lon, data):

    """
    Reorders a longitude and a data array about the central longitude 
            
    Adapted from the NCL `lonFlip
    <https://www.ncl.ucar.edu/Document/Functions/Contributed/lonFlip.shtml>`_
    function.

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

    mlon = len(lon)

    if lon.ndim != 1:
        raise ValueError("The longitude argument must be 1D")

    if mlon % 2 != 0:
        raise ValueError("The lonflip function requires that the number \n \
        of longitudes be even. mlon= %i \n =======================" % mlon)

    mlon2 = mlon/2

    if data.ndim == 1:
        temp = np.ma.empty(data.shape)
        temp[0:mlon2] = data[mlon2:]
        temp[mlon2:] = data[0:mlon2]
    elif data.ndim == 2:
        temp = np.ma.empty(data.shape)
        temp[:, 0:mlon2] = data[:, mlon2:]
        temp[:, mlon2:] = data[:, 0:mlon2]
    elif data.ndim == 3:
        temp = np.ma.empty(data.shape)
        temp[:, :, 0:mlon2] = data[:, :, mlon2:]
        temp[:, :, mlon2:] = data[:, :, 0:mlon2]
    elif data.ndim == 4:
        temp = np.ma.empty(data.shape)
        temp[:, :, :, 0:mlon2] = data[:, :, :, mlon2:]
        temp[:, :, :, mlon2:] = data[:, :, :, 0:mlon2]
    else:
        raise ValueError("Dimension %i cannot exceed 4" % data.ndim)

    tlon = np.empty(lon.shape)

    tlon[0:mlon2] = lon[mlon2:]
    tlon[mlon2:] = lon[0:mlon2]

    if lon[0] >= 0:  # (say) 0=>355
        tlon[0:mlon2] = lon[mlon2:] - 360
        tlon[mlon2:] = lon[0:mlon2]
    else:                   # (say) -180=>175
        tlon[0:mlon2] = lon[mlon2:]
        tlon[mlon2:] = lon[0:mlon2] + 360

    return tlon, temp


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


class Lambert(object):

    """
    Class for the handling of Lambert Conic Projection.
    
    It is initialised by providing 
    the bounding longitudes and latitudes of the domain to plot. It 
    may also take the resolution of the 
    :py:class:`mpl_toolkits.basemap.Basemap` class.

    :param float lonmin: minimum longitude
    :param float lonmax: maximum longitude
    :param float latmin: minimum latitude
    :param float latmax: maximum latitude
    :param str resolution: resolution of the map (cf `Basemap`)

    """

    def __init__(self, lonmin, lonmax, latmin, latmax,
                 resolution='l'):

        """ Initialisation of the Lambert class """

        # We set the map limits attributes
        self.lonmin = lonmin
        self.lonmax = lonmax
        self.latmin = latmin
        self.latmax = latmax

        # Definition of the hemisphere attribute
        # it is either NH or SH
        if (self.latmin >= 0) and (self.latmax > 0):
            self.hemisphere = 'NH'
        elif (self.latmin > -90) and (self.latmax <= 0):
            self.hemisphere = 'SH'
        else:
            raise ValueError('latmin and latmax should be the same sign. ' + 
            'Here, latmin=%f and latmax%f' % (latmin, latmax))

        # Initialisation of some map coordinates
        # that will be used for setting the map projection
        if self.hemisphere == 'NH':
            lat2 = 89.999
            lat1 = 0.001
        else:
            lat2 = -89.999
            lat1 = -0.001

        # definition of the center longitude lon0
        lon0 = 0.5*(self.lonmin + self.lonmax)

        # initialisation of a dummy bmap object
        # its lefternmost corner corresponds to lonmin and latmax
        # it is therefore not exactly the map we want
        bmap = Basemap(llcrnrlon=self.lonmin, llcrnrlat=self.latmin,
                       urcrnrlon=self.lonmax, urcrnrlat=self.latmax,
                       lon_0=lon0, lat_1=lat1, lat_2=lat2, projection='lcc',
                       resolution="c")

        # determination of the MAP coordinates of points
        # in which we are interested in. They are extracted from
        # the geographical coordinates
        if self.hemisphere == 'NH':
            xcoordmin, ycoord = bmap(self.lonmin, self.latmin)
            xcoordmax, ycoord = bmap(self.lonmax, self.latmin)
            xcoord, ycoordmin = bmap(lon0, self.latmin)
            xcoord, ycoordmax = bmap(self.lonmin, self.latmax)
        else:
            xcoordmin, ycoord = bmap(self.lonmin, self.latmax)
            xcoordmax, ycoord = bmap(self.lonmax, self.latmax)
            xcoord, ycoordmax = bmap(lon0, self.latmax)
            xcoord, ycoordmin = bmap(self.lonmin, self.latmin)

        # removing old cordinates
        del ycoord, xcoord

        # we recover the longitude/latitudes of the points that will
        # define the NEW corners of the map
        lonmin, latmin = bmap(xcoordmin, ycoordmin, inverse=True)
        lonmax, latmax = bmap(xcoordmax, ycoordmax, inverse=True)

        # we now create the new colormap
        self.bmap = Basemap(llcrnrlon=lonmin, llcrnrlat=latmin,
                            urcrnrlon=lonmax, urcrnrlat=latmax,
                            lon_0=lon0, lat_1=lat1, lat_2=lat2,
                            projection='lcc', resolution=resolution)

    def make_mask(self, mask_res=300, zorder=1000,
                  edcol='k', bgcolor='white', **kwargs):

        """
        Masks the Lambert Conic projection. 
        
        It creates a dummy NxN array, 
        which is then overlayed to the map using 
        the `imshow` function. Data inside the data 
        are plotted in transparent, 
        data outside the domains are plotted in white.
        Finally, the domain boundaries are drawn.

        :param int mask_res: size of the 2-D array used to draw the mask
        :param int zorder: the plot order at which the mask will be drawn
        :param str edcol: color of the domain boundaries
        :param str bgcolor: color of the mask
        :param dict \**kwargs: keyword argument for the boundary lines

        """
        
        xmin = self.bmap.llcrnrx
        xmax = self.bmap.urcrnrx
        ymin = self.bmap.llcrnry
        ymax = self.bmap.urcrnry

        mask_x = np.linspace(xmin, xmax, mask_res)
        mask_y = np.linspace(ymin, ymax, mask_res)
        mask_x, mask_y = np.meshgrid(mask_x, mask_y)

        lonf, latf = self.bmap(mask_x, mask_y, inverse=True)
        maskarr = np.zeros(lonf.shape)
        maskarr[(lonf <= self.lonmax) &
                (lonf >= self.lonmin) &
                (latf <= self.latmax) &
                (latf >= self.latmin)] = 1
        maskarr = np.ma.masked_where(maskarr == 1, maskarr)

        # We create a Dummy colorbar, containing on the left the bg color
        # and on the right blue colors (which will be unused)
        cmap = ListedColormap([bgcolor, 'blue'])

        # we specify that Masked values will be transparent
        cmap.set_bad(color='blue', alpha=0)

        # We define the normalisation of the imshow
        bounds = [0, 0.5, 1]
        norm = BoundaryNorm(bounds, cmap.N)

        # we draw the imshow
        self.bmap.pcolormesh(mask_x, mask_y, maskarr, 
                             cmap=cmap, zorder=zorder, norm=norm)

        lontt = np.linspace(self.lonmin, self.lonmax, mask_res)
        lattt = np.linspace(self.latmin, self.latmax, mask_res)
        lontt, lattt = np.meshgrid(lontt, lattt)
        xtt, ytt = self.bmap(lontt, lattt)

        self.bmap.plot(xtt[-1, :], ytt[-1, :], color=edcol,
                       zorder=zorder+1, **kwargs)
        self.bmap.plot(xtt[0, :], ytt[0, :], color=edcol,
                       zorder=zorder+1, **kwargs)
        self.bmap.plot(xtt[:, -1], ytt[:, -1], color=edcol,
                       zorder=zorder+1, **kwargs)
        self.bmap.plot(xtt[:, 0], ytt[:, 0], color=edcol,
                       zorder=zorder+1, **kwargs)

        for item in plt.gca().spines.keys():
            plt.gca().spines[item].set_color('w')

    def add_lc_labels(self, labelleft=True, labelright=True, spacelat=10,
                      spacelon=10, zorder=10001, nbspaces=15, **kwargs):

        """
        Add the longitude and latitude labels on masked projection.
        
        Adapted from `NCL <http://www.ncl.ucar.edu/Applications/Scripts/mptick_10.ncl>`_

        .. note:: Does not work if `usetex=True`. 
                  Hence a call to this function sets it to False

        :param bool labelleft: defines whether latitudes labels on
         the left of the map should be displayed
        :param bool labelright: defines whether latitudes on
         the right of the map should be displayed
        :param float spacelat: Latitude spacing between the labels
        :param float spacelon: Longitude spacing between the labels
        :param int zorder: plot order when the labels are drawn
        :param int nbspaces: number of whitespaces to
         add at the end (for left labels) or at
         the beginning (for right labels) of the string
        :param dict \**kwargs: Text keyword arguments (fontsize for instance)
        :return: None

        """

        mpl.rcParams['text.usetex'] = False

        self._ticklon(spacelon, zorder, **kwargs)
        self._ticklat(spacelat, labelleft,
                      labelright, zorder, nbspaces, **kwargs)

    def _ticklon(self, spacelon, zorder, **kwargs):

        """ Function for the ticking of longitude """

        rad_to_deg = 180./np.pi
        deg = r'$^\circ$'

        # Ticking with the longitudes
        lonvalues = np.arange(int(self.lonmin+spacelon),
                              int(self.lonmax), spacelon)
        nlon = len(lonvalues)

        if self.hemisphere == 'NH':
            latval = self.latmin
            stbef = '\n\n'
            staft = ''
        else:
            latval = self.latmax
            stbef = ''
            staft = '\n\n'

        for indlon in range(0, nlon):
            lontt = lonvalues[indlon]
            if self.hemisphere == 'NH':
                lon1_ndc, lat1_ndc = self.bmap(lontt-0.25, latval)
                lon2_ndc, lat2_ndc = self.bmap(lontt+0.25, latval)
            else:
                lon1_ndc, lat1_ndc = self.bmap(lontt+0.25, latval)
                lon2_ndc, lat2_ndc = self.bmap(lontt-0.25, latval)

            slope_bot = (lat1_ndc-lat2_ndc)/(lon1_ndc-lon2_ndc)
            angle = np.arctan(slope_bot)*rad_to_deg

            lon_label_bot = str(np.abs(lontt)) + deg

            if lontt < 0:
                lon_label_bot = stbef+lon_label_bot + "W"+staft
                xtt, ytt = self.bmap(lontt, latval)
                plt.text(xtt, ytt, lon_label_bot, rotation=angle,
                         zorder=zorder, va='center', ha='center',
                         **kwargs)

            elif lontt > 0:
                lon_label_bot = stbef+lon_label_bot + "E"+staft
                xtt, ytt = self.bmap(lontt, latval)
                plt.text(xtt, ytt, lon_label_bot, rotation=angle,
                         zorder=zorder, va='center',
                         ha='center', **kwargs)

            elif lontt == 0:
                lon_label_bot = stbef+lon_label_bot+staft
                xtt, ytt = self.bmap(lontt, latval)
                plt.text(xtt, ytt, lon_label_bot, rotation=angle,
                         zorder=zorder, va='center',
                         ha='center', **kwargs)

    def _ticklat(self, spacelat, labelleft,
                 labelright, zorder, nbspaces, **kwargs):

        """ Add the latitude labels """

        rad_to_deg = 180./np.pi
        space = nbspaces*r' '
        deg = r'$^\circ$'

        # Ticking with the latitudes
        latvalues = np.arange(int(self.latmin),
                              int(self.latmax)+spacelat, spacelat)
        nlat = len(latvalues)

        lon1_ndc, lat1_ndc = self.bmap(self.lonmin, self.latmin)
        lon2_ndc, lat2_ndc = self.bmap(self.lonmin, self.latmax)
        slope_lft = (lat2_ndc-lat1_ndc)/(lon2_ndc-lon1_ndc)

        lon1_ndc, lat1_ndc = self.bmap(self.lonmax, self.latmin)
        lon2_ndc, lat2_ndc = self.bmap(self.lonmax, self.latmax)
        slope_rgt = (lat2_ndc-lat1_ndc)/(lon2_ndc-lon1_ndc)

        if self.hemisphere == "NH":
            rotate_left = -90
            rotate_right = 90
        else:
            rotate_left = 90
            rotate_right = -90

        for indlat in range(0, nlat):

            lattt = latvalues[indlat]
            lat_label_rgt = space + str(np.abs(lattt)) + deg

            if lattt < 0:
                lat_label_lft = str(np.abs(lattt)) + deg+"S"+space
                lat_label_rgt = lat_label_rgt + "S"

            elif lattt > 0:
                lat_label_lft = str(np.abs(lattt)) + deg+"N"+space
                lat_label_rgt = lat_label_rgt + "N"

            elif lattt == 0:
                lat_label_lft = str(np.abs(lattt)) + deg + space

            angle = rad_to_deg * np.arctan(slope_lft) + rotate_left
            xtt, ytt = self.bmap(self.lonmin, lattt)
            if labelleft:
                plt.text(xtt, ytt, lat_label_lft, rotation=angle,
                         zorder=zorder, va='center',
                         ha='center', **kwargs)

            angle = rad_to_deg * np.arctan(slope_rgt) + rotate_right
            xtt, ytt = self.bmap(self.lonmax, lattt)
            if labelright:
                plt.text(xtt, ytt, lat_label_rgt, rotation=angle,
                         zorder=zorder, va='center',
                         ha='center', **kwargs)
