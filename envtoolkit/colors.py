
"""
Functions/classes relative to colors.
"""

import numpy as np
import os
import re

def make_percentile_cmap(array, perc):

    """
    Returns colormaps limits using percentiles. The
    minimum color value. 

    :param numpy.array array: Input array

    :param int perc : Percentile (in percentage)

    :return: A tuple containing the lower 
     and upper colorbar limits
     associated with the array percentiles.

    :rtype: tuple

    """

    # converts from ND to 1D
    array = np.ravel(array)

    # mask data where missing and extracts only unmasked
    # data
    array = np.ma.masked_where(array!=array, array)
    array = array[np.ma.getmaskarray(array)==False]

    if len(array) == 0:
        cmin, cmax = None, None

    else:
        # computes the percentiles
        cmin = np.percentile(array, perc)
        cmax = np.percentile(array, 100 - perc)

    # rounds the percentiles (unnecessary???)
    #cmin = unit_round(cmin, "floor")
    #cmax = unit_round(cmax, "ceil")

    return cmin, cmax


class FancyCmap(object):

    """
    Class for generating fancy colormaps.

    It creates the attributes `offset` 
    (containing the linear mapping segments) and 
    `col` (containing the associated colors). 

    :param str|int|numpy.array input_cmap: Colormap to use. Can be 
     a string (the name of the .svg file to read), a list 
     of named colors or a 3xN/Nx3 numpy.array (percentage of RGB values)

    :param bool smooth: If the `input_cmap` argument is a list 
     or an array, define whether the colormap 
     is smoothed (i.e gradient type colormap) or not.

    :param bool reverse: Defines whether the colors should be
     reversed.

    :param list list_offset: List of values containing 
     mapping segments. Only used for user defined colormaps. 
     If None, regular offseting is used.

    """

    def __init__(self, input_cmap, smooth=False, \
                 reverse=False, list_offset=None):

        """
        Initialisation of the :py:class:`FancyCmap` class.
        """


        self.dict_cmap = None
        self.smooth = smooth
        self.reverse = reverse
        self.list_offset = list_offset
        self.input_cmap = input_cmap

        # if cmap is a string, i.e a .svg file
        if isinstance(input_cmap, str):
            self._initfromstr()

        # if input_cmap is a list
        elif isinstance(input_cmap, list):
            self._initfromlist()

        # if input_cmap is a numpy array
        elif isinstance(input_cmap, np.ndarray):
            self._initfromarray()

        # if input_cmap is none of this
        else:
            raise ValueError('The argument of the variables must either \
            be a list, a numpy array or a string, not %s' % isinstance(input_cmap))

        # If the colormap is to be reversed, the offset is changed
        if reverse:
            self.offset = 1 - self.offset[::-1]
            self.col = self.col[::-1, :]

        # creation of the mapping directory
        self._makesegment()

    def makecmap(self, nbcol=256):

        """
        Create the colormap instance. 
        
        The output of this function will be provided in `cmap` argument
        of the :py:func:`matplotlib.pyplot.contourf` or similar functions.

        :param int nbcol: Number of colors 
         in the colormap 

        :param bool reverse: defines whether the colormap should be 
         reversed or not (equivalent to the _r for instance) 

        :return: A colormap object

        :rtype: matplotlib.colors.LinearSegmentedColormap

        """

        from matplotlib.colors import LinearSegmentedColormap

        # Matplotlib function that creates the colormap from
        # the mapping directory
        cmap = LinearSegmentedColormap('', self.dict_cmap, N=nbcol)

        return cmap

    def generate_cbticks(self, cbar):

        """
        Writes colorbar ticks for discrete colormaps. 

        It writes a label at each edges of the colorbar.

        :param matplotlib.colorbar.Colorbar cbar: The colorbar object

        :return: None

        """

        if self.smooth:
            print "Warning: This function cannot be used with smoothed colormap."
            print "Nothing has been changed"
            return

        cmin, cmax = cbar.get_clim()

        # delta_c is the range of the colormap limit
        delta_c = cmax - cmin

        # the first element of the output list is at the cmin location
        lout = [cmin]

        # now we loop over all the elements of the offset array,
        # and we append the values to the lout list
        for indoff in range(1, len(self.offset)):
            lout.append(cmin + self.offset[indoff]*delta_c)

        # we seek for single elements (in case of discrete cmaps,
        # the list will contain duplications)
        lout = np.unique(np.array(lout))

        return lout

    def _makesegment(self):

        """
        Creates the dictionnary that will be
        used to create the colormap.
        This is a dictionnary similar to those 
        in the _cm.py file. It contains three entries, 'red', 
        'green' and 'blue'. They are list of tuples, each 
        element containing the offset
        and the percentage of the color
        """

        # initialisation of the color arrays
        red = []
        green = []
        blue = []

        # we loop over all the element of the offset array
        # in the lists, we append each colors twice (cf. Matlotlib site)
        for indoff in range(0, len(self.offset)):
            red.append((self.offset[indoff],
                        self.col[indoff, 0],
                        self.col[indoff, 0]))
            green.append((self.offset[indoff],
                          self.col[indoff, 1],
                          self.col[indoff, 1]))
            blue.append((self.offset[indoff],
                         self.col[indoff, 2],
                         self.col[indoff, 2]))

        self.dict_cmap = {'red':   red,
                          'green': green,
                          'blue':  blue}

    def _initfromlist(self):

        """
        Initialises the class by using a list of named colors

        :param list input_cmap: List of named color strings
        :param bool smooth: Logical that defines
         whether the colormap should me smoothed out
         (gradient type) or not (discrete type)

        """

        from matplotlib.colors import colorConverter

        # initialisation of the output list
        # that will contain RGB colormaps
        col = []

        # we loop over all the colors contained in the input_cmap
        for color in self.input_cmap:

            # we try to convert the named color into RGB values
            # if successful, we append the RGB output into the col list
            try:
                rgb = colorConverter.to_rgb(color)
                col.append(rgb)
            except:
                raise ValueError('Colorname %s is not a valid name' % color)

        # we transform the output list into an array
        col = np.array(col)

        # we call the function that creates the offset array from a color array
        self._makecoloffset(col)

    def _makecoloffset(self, col):

        """
        Defines the offset and color 
        attributes.
        
        Used only for user-defined colormaps. 

        :param list col: List of colors

        """


        # Check that if the input list is defined,
        # it has a consistent number of elements
        if self.list_offset is not None:
            if len(self.list_offset) != len(col) + 1:
                raise IOError("The number of elements within \
                        the 'list_offset' list is not good'. Should \
                        be %d" % len(col) + 1)

        # here, we define the offset depending on whether a smoothing of
        # the colormap is applied
        # if smooth=True, the offset will contain elements which are
        # all differents
        # if 4 colors, it will move from c1 to c2 between 0 and 33%,
        # from c2 to c3 between 33 and 66%, from c3 to c4
        # between 66 and 100%
        if self.smooth:

            if self.list_offset is None:
                self.offset = np.linspace(0, 1, len(col[:, 0]))
            else:
                offset = np.array(self.list_offset)
                offset = 0.5*(offset[1:] + offset[:-1])
                offset[0] = 0
                offset[-1] = 1
                self.offset = offset
            self.col = col
            return

        else:

            # if not(smooth), col. elements and off. work as pairs (fs2008.svg)
            # the same color is given two different offsets
            # and then the same offset is given two different colors

            # we create an offset array as if smooth=True
            if self.list_offset is None:
                offset = np.linspace(0, 1, len(col[:, 0]) + 1)
            else:
                offset = self.list_offset

            offint = []
            colint = []

            for indoff in range(0, len(offset)-1):

                # color elements and offsets work as pair
                colint.append(col[indoff, :])
                colint.append(col[indoff, :])
                offint.append(offset[indoff])
                offint.append(offset[indoff+1])

            colint = np.array(colint)
            offint = np.array(offint)
            self.offset = offint
            self.col = colint

            return

    def _initfromstr(self):

        """
        Initialises the class by using
        a string (i.e. name of the .svg file in the $SVG directory)
        """


        # try to look for the SVG environment.
        # if not defined, error is raised
        try:
            dirin = os.environ['SVG']
        except:
            raise EnvironmentError('The SVG environment, in which the .svg \
            files are created, needs to be defined!')

        output = []  # init of the offset file

        try:
            # the file is opened
            filein = open(dirin+self.input_cmap.replace('.svg', '') + '.svg', 'r')
        except:
            raise IOError('filename %s.svg is not defined' % (dirin+self.input_cmap))

        alllines = filein.readlines()  # we read the lines
        filein.close()

        # regular expression which is read
        pattern = r' *<stop offset="([0-9]{1,3}\.[0-9]{2})%" ' + \
                  r'*stop-color="rgb\(([0-9]{1,3}), ?([0-9]{1,3}), ?([0-9]{1,3})\)"'

        valid = re.compile(pattern)

        for line in alllines:

            # if the line l matches the pattern, we recover the variables
            if valid.match(line):

                # we append the variable in the output list
                output.append(list(valid.match(line).groups()))

        output = np.array(output)  # we convert the list into an array

        # we recover the offset percentage
        self.offset = output[:, 0].astype(np.float)/100.

        # we recover the col and convert them into percentage
        self.col = output[:, 1:].astype(np.int)/255.

    def _initfromarray(self):

        """
        Initialises the class by using a numpy array
        of colors (in RGB percentage or in 255 values)

        """


        # if it is an array, we recover it
        if (3 not in self.input_cmap.shape) and (self.input_cmap.ndim == 2):
            raise ValueError('array must eithe be Nx3 or 3xN, not %s' % str(input_cmap.shape))

        else:

            col = np.empty(self.input_cmap.shape)
            col[:] = self.input_cmap[:]

            if col.shape[0] == 3:
                col = np.transpose(col)  # transpose the array if needed

            if col.max() > 1:
                print('Values are greater than 1. We suspect ' +
                        'that you used [0, 255] values. We therfore ' +
                        'have divided by 255')
                col = col/255.

            self._makecoloffset(col)

        return


def make_boundary_norm(ncolors, boundaries):

    """ Normalisation of colormaps for "discrete" mapping. If

    .. code-block:: python

        boundary[i] <= value < boundary[i+1]

    then the data is mapped to color(i) 
    (cf :py:class:`matplotlib.colors.BoundaryNorm` documentation).

    :param int ncolors: Number of colors to display on the colormap
    :param numpy.array boundaries: Array containing the colorbar boundaries (i.e 
     the values of the colorbar edges). *Must have (ncolors+1) elements*.

    .. warning:: When using the output norm in a pcolormesh/imshow function, the
       colormap that is used must have ncolors colors. This is achieved either
       by defining your own colormap with the :py:class:`FancyCmap` class or by using
       the :py:func:`subspan_default_cmap` function on Matplotlib default colormaps.

    :rtype: :py:class:`matplotlib.colors.Normalize`

    :return: A norm instance that contains the discrete mapping

    """

    import matplotlib.colors

    norm = matplotlib.colors.BoundaryNorm(boundaries, ncolors=ncolors)

    if len(boundaries) != ncolors+1:
        raise IOError('The number of elements in the boundaries ' + \
        ' array is not good. %d expected, %d received' % (ncolors+1, len(boundaries)))

    return norm


def subspan_default_cmap(cmapname, nbcol):

    """
        Function that allows to generate a colormap containing *nbcol* colors
        from the default colormap data mapping dictionnaries
        (cf the `datad` dictionnary in the :py:mod:`matplotlib._cm` module)

        :param str cmapname: Name of the matplotlib colormap
        :param int nbcol: Number of colors to extract from the *cmapname* colormap
        :return: A colormap object
        :rtype: matplotlib.colors.LinearSegmentedColormap

    """

    from matplotlib._cm import datad
    from matplotlib.colors import LinearSegmentedColormap

    if not(isinstance(nbcol, int)):
        raise IOError("The nbcol must be an integer")

    if not(isinstance(cmapname, str)):
        raise IOError("The 'cmapname' must be a string")

    try:
        dictout = datad[cmapname]
    except:
        raise IOError("The %s colormap does not exist in Matplotlib" % cmapname)

    # Matplotlib function that creates the colormap from
    # the mapping directory
    cmap = LinearSegmentedColormap('', dictout, N=nbcol)

    return cmap
