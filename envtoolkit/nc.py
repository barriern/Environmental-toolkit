
""" Module that handles functions relative to NetCDF files """


def copy_ncstruct(fin, fout):

    """ Copy the internal structure of a NetCDF file. 

    It copies the dimension names, variable names and attributes
    from a file to another file. 

    :param netCDF4.Dataset fin: Source file
    :param netCDF4.Dataset fout: Destination file
    
    .. note :: The file contents are not copied! Only the structure of the file

    .. code-block :: python 

        from netCDF4 import Dataset
        import envtoolkit.nc

        fin = Dataset("source_file.nc", "r")
        fout = Dataset("dest_file.nc", "w")

        # copy the dimensions/variables/attributes of
        # the source file into the dest file
        envtoolkit.nc.copy_ncstruct(fin, fout)
        
        # Fill in the variable of the destination file

        # closing the files
        fin.close()
        fout.close()

    :param netCDF4.Dataset fin: the source file

    :param netCDF4.Dataset fout: the destination file
    """

    # Working on the dimensions
    dim = fin.dimensions
    for dimname, dimval in zip(dim.keys(), dim.values()):
        if dimval.isunlimited():
            # If record dimension, then we also set
            # the output dimension to unlimited
            fout.createDimension(dimname, None)
        else:
            # len(v) is the way to recover the value of the dimension
            # cf. netCDF4 documentation
            fout.createDimension(dimname, len(dimval))

    # Working on the variables
    var = fin.variables
    for varname, varval in zip(var.keys(), var.values()):
        # we create a variable of the same
        # type and dimensions as in the input file
        fout.createVariable(varname, varval.datatype, varval.dimensions)
        vout = fout.variables[varname]

        # we loop over the variable attributes
        # varval.ncattrs -> recover the attribute names ("units", etc)
        for attr in varval.ncattrs():
            vout.setncattr(attr, getattr(varval, attr))

    # We do the same thing for the global (file) attributes
    for attr in fin.ncattrs():
        fout.setncattr(attr, getattr(fin, attr))


def extract_date(fin, timevar_name="time", units=None, calendar='gregorian', timeoff=0):

    """ Converts the time array of a NetCDF file into a date.

    :param netCDF4.Dataset fin: the NetCDF file
    :param str timevar_name: the name of the time variable
    :param str units: the time units (used only if no time units in the
     the file)
    :param str calendar: the time calendar (used only if no time units in the 
     the file)
    :param float timeoff: A time offset that is *added* to the time array

    :return: an array of :py:class:`datetime.datetime` object
    :rtype: numpy.array

    """

    from netcdftime import utime

    try:
        # try to extract the time variable from the file
        timevar = fin.variables[timevar_name]
    except:
        raise IOError("The %s variable does not exist in the input file" % timevar_name)

    # check the existence of the units attribute
    if "units" in timevar.ncattrs():
        units = getattr(timevar, 'units')
    else:
        # if no time units, we check that the unit argument is set, and that it is in right format
        if units is None:
            raise IOError("No time units provided in the file. Provide the units as an argument")
        else:
            if not isinstance(units, str):
                raise ValueError("The units argument must be a string")
    
    # check the existence of the calendar attribute
    if "calendar" in timevar.ncattrs():
        calendar = getattr(timevar, 'calendar')
    else:
        # if no time calendar, we check that the unit argument is set, and that it is in right format
        if calendar is None:
            print('No calendar provided. Gregorian calendar is used')
            calendar = 'gregorian'
        else:
            if not isinstance(calendar, str):
                raise ValueError("The calendar argument must be a string")

    # creates the utime object for conversion
    cdftime = utime(units, calendar=calendar)

    # conversion of numerics into date
    date = cdftime.num2date(timevar[:] + timeoff)

    return date
