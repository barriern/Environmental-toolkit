import os
from setuptools import setup, find_packages

setup(
    name="envtoolkit",
    version="1.0.0",
    author="Nicolas Barrier",
    author_email="nicolas.barrier@ird.fr",
    description=("Python package dedicated to the Earth Scientists. It contains a module for the processing of time-series (filtering, anomalies), a module for spatial analysis, a module with color-specific functions, a module for NetCDF processing and a module for spectral analysis"),
    license=None,  
    packages=find_packages(),
    requires=['numpy',
              'netCDF4',
              'matplotlib',
              're',
              'mpl_toolkits.basemap',
              'datetime.datetime',
              'spectrum'],
    )
           
