import os
from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="envtoolkit",
    version="1.0.1",
    author="Nicolas Barrier",
    author_email="nicolas.barrier@ird.fr",
    description=("Python package dedicated to the Earth Scientists. It contains a module for the processing of time-series (filtering, anomalies), a module for spatial analysis, a module with color-specific functions, a module for NetCDF processing and a module for spectral analysis"),
    long_description_content_type="text/markdown",
    packages=find_packages(),
    url="https://github.com/barriern/Environmental-toolkit",
    requires=['numpy',
              'netCDF4',
              'matplotlib',
              're',
              'mpl_toolkits.basemap',
              'datetime.datetime',
              'spectrum'],

    classifiers = [
        #"Development Status :: 5 - Production/Stable",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Visualization",
        "Topic :: Scientific/Engineering :: Physics",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Unix",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
    ],

)

           
