import os
from setuptools import setup, find_packages

__docformat__ = "restructuredtext en"

with open("README.md", "r") as fh:
    long_description = fh.read()

VERSION_FILE = 'VERSION'
with open(VERSION_FILE) as fv:
    version = fv.read().strip()

setup(
    name="envtoolkit",
    version=version,
    author="Nicolas Barrier",
    author_email="nicolas.barrier@ird.fr",
    description=("Python package dedicated to the Earth Scientists."),
    long_description_content_type="text/markdown",
    packages=find_packages(),
    url="https://github.com/barriern/Environmental-toolkit",
    requires=['numpy',
              'netCDF4',
              'matplotlib',
              're',
              'datetime.datetime',
              'spectrum'],

    long_description=long_description,

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

           
