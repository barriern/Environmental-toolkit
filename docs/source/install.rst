
.. _install:

Getting the *envtoolkit* library
====================================

.. _install-req:

Module requierements
***********************************

The *envtoolkit* library requieres the following packages:
    - :py:mod:`numpy` (for maths)
    - :py:mod:`matplotlib.pyplot` (for plots)
    - :py:mod:`matplotlib.colors` (for colors)
    - :py:mod:`matplotlib.path` (for polygons)
    - :py:mod:`netCDF4` (http://unidata.github.io/netcdf4-python/)
    - :py:mod:`spectrum` (http://pythonhosted.org/spectrum/)

.. _install-down:

Downloading the source
*******************************************

The *envtoolkit* library is available on a  `GitHub repository <https://github.com/barriern/Environmental-toolkit>`_. It can be downladed by using `Git <https://git-scm.com/>`_ as follows:

.. code-block:: none
    
    > git clone https://github.com/barriern/Environmental-toolkit

The last version of the library is obtained as follows:

.. code-block:: none
    
    > cd Environmental-toolkit
    > git pull

.. _install-inst:

Module installation
*********************************

The module installation is achieved as follows:

.. code-block:: bash

    > cd Environmental-toolkit
    > python setup.py install

.. _install-doc:

Generating the documentation
*******************************************

The *envtoolkit* documentation can be built offline by using the `Sphinx <http://www.sphinx-doc.org/>`_ library. This is achieved as follows:

.. code-block:: bash

    > cd docs
    > make html
    > make latexpdf

The *HTML* documentation is built in the *docs/build/html* directory (*index.html* file), while the *LaTex* documentation is built in the *docs/build/latex* directory (*envtoolkit.pdf* file).
