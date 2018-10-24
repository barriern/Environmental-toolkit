
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
    - :py:mod:`mpl_toolkits.basemap` (for maps)
    - :py:mod:`netCDF4` (`documentation here <http://unidata.github.io/netcdf4-python/>`_)
    - :py:mod:`spectrum` (`documentation here <http://pythonhosted.org/spectrum/>`_)

.. _install-down:

Downloading the source
*******************************************

The *envtoolkit* library is available via the `Sourcesup platform <https://sourcesup.renater.fr/>`_. It can be downladed by using `Subversion <https://subversion.apache.org>`_ as follows:

.. code-block:: none
    
    > mkdir envtoolkit
    > cd envtoolkit
    > svn checkout https://subversion.renater.fr/envtoolkit/trunk

The last version of the library is obtained as follows:

.. code-block:: none
    
    > cd envtoolkit
    > svn update

.. _install-inst:

Module installation
*********************************

The module installation is achieved as follows:

.. code-block:: bash

    > cd trunk
    > python setup.py install --home=/install/directory

The *home* option indicates the location where the library will be installed. Note that the install directory must be in the *PYTHONPATH*. If it is not the case, you should edit your *.bashrc* file and add the following line:

.. code-block:: bash

    export PYTHONPATH=$PYTHONPATH:/install/directory/lib/python  # for BASH/KSH

.. code-block:: csh

    setenv PYTHONPATH ${PYTHONPATH}:/install/directory/lib/python # for CSH/TCSH

.. _install-doc:

Generating the documentation
*******************************************

The *envtoolkit* documentation can be built offline by using the `Sphinx <http://www.sphinx-doc.org/>`_ library. This is achieved as follows:

.. code-block:: bash

    > cd trunk/docs
    > make html
    > make latexpdf

Note that the environment variable *SPHINXBUILD*, containing the name of the Sphinx executable file, must be defined. 

.. code-block:: bash

    export SPHINXBUILD=sphinx-build-2.7  # for BASH/KSH

.. code-block:: csh

    setenv SPHINXBUILD sphinx-build-2.7  # for CSH/TCSH

The *HTML* documentation is built in the *trunk/docs/build/html* directory (*index.html* file), while the *LaTex* documentation is built in the *trunk/docs/build/latex* directory (*envtoolkit.pdf* file).

.. note::

    An HTML version of this documentation is `available here <http://mio.pytheas.univ-amu.fr/~barrier.n/>`_ (see the *Python/NCL* tab).
