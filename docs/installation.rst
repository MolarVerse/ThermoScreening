Installation
============

The Python package installs with pip, but the calculation **backends** (DFTB+,
``modes``, ``xtb``, ``tblite``) are compiled programs distributed via
conda-forge and cannot be installed with pip.

Package only
------------

.. code-block:: bash

   python -m pip install .
   # for development and tests
   python -m pip install -e ".[test,lint]"

Full environment (recommended)
------------------------------

The bundled Conda environment installs the package **and** every backend:

.. code-block:: bash

   conda env create -f environment.yml
   conda activate thermoscreening

After this, ``thermo doctor`` should report every backend as found.

Backends
--------

.. list-table::
   :header-rows: 1
   :widths: 20 35 45

   * - Backend
     - Needed for
     - Install
   * - ``dftb+``, ``modes``
     - the DFTB+ engine
     - ``conda install -c conda-forge dftbplus``
   * - ``xtb``
     - ``--engine xtb-cli`` (native xtb)
     - ``conda install -c conda-forge xtb``
   * - ``tblite``
     - ``--engine xtb`` (in-process GFN-xTB)
     - ``conda install -c conda-forge tblite-python``

DFTB+ and Slater-Koster setup
-----------------------------

The DFTB+ engine needs Slater-Koster (``.skf``) parameter files, located via the
``DFTB_PREFIX`` environment variable. ThermoScreening can download the supported
sets and the GBSA solvation parameters for you:

.. code-block:: bash

   thermo setup-dftb --parameter-set 3ob      # or: --parameter-set mio
   thermo setup-dftb --solvent water          # a GBSA solvent parameter file
   export DFTB_PREFIX="$HOME/.local/share/thermoscreening/slakos/3ob-3-1/"

``thermo doctor`` checks the executables, ``DFTB_PREFIX``, and the optional xTB
backends, and prints the exact install command for anything missing.
