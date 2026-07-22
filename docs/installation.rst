Installation
============

The Python package installs with pip, but the calculation **backends** (DFTB+,
``modes``, ``xtb``, ``tblite``) are compiled programs distributed via
conda-forge and cannot be installed with pip.

Package only
------------

.. code-block:: bash

   python -m pip install thermoscreening

Full environment (recommended)
------------------------------

Create an environment containing every backend, then install the released
Python package:

.. code-block:: bash

   conda create -n thermoscreening -c conda-forge \
       python=3.12 dftbplus xtb tblite-python pip
   conda activate thermoscreening
   python -m pip install thermoscreening

For development from a checkout, ``conda env create -f environment.yml``
installs the same backends and an editable package with the test tools.

After installation, ``thermo doctor`` reports every backend and lists the usable
engines. Use ``thermo doctor --engine dftb+`` (or ``xtb``/``xtb-cli``) when a
specific backend is required.

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

The DFTB+ engine needs Slater-Koster (``.skf``) parameter files.
ThermoScreening downloads supported sets to an automatically discovered
user-local directory:

.. code-block:: bash

   thermo setup-dftb --parameter-set 3ob      # or: --parameter-set mio
   thermo setup-dftb --solvent water          # a GBSA solvent parameter file

Use ``DFTB_PREFIX`` only to override that directory with a custom installation.
``thermo doctor`` verifies that executables can start, checks the selected
parameter set, and reports the optional xTB backends.
