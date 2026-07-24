QM9 thermochemistry benchmark
=============================

The `QM9 dataset <https://doi.org/10.6084/m9.figshare.978904_D12>`_
contains optimized geometries, harmonic frequencies, and thermochemical
properties for 133,885 neutral CHONF molecules. The associated
`Scientific Data paper <https://doi.org/10.1038/sdata.2014.22>`_ reports that
the calculations used B3LYP/6-31G(2df,p) and provides ``U0``, ``U``, ``H``,
``G``, and ``Cv`` at 298.15 K.

Run the benchmark from a source checkout:

.. code-block:: bash

   python scripts/validate_qm9.py

The command downloads the original 86 MB archive, verifies its MD5 checksum,
and streams a pinned 1,000-molecule sample without extracting the full dataset.
The sample contains the first 100 molecules, 450 evenly spaced records, and 450
seeded random records.

Method
------

QM9 stores the zero-point energy and internal energy at 0 K, so the electronic
energy is reconstructed as ``U0 - ZPVE``. ThermoScreening then recomputes all
thermal quantities from the published geometry and frequencies.

The XYZ records do not store the rotational symmetry number used by Gaussian.
The pinned source settings use ``sigma=2`` for records 3, 4, and 23 and
``sigma=1`` for the other sampled records. This avoids treating a symmetry
metadata difference as a thermochemistry error.

Reference result
----------------

Energy differences are reported in microhartree. The source energies are
printed to six decimal places and ``Cv`` to three decimal places.

.. list-table::
   :header-rows: 1

   * - Property
     - Mean absolute error
     - Maximum absolute error
   * - ``U0``
     - 0.246147 microhartree
     - 0.497836 microhartree
   * - ``U``
     - 0.401046 microhartree
     - 1.336198 microhartree
   * - ``H``
     - 0.404191 microhartree
     - 1.382253 microhartree
   * - ``G``
     - 1.622078 microhartree
     - 3.211252 microhartree
   * - ``Cv``
     - 0.000259 cal/(mol K)
     - 0.000585 cal/(mol K)

The maximum Gibbs-energy difference is approximately 0.0020 kcal/mol. The
benchmark therefore reproduces the published quantities to the precision
available from QM9's rounded coordinates, frequencies, and energies.

Near-linear geometry
--------------------

QM9 record 25 is a slightly bent cyanogen geometry with six vibrational modes.
Its smallest-to-largest principal-moment ratio is about ``9.7e-8``. A looser
linearity tolerance treated it as exactly linear and expected seven modes.
The external benchmark identified this edge case; the regression test now
requires the geometry to be handled as nonlinear, matching the QM9
thermochemistry.

Interpretation
--------------

This is an independent numerical validation of the ideal-gas rigid-rotor
harmonic-oscillator implementation over varied molecular sizes and structures.
It does not measure the predictive accuracy of B3LYP or compare computed
thermochemistry with experiment.
