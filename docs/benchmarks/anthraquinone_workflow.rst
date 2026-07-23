Anthraquinone reference benchmark
=================================

ThermoScreening can reproduce the DFTB/3ob/COSMO redox-potential tables from
Kröll *et al.*, `Phys. Chem. Chem. Phys. 24, 2022
<https://doi.org/10.1039/D2CP01717B>`_. The raw calculations and processed
tables are published in the immutable `Zenodo record 20796838
<https://doi.org/10.5281/zenodo.20796838>`_.

Run the benchmark from a source checkout:

.. code-block:: bash

   python scripts/validate_anthraquinone_zenodo.py

The command:

- downloads the raw DFTB/3ob/COSMO-DMF outputs and published plot tables;
- verifies both archives against pinned SHA-256 checksums;
- processes 13 molecules in three charge states with :func:`run_thermo`;
- uses the archived analysis settings, including ``spin=1`` for monoanions and
  ``symmetry_number=1`` because the historical workflow disabled symmetry
  analysis;
- calibrates the parent anthraquinone result against each published table; and
- compares 13 first-reduction and 12 overall two-electron potentials.

Reference result
----------------

.. list-table::
   :header-rows: 1

   * - Metric
     - Result
     - Acceptance limit
   * - Mean absolute error
     - 0.021320 mV
     - 0.05 mV
   * - Maximum absolute error
     - 0.334662 mV
     - 0.5 mV

The largest difference is the first reduction of molecule 7
(1,4-diaminoanthraquinone). The remaining discrepancy is below the precision
relevant to the published millivolt-scale table.

Interpretation
--------------

This benchmark verifies that ThermoScreening reproduces the archived
thermochemistry and redox post-processing when given the same raw electronic
energies, frequencies, geometries, spin convention, symmetry number, and
reference calibration. It does not establish experimental predictive accuracy.

The archived ``DFTB_2e`` column is the calibrated neutral-to-dianion
two-electron average. It is reproduced with ``n_electrons=2`` and must not be
interpreted as the stepwise anion-to-dianion ``E2`` reported by the current
three-state redox workflow. The current workflow calibrates ``E1`` and ``E2``
separately and reports ``E2e`` as their arithmetic mean.

The ordinary test suite also retains a compact synthetic anthraquinone fixture.
That fixture checks workflow invariants without network access; the Zenodo
command is the quantitative external-data benchmark.
