API reference
=============

High-level pipeline
-------------------

.. autofunction:: ThermoScreening.thermo.api.dftbplus_thermo
.. autofunction:: ThermoScreening.thermo.api.xtb_thermo
.. autofunction:: ThermoScreening.thermo.api.xtb_cli_thermo
.. autofunction:: ThermoScreening.thermo.api.xtb_fukui_indices
.. autofunction:: ThermoScreening.thermo.api.orca_thermo
.. autofunction:: ThermoScreening.thermo.api.cclib_thermo
.. autofunction:: ThermoScreening.thermo.api.pyscf_thermo
.. autofunction:: ThermoScreening.thermo.api.run_thermo
.. autofunction:: ThermoScreening.thermo.api.execute

Screening
---------

.. autofunction:: ThermoScreening.thermo.screening.screen
.. autofunction:: ThermoScreening.thermo.screening.rank_by_gibbs
.. autoclass:: ThermoScreening.thermo.screening.ScreeningJob
   :members:

Conformer generation
--------------------

.. autofunction:: ThermoScreening.thermo.conformers.generate
.. autofunction:: ThermoScreening.thermo.conformers.generate_thermo_ensemble
.. autofunction:: ThermoScreening.thermo.conformers.write_conformers

Reactions and redox
-------------------

.. autofunction:: ThermoScreening.thermo.reactions.reaction_free_energy
.. autofunction:: ThermoScreening.thermo.reactions.reduction_potential
.. autofunction:: ThermoScreening.thermo.reactions.calibrate_reduction_reference

Dataset redox screening
-----------------------

.. autoclass:: ThermoScreening.thermo.redox_screening.DeltaRedoxModel
   :members: fit, predict, correction_coefficients
.. autofunction:: ThermoScreening.thermo.redox_screening.canonical_structure_identity
.. autofunction:: ThermoScreening.thermo.redox_screening.audit_redox_dataset
.. autofunction:: ThermoScreening.thermo.redox_screening.balanced_group_folds
.. autofunction:: ThermoScreening.thermo.redox_screening.grouped_delta_validation
.. autofunction:: ThermoScreening.thermo.redox_screening.pareto_front

Acid dissociation (pKa)
------------------------

.. autofunction:: ThermoScreening.thermo.pka.pKa
.. autofunction:: ThermoScreening.thermo.pka.calibrate_proton_reference

Transition states and kinetics
-------------------------------

.. autofunction:: ThermoScreening.thermo.kinetics.eyring_rate_constant
.. autofunction:: ThermoScreening.thermo.kinetics.wigner_tunneling_correction

Conformer ensembles
-------------------

.. autofunction:: ThermoScreening.thermo.ensemble.boltzmann_weights
.. autofunction:: ThermoScreening.thermo.ensemble.ensemble_free_energy
.. autofunction:: ThermoScreening.thermo.ensemble.lowest_gibbs
.. autoclass:: ThermoScreening.thermo.ensemble.EnsembleThermo
   :members:

Thermochemistry core
--------------------

.. autoclass:: ThermoScreening.thermo.system.System

.. autoclass:: ThermoScreening.thermo.thermo.Thermo
   :members: total_energy, total_enthalpy, total_gibbs_free_energy,
             total_entropy, total_heat_capacity, total_EeGtot, electronic_energy,
             temperature_scan, imaginary_mode_wavenumber

Coordinate and frequency readers
--------------------------------

.. autofunction:: ThermoScreening.thermo.api.read_xyz
.. autofunction:: ThermoScreening.thermo.api.read_gen
.. autofunction:: ThermoScreening.thermo.api.read_coord
.. autofunction:: ThermoScreening.thermo.api.read_vibrational
.. autofunction:: ThermoScreening.calculator.orca.read_orca_hess
.. autofunction:: ThermoScreening.calculator.qm.read_cclib

Backend setup
-------------

.. autofunction:: ThermoScreening.cli.dftb_setup.install_slakos
.. autofunction:: ThermoScreening.cli.dftb_setup.install_gbsa_param
.. autofunction:: ThermoScreening.cli.dftb_setup.check_dftb_setup
