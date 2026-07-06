API reference
=============

High-level pipeline
-------------------

.. autofunction:: ThermoScreening.thermo.api.dftbplus_thermo
.. autofunction:: ThermoScreening.thermo.api.xtb_thermo
.. autofunction:: ThermoScreening.thermo.api.xtb_cli_thermo
.. autofunction:: ThermoScreening.thermo.api.run_thermo
.. autofunction:: ThermoScreening.thermo.api.execute

Screening
---------

.. autofunction:: ThermoScreening.thermo.screening.screen
.. autoclass:: ThermoScreening.thermo.screening.ScreeningJob
   :members:

Conformer generation
--------------------

.. autofunction:: ThermoScreening.thermo.conformers.generate
.. autofunction:: ThermoScreening.thermo.conformers.write_conformers

Reactions and redox
-------------------

.. autofunction:: ThermoScreening.thermo.reactions.reaction_free_energy
.. autofunction:: ThermoScreening.thermo.reactions.reduction_potential

Thermochemistry core
--------------------

.. autoclass:: ThermoScreening.thermo.system.System

.. autoclass:: ThermoScreening.thermo.thermo.Thermo
   :members: total_energy, total_enthalpy, total_gibbs_free_energy,
             total_entropy, total_heat_capacity, total_EeGtot, electronic_energy

Coordinate and frequency readers
--------------------------------

.. autofunction:: ThermoScreening.thermo.api.read_xyz
.. autofunction:: ThermoScreening.thermo.api.read_gen
.. autofunction:: ThermoScreening.thermo.api.read_coord
.. autofunction:: ThermoScreening.thermo.api.read_vibrational

Backend setup
-------------

.. autofunction:: ThermoScreening.cli.dftb_setup.install_slakos
.. autofunction:: ThermoScreening.cli.dftb_setup.install_gbsa_param
.. autofunction:: ThermoScreening.cli.dftb_setup.check_dftb_setup
