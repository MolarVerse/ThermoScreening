Usage
=====

Command line
------------

The ``thermo`` command has four subcommands.

``thermo doctor``
    Report which backends (dftb+, modes, DFTB_PREFIX + a Slater-Koster file, and
    the optional xtb/tblite) are available.

``thermo setup-dftb``
    Download a Slater-Koster parameter set (``--parameter-set 3ob|mio``) or a
    GBSA solvent parameter file (``--solvent <name>``).

``thermo conformers``
    Generate conformers from a SMILES (RDKit ETKDG) and write them as ``.xyz``
    files, ready to screen:

    .. code-block:: bash

       thermo conformers "CCCC" --out-dir butane --max-conformers 10
       thermo screen butane --engine xtb-cli

``thermo screen``
    Run the thermochemistry pipeline over a directory of ``.xyz``/``.gen``
    structures or a ``.csv`` manifest, writing ``<out>.csv`` and ``<out>.json``:

    .. code-block:: bash

       # DFTB+ (3ob) with water solvation and quasi-RRHO entropy
       thermo screen molecules/ --parameter-set 3ob --solvent water --quasi-rrho

       # native xtb, radical anions in solution
       thermo screen molecules/ --engine xtb-cli --solvent water

       # resume an interrupted screen
       thermo screen molecules/ -o results --resume

Python API
----------

A single molecule through DFTB+:

.. code-block:: python

   from ase.build import molecule
   from ThermoScreening.thermo.api import dftbplus_thermo
   from ThermoScreening.calculator.dftbplus import dftb_3ob_parameters

   thermo = dftbplus_thermo(molecule("H2O"), **dftb_3ob_parameters)
   print(thermo.total_entropy("cal/(mol*K)"))       # ~45 cal/mol/K
   print(thermo.total_gibbs_free_energy("H"))

The GFN-xTB engines (no Slater-Koster files needed):

.. code-block:: python

   from ThermoScreening.thermo.api import xtb_thermo, xtb_cli_thermo

   gas = xtb_thermo(molecule("H2O"))                       # tblite, gas phase
   solv = xtb_cli_thermo(molecule("H2O"), solvent="water") # native xtb + ALPB

Conformers, then screen the ensemble:

.. code-block:: python

   from ThermoScreening.thermo import generate_conformers, write_conformers, screen

   conformers = generate_conformers("CCCCO", max_conformers=10)
   write_conformers(conformers, "butanol_confs")
   results = screen("butanol_confs", engine="xtb-cli")

Reactions and redox
-------------------

Once you have a ``Thermo`` object per species (from any engine, under the same
conditions), combine them into reaction free energies and reduction potentials:

.. code-block:: python

   from ase.build import molecule
   from ThermoScreening.thermo.api import xtb_cli_thermo
   from ThermoScreening.thermo.conformers import generate
   from ThermoScreening.thermo import reaction_free_energy, reduction_potential

   # reaction free energy (stoichiometry via (coefficient, Thermo) tuples)
   h2  = xtb_cli_thermo(molecule("H2"))
   o2  = xtb_cli_thermo(molecule("O2"), spin=1.0)   # S = 1, triplet
   h2o = xtb_cli_thermo(molecule("H2O"))
   dG = reaction_free_energy([(2, h2), o2], [(2, h2o)], unit="kcal")   # 2 H2 + O2 -> 2 H2O

   # one-electron reduction potential (Ox + e- -> Red), both in solvent
   mol = generate("O=C1C=CC(=O)C=C1", max_conformers=1)[0]   # benzoquinone
   neutral = xtb_cli_thermo(mol, charge=0,  solvent="water")
   anion   = xtb_cli_thermo(mol, charge=-1, solvent="water")   # radical anion, auto open-shell
   E = reduction_potential(neutral, anion)          # vs SHE (default reference 4.44 V)
   E_abs = reduction_potential(neutral, anion, reference_potential=0.0)

.. warning::

   ``reaction_free_energy`` / ``reduction_potential`` are exact given the input
   energies, but the *accuracy* is set by the underlying method. GFN-xTB and DFTB
   give poor **absolute** electron affinities and redox potentials (benzoquinone's
   GFN2 EA is ~7 eV vs ~1.9 eV experimental). Use them for **relative** trends
   across similar species, with a higher-accuracy method, or with a
   ``reference_potential`` calibrated against experiment.

End-to-end example
------------------

The ``examples/`` directory in the repository contains a complete DFTB+
thermochemistry input (``thermo.in``) together with the optimised geometry,
frequencies and energy it consumes. Run it directly with:

.. code-block:: bash

   thermo examples/thermo.in
