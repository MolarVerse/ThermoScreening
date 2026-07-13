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

       # DFTB+ (3ob) with water solvation, D3(BJ) dispersion and quasi-RRHO entropy
       thermo screen molecules/ --parameter-set 3ob --solvent water --dispersion d3-bj --quasi-rrho

       # native xtb, radical anions in solution
       thermo screen molecules/ --engine xtb-cli --solvent water

       # resume an interrupted screen
       thermo screen molecules/ -o results --resume

       # run 4 molecules at a time (each in its own process/directory)
       thermo screen molecules/ --jobs 4

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

Conformer ensembles
-------------------

A flexible molecule is a Boltzmann ensemble of conformers, not a single
structure. Generate and compute each conformer, then combine their free energies
into ensemble properties:

.. code-block:: python

   from ThermoScreening.thermo.api import xtb_cli_thermo
   from ThermoScreening.thermo.conformers import generate
   from ThermoScreening.thermo import (
       boltzmann_weights, ensemble_free_energy, lowest_gibbs,
   )

   conformers = generate("CCCCO", max_conformers=10)          # n-butanol
   thermos = [xtb_cli_thermo(c, solvent="water") for c in conformers]

   weights = boltzmann_weights(thermos)          # populations at 298.15 K
   G = ensemble_free_energy(thermos, unit="kcal")  # Boltzmann-averaged free energy
   best = lowest_gibbs(thermos)                   # the single lowest-G conformer

The ensemble free energy lies at or below the lowest conformer's by the mixing
(conformational) entropy; use ``lowest_gibbs`` when you instead want the single
dominant structure to carry forward (e.g. into :func:`reaction_free_energy`).

``pKa``, ``calibrate_proton_reference``, ``reduction_potential``, and
``reaction_free_energy`` all take whatever they're given and call
``.total_EeGtot()`` on it -- they don't care whether that's a single ``Thermo``
or something else with the same method. ``EnsembleThermo`` wraps a conformer
ensemble's Boltzmann free energy behind that same method, so a flexible acid
or base can be passed to ``pKa`` as an ensemble directly, instead of picking
one (possibly non-representative) conformer:

.. code-block:: python

   from ThermoScreening.thermo.api import xtb_cli_thermo
   from ThermoScreening.thermo.conformers import generate
   from ThermoScreening.thermo import EnsembleThermo, pKa

   def ensemble_thermo(smiles, charge):
       conformers = generate(smiles, max_conformers=10)
       thermos = [xtb_cli_thermo(c, charge=charge, solvent="water") for c in conformers]
       return EnsembleThermo(thermos)

   acid = ensemble_thermo("OCCCC(=O)O", 0)     # 4-hydroxybutanoic acid
   base = ensemble_thermo("OCCCC(=O)[O-]", -1)

   p_ka = pKa(acid, base)

In practice, some conformers -- especially of a charged species -- optimise
onto a saddle point rather than a true minimum instead of failing outright.
``generate_thermo_ensemble`` wraps the loop above with automatic retries: it
generates conformers, computes ``Thermo`` for each, skips any that raise
(imaginary frequencies, or the engine's process failing outright), and
reseeds and retries until enough succeed:

.. code-block:: python

   from ThermoScreening.thermo.api import xtb_cli_thermo
   from ThermoScreening.thermo import generate_thermo_ensemble, EnsembleThermo, pKa

   acid_thermos = generate_thermo_ensemble(
       "OCCCC(=O)O", xtb_cli_thermo, charge=0, max_conformers=10, solvent="water",
   )
   base_thermos = generate_thermo_ensemble(
       "OCCCC(=O)[O-]", xtb_cli_thermo, charge=-1, max_conformers=10, solvent="water",
   )

   p_ka = pKa(EnsembleThermo(acid_thermos), EnsembleThermo(base_thermos))

DFT-quality thermochemistry (ORCA)
----------------------------------

GFN-xTB/DFTB absolute energies are semiempirical. To get DFT-accurate
thermochemistry (e.g. for reliable reaction/redox free energies), feed an ORCA
frequency calculation's ``.hess`` file straight in:

.. code-block:: python

   from ThermoScreening.thermo.api import orca_thermo
   from ThermoScreening.thermo import reduction_potential

   ox  = orca_thermo("oxidized.hess")   # geometry, frequencies and energy read
   red = orca_thermo("reduced.hess")    # from the companion .property.txt file
   E = reduction_potential(ox, red)     # now on DFT energetics

The energy is read from the ``<basename>.property.txt`` file ORCA writes
alongside the ``.hess`` (both must be kept together); pass ``energy=`` to
override it (e.g. a higher-level single-point), or when only a bare ``.hess``
is available. ``read_orca_hess`` exposes the parsed geometry/frequencies/energy
directly if you need them.

For **Gaussian, Turbomole, Psi4, NWChem** (and ORCA) in one call, install the
``qm`` extra (``pip install thermoscreening[qm]``) and use ``cclib_thermo``,
which auto-detects the program via `cclib <https://cclib.github.io/>`_:

.. code-block:: python

   from ThermoScreening.thermo.api import cclib_thermo

   thermo = cclib_thermo("freq.log")        # Gaussian, ORCA, Psi4, NWChem, ...
   print(thermo.total_gibbs_free_energy())  # energy read from the output, in Hartree

   # Turbomole splits a job's output across many files instead of one logfile;
   # pass every relevant file's path as a list (cclib's multi-file mode)
   ts_thermo = cclib_thermo(["control", "coord", "aoforce.out"])

``read_cclib`` returns the parsed ``(atoms, frequencies, energy)`` if you want
them directly; ``energy=`` overrides the parsed energy.

**PySCF** results live in memory, so pass the mean-field object directly (with
the ``pyscf`` extra, ``pip install thermoscreening[pyscf]``):

.. code-block:: python

   from ThermoScreening.thermo.api import pyscf_thermo

   mf = mol.RKS(xc="b3lyp").run()
   hessian = mf.Hessian().kernel()
   thermo = pyscf_thermo(mf, hessian=hessian)   # or frequencies=<cm^-1 array>

The geometry and energy (``mf.e_tot``) come from the mean field; frequencies are
taken from ``frequencies=`` or derived from ``hessian=`` via ``pyscf.hessian``.

.. note::

   ``orca_thermo`` and ``cclib_thermo`` are tested against genuine program
   output (real ORCA 6.1.1, Gaussian 16 and Turbomole 7.2 calculations), not
   just synthetic files, and any file cclib or ``read_orca_hess`` cannot parse
   raises a ``TSValueError`` rather than an arbitrary exception from the
   underlying parser.

Temperature scans
-----------------

The electronic energy, geometry and frequencies do not depend on temperature,
so a single (expensive) DFTB+/Hessian run can be reused across many
temperatures:

.. code-block:: python

   thermo = dftbplus_thermo(molecule("H2O"), **dftb_3ob_parameters)
   scan = thermo.temperature_scan([250, 273.15, 298.15, 350, 400])
   gibbs = [(t._temperature, t.total_gibbs_free_energy("H")) for t in scan]

Each entry is a fully-evaluated ``Thermo`` at that temperature; the original
object is left unchanged.

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

Reactivity site prediction (Fukui indices)
-------------------------------------------

``reaction_free_energy`` and ``reduction_potential`` tell you *whether* a
reaction is favourable for the molecule as a whole; they say nothing about
*where* on the molecule it happens, or whether some other site might react
first (an unintended side reaction competing with the one you're modelling).
``xtb_fukui_indices`` runs a single ``xtb --vfukui`` calculation on an
already-optimised geometry and reports, per atom, how susceptible that site
is to gaining electron density (``f_plus``, i.e. reduction), losing it
(``f_minus``, oxidation), or either (``f_zero``, radical attack):

.. code-block:: python

   import ase.io
   from ThermoScreening.thermo.api import xtb_cli_thermo, xtb_fukui_indices
   from ThermoScreening.thermo.conformers import generate

   mol = generate("O=C1C=CC(=O)C=C1", max_conformers=1)[0]   # benzoquinone
   xtb_cli_thermo(mol, directory="opt")   # optimise; xtbopt.xyz is written to "opt"

   optimized = ase.io.read("opt/xtbopt.xyz")
   rows = xtb_fukui_indices(optimized)
   most_oxidisable = max(rows, key=lambda row: row[2])   # (symbol, f_plus, f_minus, f_zero)

Each of ``f_plus``/``f_minus``/``f_zero`` sums to ~1 over all atoms, so a
value well above the 1/n_atoms baseline flags that site as the dominant one
for that kind of reactivity.

Acid dissociation (pKa)
------------------------

Proton-coupled redox (e.g. hydroquinone/semiquinone/quinone protonation
states) needs a pKa alongside the reduction potential. ``pKa`` uses the
"direct method" thermodynamic cycle -- ``HA(soln) -> A-(soln) + H+(soln)`` --
combining the acid and conjugate base's computed free energies with a
literature reference free energy for the (uncomputable) aqueous proton:

.. code-block:: python

   from ThermoScreening.thermo.api import xtb_cli_thermo
   from ThermoScreening.thermo.conformers import generate
   from ThermoScreening.thermo import pKa, calibrate_proton_reference

   def acid_base_pair(acid_smiles, base_smiles):
       # the acid and its conjugate base are different structures (one fewer
       # H), not the same structure at a different charge (that would be
       # reduction, not deprotonation)
       acid = generate(acid_smiles, max_conformers=1)[0]
       base = generate(base_smiles, max_conformers=1)[0]
       return (
           xtb_cli_thermo(acid, charge=0, solvent="water"),
           xtb_cli_thermo(base, charge=-1, solvent="water"),
       )

   phenol, phenolate = acid_base_pair("Oc1ccccc1", "[O-]c1ccccc1")
   hq, hq_anion = acid_base_pair("Oc1ccc(O)cc1", "[O-]c1ccc(O)cc1")  # hydroquinone

   p_ka = pKa(hq, hq_anion)  # -99.45 -- see the warning below
   phenol_p_ka = pKa(phenol, phenolate)  # -100.46 vs. experimental 9.99

.. warning::

   The default proton reference (``PROTON_AQUEOUS_FREE_ENERGY_KCAL``) is a
   literature constant derived assuming DFT/ab-initio-quality absolute
   energies. With **GFN-xTB or DFTB it is not usable at all uncalibrated**:
   verified with the real xtb-cli calculation above, ``phenol_p_ka`` comes out
   at **-100.46** against phenol's experimental 9.99 --
   an error of ~110 pKa units, not the "several units" DFT+continuum-solvent
   methods show (Ho & Coote, *Theor. Chem. Acc.* **2010**, *125*, 3).
   Semiempirical tight-binding methods do not preserve an absolute
   ab-initio/experimental energy scale, so their energies cannot be mixed
   directly with this constant.

   ``calibrate_proton_reference`` against **one** reference compound absorbs
   this offset completely (it is an additive constant, not
   molecule-specific) and is **required** for these engines, not just
   recommended for extra accuracy:

   .. code-block:: python

      ref_g = calibrate_proton_reference(phenol, phenolate, experimental_pKa=9.99)
      p_ka = pKa(hq, hq_anion, reference_free_energy=ref_g)  # 11.00 vs. exp 10.35

Transition states and rate constants
-------------------------------------

By default ``Thermo`` requires a minimum (no imaginary frequencies). Pass
``transition_state=True`` to the ``*_thermo`` functions to evaluate a
first-order saddle point instead: its one required imaginary (reaction
coordinate) mode is excluded from the vibrational thermochemistry rather than
raising, and exposed via ``Thermo.imaginary_mode_wavenumber()``. This works
with any engine that imports an externally-computed structure -- ``orca_thermo``,
``cclib_thermo`` and ``pyscf_thermo`` -- since a DFTB+/xtb geometry
*optimization* is a minimizer and cannot itself locate a saddle point:

.. code-block:: python

   from ThermoScreening.thermo.api import orca_thermo
   from ThermoScreening.thermo import eyring_rate_constant, wigner_tunneling_correction

   reactant = orca_thermo("reactant.hess")
   ts = orca_thermo("ts.hess", transition_state=True)

   nu_imag = ts.imaginary_mode_wavenumber()          # cm^-1, negative
   kappa = wigner_tunneling_correction(nu_imag, temperature=298.15)
   k = eyring_rate_constant([reactant], ts, temperature=298.15, kappa=kappa)

``eyring_rate_constant`` uses the same reactant convention as
``reaction_free_energy`` (a list of ``Thermo`` or ``(coefficient, Thermo)``
entries), so a bimolecular reaction is ``eyring_rate_constant([a, b], ts)``.
For a single reactant the result is a first-order rate constant in s\ :sup:`-1`;
for multiple reactants it is the pseudo-first-order TST rate (no standard-state
correction is applied). ``wigner_tunneling_correction`` is a small, first-order
tunneling estimate -- valid only when it stays close to 1; for deep tunneling
use a more complete treatment.

End-to-end example
------------------

The ``examples/`` directory in the repository contains a complete DFTB+
thermochemistry input (``thermo.in``) together with the optimised geometry,
frequencies and energy it consumes. Run it directly with:

.. code-block:: bash

   thermo examples/thermo.in
