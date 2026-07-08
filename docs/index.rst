ThermoScreening
===============

.. image:: https://img.shields.io/pypi/v/thermoscreening.svg
   :target: https://pypi.org/project/thermoscreening/
   :alt: PyPI version

.. image:: https://img.shields.io/pypi/pyversions/thermoscreening.svg
   :target: https://pypi.org/project/thermoscreening/
   :alt: Supported Python versions

.. image:: https://img.shields.io/pypi/l/thermoscreening.svg
   :target: https://github.com/MolarVerse/ThermoScreening/blob/main/LICENSE
   :alt: License

.. image:: https://github.com/MolarVerse/ThermoScreening/actions/workflows/python-app.yml/badge.svg
   :target: https://github.com/MolarVerse/ThermoScreening/actions/workflows/python-app.yml
   :alt: Build status

.. image:: https://codecov.io/gh/MolarVerse/ThermoScreening/graph/badge.svg?token=KhrG0zVZmS
   :target: https://codecov.io/gh/MolarVerse/ThermoScreening
   :alt: Code coverage

ThermoScreening computes molecular **thermochemistry** and screens sets of
molecules. It drives DFTB+ and GFN-xTB (or imports DFT results from ORCA,
Gaussian, Turbomole and PySCF), runs geometry optimisation, Hessian and
normal-mode analysis, and evaluates the ideal-gas thermochemistry --
entropy, enthalpy, Gibbs free energy, heat capacity -- with implicit
solvation, quasi-RRHO, conformer ensembles, reaction/redox free energies,
and transition-state kinetics.

.. code-block:: python

   from ase.build import molecule
   from ThermoScreening.thermo.api import dftbplus_thermo
   from ThermoScreening.calculator.dftbplus import dftb_3ob_parameters

   thermo = dftbplus_thermo(molecule("H2O"), **dftb_3ob_parameters)
   print(thermo.total_entropy("cal/(mol*K)"))   # ~45 cal/mol/K
   print(thermo.total_gibbs_free_energy("H"))   # Hartree

.. code-block:: bash

   $ thermo screen molecules/ --parameter-set 3ob --solvent water --quasi-rrho --jobs 4

New to ThermoScreening? Start with :doc:`installation` and :doc:`usage`.

Features
--------

.. grid:: 1 2 2 3
   :gutter: 3

   .. grid-item-card:: :octicon:`cpu` Engines
      :link: configuration
      :link-type: doc

      DFTB+ (3ob / mio), GFN-xTB via tblite, and the native ``xtb`` binary
      (open-shell radicals, charge, implicit solvation).

      Or import DFT-quality results from **ORCA**, **Gaussian**,
      **Turbomole** (via cclib) and **PySCF**.

   .. grid-item-card:: :octicon:`beaker` Thermochemistry
      :link: usage
      :link-type: doc

      Entropy, enthalpy, Gibbs free energy and heat capacity, validated
      against ASE ``IdealGasThermo`` and native xtb. GBSA/ALPB implicit
      solvation, spin polarisation, and Grimme's quasi-RRHO entropy.

   .. grid-item-card:: :octicon:`git-branch` Reactions & redox
      :link: usage
      :link-type: doc

      Reaction free energies and one-electron reduction potentials
      (vs. SHE or a calibrated reference) from any two ``Thermo`` objects.

   .. grid-item-card:: :octicon:`flame` Kinetics
      :link: usage
      :link-type: doc

      Transition-state thermochemistry, Eyring rate constants, and a
      Wigner tunneling correction from a saddle point's imaginary mode.

   .. grid-item-card:: :octicon:`list-unordered` Batch screening
      :link: usage
      :link-type: doc

      Screen a directory or CSV manifest in parallel (``--jobs``), with
      per-molecule error isolation, ``--resume``, and ranking by Gibbs
      free energy.

   .. grid-item-card:: :octicon:`versions` Conformers & ensembles
      :link: usage
      :link-type: doc

      Generate conformers from SMILES (RDKit ETKDG) and combine them into
      Boltzmann-weighted ensemble thermochemistry.

.. toctree::
   :maxdepth: 2
   :caption: Contents
   :hidden:

   installation
   usage
   configuration
   api

Indices
-------

* :ref:`genindex`
* :ref:`modindex`
