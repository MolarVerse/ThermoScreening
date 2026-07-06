ThermoScreening
===============

ThermoScreening computes molecular thermochemistry and screens sets of molecules.
It drives DFTB+ and GFN-xTB, runs geometry optimisation, Hessian and normal-mode
analysis, and evaluates the ideal-gas thermochemistry (entropy, enthalpy, Gibbs
free energy, heat capacity), with implicit solvation, quasi-RRHO, and conformer
generation.

Features
--------

- **Engines**: DFTB+ (3ob / mio parameter sets), GFN-xTB via tblite
  (``--engine xtb``), and the native ``xtb`` binary (``--engine xtb-cli``,
  which adds implicit solvation of charged radicals).
- **Thermochemistry** validated against ASE ``IdealGasThermo`` and native xtb.
- **Implicit solvation** (GBSA/ALPB), **spin polarisation** for open-shell
  species, and **Grimme quasi-RRHO** vibrational entropy.
- **Batch screening** with per-molecule error isolation and ``--resume``.
- **Conformer generation** from SMILES (RDKit ETKDG).

.. toctree::
   :maxdepth: 2
   :caption: Contents

   installation
   usage
   configuration
   api

Indices
-------

* :ref:`genindex`
* :ref:`modindex`
