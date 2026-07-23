Gas-phase entropy accuracy
==========================

Standard molar entropies S°(298.15 K, 1 bar) from ThermoScreening compared with
experiment, for both engines. These are method accuracies (GFN2-xTB and
DFTB3/3ob), not tool errors. The thermochemistry implementation itself matches
ASE ``IdealGasThermo`` to within 0.0003 cal/(mol K).

.. list-table::
   :header-rows: 1

   * - Molecule
     - Symmetry number
     - xTB (GFN2)
     - DFTB+ 3ob
     - Experiment
   * - H₂O
     - 2
     - 45.08
     - 45.07
     - 45.1
   * - CH₄
     - 12
     - 44.44
     - 44.47
     - 44.5
   * - NH₃
     - 3
     - 45.96
     - 46.00
     - 46.0
   * - N₂
     - 2
     - 45.79
     - 45.76
     - 45.8
   * - CO₂
     - 2
     - 51.23
     - 51.53
     - 51.1
   * - CO
     - 1
     - 47.23
     - 47.14
     - 47.2
   * - CH₃OH
     - 1
     - 56.63
     - 57.00
     - 57.3
   * - C₂H₆
     - 6
     - 54.38
     - 54.54
     - 54.8

Values are in cal/(mol K).

Deviation from experiment
-------------------------

.. list-table::
   :header-rows: 1

   * - Engine
     - MAD
     - RMSD
     - Maximum absolute deviation
     - Mean signed deviation
   * - xTB (GFN2)
     - 0.17
     - 0.28
     - 0.67 (CH₃OH)
     - -0.13
   * - DFTB+ 3ob
     - 0.14
     - 0.21
     - 0.43 (CO₂)
     - -0.04

Both engines reproduce experimental gas-phase entropies to well under
1 cal/(mol K).

Notes
-----

- Symmetry numbers are detected correctly for all eight molecules.
- The values use ``pressure=100000`` (1 bar). The default pressure of 101325 Pa
  raises the entropy by 0.026 cal/(mol K).
- Residual error for methanol and ethane comes from treating low-frequency
  torsions as harmonic oscillators; ``quasi_rrho=True`` softens this.
- The Linux release gate runs a mandatory DFTB+ calculation. Broader xTB and
  solvent checks remain optional integration tests.
