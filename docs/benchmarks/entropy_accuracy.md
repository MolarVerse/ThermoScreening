# Gas-phase entropy accuracy

Standard molar entropies S°(298.15 K, 1 bar) from ThermoScreening compared with
experiment, for both engines. These are *method* accuracies (GFN2-xTB and
DFTB3/3ob), not tool errors — the tool's thermochemistry itself matches ASE
`IdealGasThermo` to < 0.0003 cal/mol/K (a non-circular cross-check).

| Molecule | σ | xTB (GFN2) | DFTB+ 3ob | Experiment |
|----------|---|-----------|-----------|------------|
| H₂O      | 2 | 45.08 | 45.07 | 45.1 |
| CH₄      | 12| 44.44 | 44.47 | 44.5 |
| NH₃      | 3 | 45.96 | 46.00 | 46.0 |
| N₂       | 2 | 45.79 | 45.76 | 45.8 |
| CO₂      | 2 | 51.23 | 51.53 | 51.1 |
| CO       | 1 | 47.23 | 47.14 | 47.2 |
| CH₃OH    | 1 | 56.63 | 57.00 | 57.3 |
| C₂H₆     | 6 | 54.38 | 54.54 | 54.8 |

(cal/mol/K)

**Deviation from experiment**

| Engine     | MAD  | RMSD | max\|dev\|      | mean signed |
|------------|------|------|----------------|-------------|
| xTB (GFN2) | 0.17 | 0.28 | 0.67 (CH₃OH)   | −0.13 |
| DFTB+ 3ob  | 0.14 | 0.21 | 0.43 (CO₂)     | −0.04 |

Both engines reproduce experimental gas-phase entropies to well under
1 cal/mol/K.

## Notes

- **Symmetry numbers** are auto-detected correctly for all eight molecules
  (H₂O = 2, CH₄ = 12, C₂H₆ = 6, …).
- **Standard-state convention.** The values above use `pressure=100000` (1 bar,
  the S° convention). The tool's default `pressure=101325` (1 atm) raises the
  entropy by `R·ln(101325/100000) = 0.026 cal/mol/K` — pass `pressure=100000`
  to compare with tabulated S°.
- **Residual error** on CH₃OH and C₂H₆ comes from treating low-frequency
  torsions as harmonic oscillators; `quasi_rrho=True` softens this.
- These numbers require the xtb/DFTB+ engines to regenerate optimized
  geometries and frequencies, so they live as skippable integration tests
  (`tests/calculator/test_xtb.py::test_xtb_thermo_entropy_vs_experiment`,
  `tests/calculator/test_dftbplus.py::TestDftbplus::test_dftbplus_entropy_vs_experiment`)
  rather than in the CI-safe suite.
