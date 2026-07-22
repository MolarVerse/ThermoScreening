[![Python application](https://github.com/MolarVerse/ThermoScreening/actions/workflows/python-app.yml/badge.svg)](https://github.com/MolarVerse/ThermoScreening/actions/workflows/python-app.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21473884.svg)](https://doi.org/10.5281/zenodo.21473884)
[![codecov](https://codecov.io/gh/MolarVerse/ThermoScreening/graph/badge.svg?token=KhrG0zVZmS)](https://codecov.io/gh/MolarVerse/ThermoScreening)
[![Docs](https://github.com/MolarVerse/ThermoScreening/actions/workflows/docs.yml/badge.svg)](https://molarverse.github.io/ThermoScreening/)

# ThermoScreening

ThermoScreening calculates thermochemical properties for molecular systems and provides a foundation for screening molecule sets. It currently supports thermochemistry workflows from DFTB+ inputs and exposes Python APIs for reading coordinates, parsing vibrational data, and running thermodynamic post-processing.

## Documentation

The documentation is published at **https://molarverse.github.io/ThermoScreening/**.

It (installation, usage, configuration, and the API reference) is built with
Sphinx from the `docs/` directory, and you can also build it locally:

```bash
python -m pip install -e ".[docs]"
python -m sphinx -b html docs docs/_build/html   # open docs/_build/html/index.html
```

## Features

- Thermochemistry calculations for molecular systems
- DFTB+ geometry optimization, Hessian, and normal-mode integration
- Stepwise and overall reference-calibrated reduction potentials
- Readers for DFTB+ `.gen`, XYZ, and vibrational frequency files
- Runtime type checking for public API calls
- Test coverage for parsing, thermochemistry, and optional DFTB+ execution paths

## Installation

Install the released package:

```bash
python -m pip install thermoscreening
```

For a complete environment with DFTB+, native xTB, and tblite:

```bash
conda create -n thermoscreening -c conda-forge \
    python=3.12 dftbplus xtb tblite-python pip
conda activate thermoscreening
python -m pip install thermoscreening
```

From a source checkout, the equivalent development environment is:

```bash
conda env create -f environment.yml
conda activate thermoscreening
```

Run `thermo doctor` to report the usable engines, or check one explicitly with
`thermo doctor --engine dftb+`.

## DFTB+ Setup

DFTB+ calculations require two external pieces:

1. The `dftb+` and `modes` executables on `PATH`.
2. Slater-Koster parameter files downloaded separately from DFTB.org.

Install DFTB+ with Conda if it is not already available:

```bash
conda install -c conda-forge dftbplus
```

Download the default `3ob-3-1` Slater-Koster files into the automatically
discovered user-local directory:

```bash
thermo setup-dftb
```

Verify the setup with:

```bash
thermo doctor
```

ThermoScreening does not vendor Slater-Koster files. For custom installations,
override the discovered directory with `DFTB_PREFIX` or pass `slako_dir`
explicitly:

```python
from ThermoScreening.thermo.api import dftbplus_thermo

thermo = dftbplus_thermo(
    atoms,
    slako_dir="/path/to/3ob-3-1/",
)
```

The bundled DFTB+ parameters were removed from the repository because they are large, independently licensed scientific data. Keeping them external makes the package smaller and keeps parameter-set licensing explicit.

## Usage

Run thermochemistry from an input file with the command-line entry point:

```bash
thermo run path/to/thermo.in
```

The historic `thermo path/to/thermo.in` form remains supported.

Run a complete reference-calibrated redox screen locally or as a Slurm array:

```bash
thermo redox molecules.csv --engine xtb-cli --solvent acetonitrile -o redox
thermo slurm --tasks 32 --submit -- redox molecules.csv -o redox
```

Use the Python API when integrating ThermoScreening into another workflow:

```python
from ThermoScreening.thermo.api import run_thermo

thermo = run_thermo(
    vibrational_frequencies,
    coord_file="geo_opt.xyz",
    temperature=298.15,
    pressure=101325,
    energy=electronic_energy,
    engine="dftb+",
)

print(thermo.total_gibbs_free_energy())
```

## Testing

Run the full test suite:

```bash
python -m pytest -q
```

Run linting:

```bash
python -m pylint ThermoScreening
```

DFTB+ integration tests run only when the executables are available and `DFTB_PREFIX` points to a valid Slater-Koster directory. Otherwise they are skipped so the pure-Python test suite remains portable.

## Citing

If you use ThermoScreening in research, cite the
[archived software release](https://doi.org/10.5281/zenodo.21473884).
Machine-readable citation metadata is available in [`CITATION.cff`](CITATION.cff),
which also powers GitHub's **Cite this repository** feature. Each GitHub release
is archived by Zenodo and receives a version-specific DOI.

## Roadmap

Planned work is tracked in GitHub issues rather than in this README. The tool
supports the DFTB+, GFN-xTB (tblite) and native-xtb engines, implicit solvation,
quasi-RRHO, batch screening with resume, and RDKit conformer generation; see the
issue tracker for further enhancements.

## License

ThermoScreening source code is licensed under the GNU Lesser General Public License v2.1 or later. See [LICENSE](LICENSE).
