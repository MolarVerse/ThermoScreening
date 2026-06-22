[![Python application](https://github.com/MolarVerse/ThermoScreening/actions/workflows/python-app.yml/badge.svg)](https://github.com/MolarVerse/ThermoScreening/actions/workflows/python-app.yml)
[![codecov](https://codecov.io/gh/MolarVerse/ThermoScreening/graph/badge.svg?token=KhrG0zVZmS)](https://codecov.io/gh/MolarVerse/ThermoScreening)

# ThermoScreening

ThermoScreening calculates thermochemical properties for molecular systems and provides a foundation for screening molecule sets. It currently supports thermochemistry workflows from DFTB+ inputs and exposes Python APIs for reading coordinates, parsing vibrational data, and running thermodynamic post-processing.

## Features

- Thermochemistry calculations for molecular systems
- DFTB+ geometry optimization, Hessian, and normal-mode integration
- Readers for DFTB+ `.gen`, XYZ, and vibrational frequency files
- Runtime type checking for public API calls
- Test coverage for parsing, thermochemistry, and optional DFTB+ execution paths

## Installation

Install the package from a checkout:

```bash
python -m pip install .
```

For development and tests:

```bash
python -m pip install -e ".[test,lint]"
```

## DFTB+ Setup

DFTB+ calculations require two external pieces:

1. The `dftb+` and `modes` executables on `PATH`.
2. Slater-Koster parameter files downloaded separately from DFTB.org.

ThermoScreening does not vendor Slater-Koster files. Point the calculator to a parameter directory in one of two ways:

```bash
export DFTB_PREFIX=/path/to/3ob-3-1/
```

or pass `slako_dir` explicitly:

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
thermo path/to/thermo.in
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

## Roadmap

Planned work is tracked in GitHub issues rather than in this README. Current roadmap areas include additional engines, conformer generation, broader test coverage, documentation, and batch screening workflows.

## License

ThermoScreening source code is licensed under the GNU Lesser General Public License v2.1 or later. See [LICENSE](LICENSE).
