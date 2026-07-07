import math

import pytest
from ase.units import Bohr

from ThermoScreening.calculator.orca import read_orca_hess
from ThermoScreening.thermo.api import orca_thermo
from ThermoScreening.exceptions import TSValueError

# a minimal water ORCA .hess (coordinates in Bohr, frequencies in cm^-1)
_HESS = """\
$orca_hessian_file

$act_energy
     -76.400000

$vibrational_frequencies
9
  0     0.000000
  1     0.000000
  2     0.000000
  3     0.000000
  4     0.000000
  5     0.000000
  6  1600.000000
  7  3700.000000
  8  3800.000000

$atoms
3
O    15.999000    0.000000    0.000000    0.221722
H     1.008000    0.000000    1.430901   -0.886569
H     1.008000    0.000000   -1.430901   -0.886569

$end
"""


def _write_hess(tmp_path, text=_HESS, name="water.hess"):
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return str(path)


def test_read_orca_hess_geometry_frequencies_energy(tmp_path):
    atoms, freqs, energy = read_orca_hess(_write_hess(tmp_path))

    assert list(atoms.get_chemical_symbols()) == ["O", "H", "H"]
    # coordinates converted Bohr -> Angstrom
    assert atoms.positions[0, 2] == pytest.approx(0.221722 * Bohr)
    assert atoms.positions[1, 1] == pytest.approx(1.430901 * Bohr)
    assert len(freqs) == 9
    assert list(freqs[-3:]) == [1600.0, 3700.0, 3800.0]
    assert energy == pytest.approx(-76.4)


def test_read_orca_hess_missing_atoms(tmp_path):
    path = _write_hess(tmp_path, "$vibrational_frequencies\n1\n0 0.0\n$end\n", "b.hess")
    with pytest.raises(TSValueError, match=r"No \$atoms block"):
        read_orca_hess(path)


def test_read_orca_hess_missing_frequencies(tmp_path):
    path = _write_hess(tmp_path, "$atoms\n1\nH 1.008 0 0 0\n$end\n", "b.hess")
    with pytest.raises(TSValueError, match=r"No \$vibrational_frequencies block"):
        read_orca_hess(path)


def test_read_orca_hess_malformed_atom_count(tmp_path):
    text = _HESS.replace("$atoms\n3\n", "$atoms\n5\n")  # declares 5, lists 3
    with pytest.raises(TSValueError, match="declares 5 atoms"):
        read_orca_hess(_write_hess(tmp_path, text, "b.hess"))


def test_orca_thermo_uses_file_energy(tmp_path):
    thermo = orca_thermo(_write_hess(tmp_path))
    assert thermo.electronic_energy() == pytest.approx(-76.4)
    assert math.isfinite(thermo.total_EeGtot())


def test_orca_thermo_energy_override(tmp_path):
    thermo = orca_thermo(_write_hess(tmp_path), energy=-77.0)
    assert thermo.electronic_energy() == pytest.approx(-77.0)


def test_orca_thermo_requires_energy(tmp_path):
    # remove the $act_energy block (disable the tag) and pass no energy
    text = _HESS.replace("$act_energy", "$x_act_energy")
    with pytest.raises(TSValueError, match="No energy"):
        orca_thermo(_write_hess(tmp_path, text, "noenergy.hess"))
