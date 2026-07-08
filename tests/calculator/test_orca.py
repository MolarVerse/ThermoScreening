import math
import shutil
from pathlib import Path

import pytest
from ase.units import Bohr

from ThermoScreening.calculator.orca import read_orca_hess
from ThermoScreening.thermo.api import orca_thermo
from ThermoScreening.exceptions import TSValueError

_REAL_HESS = (
    Path(__file__).resolve().parents[1] / "data" / "calculator" / "orca" / "water_freq.hess"
)

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


def test_read_orca_hess_malformed_atom_coordinate(tmp_path):
    text = (
        "$atoms\n1\nO 15.999 notanumber 0.0 0.0\n"
        "$vibrational_frequencies\n1\n0 0.0\n$end\n"
    )
    with pytest.raises(TSValueError, match=r"Malformed \$atoms block"):
        read_orca_hess(_write_hess(tmp_path, text, "b.hess"))


def test_read_orca_hess_malformed_frequency(tmp_path):
    text = (
        "$atoms\n1\nH 1.008 0.0 0.0 0.0\n"
        "$vibrational_frequencies\n1\n0 notanumber\n$end\n"
    )
    with pytest.raises(TSValueError, match=r"Malformed \$vibrational_frequencies block"):
        read_orca_hess(_write_hess(tmp_path, text, "b.hess"))


def test_read_orca_hess_frequency_count_mismatch(tmp_path):
    text = (
        "$atoms\n1\nH 1.008 0.0 0.0 0.0\n"
        "$vibrational_frequencies\n5\n0 0.0\n$end\n"  # declares 5, lists 1
    )
    with pytest.raises(TSValueError, match="declares 5 entries"):
        read_orca_hess(_write_hess(tmp_path, text, "b.hess"))


def test_orca_thermo_uses_file_energy(tmp_path):
    thermo = orca_thermo(_write_hess(tmp_path))
    assert thermo.electronic_energy() == pytest.approx(-76.4)
    assert math.isfinite(thermo.total_EeGtot())


def test_orca_thermo_energy_override(tmp_path):
    thermo = orca_thermo(_write_hess(tmp_path), energy=-77.0)
    assert thermo.electronic_energy() == pytest.approx(-77.0)


def test_orca_thermo_transition_state(tmp_path):
    # mode 6 (the softest "vibrational" mode) is imaginary -> a TS
    ts_hess = _HESS.replace("  6  1600.000000", "  6  -300.000000")
    thermo = orca_thermo(_write_hess(tmp_path, ts_hess, "ts.hess"), transition_state=True)
    assert thermo.imaginary_mode_wavenumber() == pytest.approx(-300.0)


_CO2_HESS = """\
$act_energy
    -188.500000

$vibrational_frequencies
9
  0     0.000000
  1     0.000000
  2     0.000000
  3     0.000000
  4     0.000000
  5   667.000000
  6   667.000000
  7  1333.000000
  8  2349.000000

$atoms
3
C    12.011000    0.000000    0.000000    0.000000
O    15.999000    0.000000    0.000000    2.196000
O    15.999000    0.000000    0.000000   -2.196000

$end
"""


def test_orca_thermo_linear_molecule(tmp_path):
    # CO2 is linear -> 5 zero modes, dof = 3N-5 = 4; exercises the
    # order-sensitive frequency_dof / linear-rotor path end to end
    thermo = orca_thermo(_write_hess(tmp_path, _CO2_HESS, "co2.hess"))
    assert thermo.electronic_energy() == pytest.approx(-188.5)
    assert math.isfinite(thermo.total_EeGtot())


def test_orca_thermo_requires_energy(tmp_path):
    # remove the $act_energy block (disable the tag) and pass no energy
    text = _HESS.replace("$act_energy", "$x_act_energy")
    with pytest.raises(TSValueError, match="No energy"):
        orca_thermo(_write_hess(tmp_path, text, "noenergy.hess"))


# --- fixtures captured from a real ORCA 6.1.1 run (water, HF/STO-3G, Freq) --- #
#
# $act_energy in the .hess is 0.000000 for this ordinary (non-scan) job -- it
# is a relaxed-surface-scan field, not a general energy field -- so these
# fixtures are the regression test for that: the real energy only lives in the
# companion .property.txt ($Single_Point_Data / &FinalEnergy).

_REAL_ENERGY = -74.963023139027911  # ORCA's own "FINAL SINGLE POINT ENERGY"


def test_read_orca_hess_real_orca_output():
    atoms, freqs, energy = read_orca_hess(str(_REAL_HESS))

    assert list(atoms.get_chemical_symbols()) == ["O", "H", "H"]
    assert len(freqs) == 9
    assert list(freqs[:6]) == [0.0] * 6  # trans/rot modes
    assert list(freqs[6:]) == pytest.approx(
        [2043.0843558341605, 4488.043025859588, 4790.282399478668]
    )
    # the fix under test: $act_energy alone is 0.0 for this job (verified in
    # the captured file); the real energy must come from .property.txt
    assert energy == pytest.approx(_REAL_ENERGY)


def test_orca_thermo_real_orca_output():
    thermo = orca_thermo(str(_REAL_HESS))

    assert thermo.electronic_energy() == pytest.approx(_REAL_ENERGY)
    # gas-phase water's experimental standard entropy is ~45 cal/(mol K); HF/
    # STO-3G is a coarse method but should land in the right ballpark
    assert thermo.total_entropy("cal/(mol*K)") == pytest.approx(45.0, abs=5.0)
    assert math.isfinite(thermo.total_EeGtot())


def test_read_orca_hess_zero_act_energy_without_property_file_is_none(tmp_path):
    # the exact real-world bug: a .hess with a literal 0.0 $act_energy (as
    # ORCA writes for an ordinary job) and no companion .property.txt must not
    # be silently read as a real zero-Hartree energy
    text = _HESS.replace("-76.400000", "0.000000")
    _, _, energy = read_orca_hess(_write_hess(tmp_path, text, "zero.hess"))
    assert energy is None


def test_read_orca_hess_property_file_overrides_act_energy(tmp_path):
    # copy the real fixture pair under a different name and perturb
    # $act_energy in the .hess -- the .property.txt value must still win
    hess_path = tmp_path / "renamed.hess"
    shutil.copy(_REAL_HESS, hess_path)
    shutil.copy(_REAL_HESS.with_suffix(".property.txt"), tmp_path / "renamed.property.txt")

    text = hess_path.read_text(encoding="utf-8").replace("0.000000", "-999.0", 1)
    hess_path.write_text(text, encoding="utf-8")

    _, _, energy = read_orca_hess(str(hess_path))
    assert energy == pytest.approx(_REAL_ENERGY)  # not -999.0


def test_read_orca_hess_falls_back_to_nonzero_act_energy_without_property_file(tmp_path):
    # a genuinely nonzero $act_energy (e.g. a relaxed scan) with no companion
    # .property.txt is still used -- this is the pre-existing fallback path
    _, _, energy = read_orca_hess(_write_hess(tmp_path))  # default _HESS: -76.4
    assert energy == pytest.approx(-76.4)


def test_read_orca_hess_property_file_without_final_energy_falls_back(tmp_path):
    # a .property.txt with no &FinalEnergy entry (e.g. an unusual job type)
    # falls back to the .hess's own (nonzero) $act_energy
    hess_path = Path(_write_hess(tmp_path))
    hess_path.with_suffix(".property.txt").write_text(
        "$Some_Other_Data\n   &Nothing 1\n$End\n", encoding="utf-8"
    )
    _, _, energy = read_orca_hess(str(hess_path))
    assert energy == pytest.approx(-76.4)
