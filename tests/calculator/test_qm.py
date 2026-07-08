import math
import types
from pathlib import Path

import numpy as np
import pytest

import cclib.io

from ThermoScreening.calculator import qm
from ThermoScreening.calculator.qm import (
    read_cclib,
    _best_energy_ev,
    _normalize_source,
    _describe_source,
)
from ThermoScreening.thermo.api import cclib_thermo
from ThermoScreening.exceptions import TSValueError

_H_TO_EV = 27.211386245988

_REAL_TURBOMOLE_FILES = sorted(
    str(p) for p in
    (Path(__file__).resolve().parents[1] / "data" / "calculator" / "turbomole").glob("*")
    if p.suffix != ".md"
)

_GAUSSIAN_DATA_DIR = Path(__file__).resolve().parents[1] / "data" / "calculator" / "gaussian"
_REAL_GAUSSIAN_WATER = str(_GAUSSIAN_DATA_DIR / "water_neutral_opt_freq.out")
_REAL_GAUSSIAN_UNPARSEABLE = str(_GAUSSIAN_DATA_DIR / "mp2_avdz_freq_tight.log")

# water: geometry (Angstrom), three real modes, SCF energy in eV (~ -76.4 Ha)
_WATER = dict(
    atomnos=np.array([8, 1, 1]),
    atomcoords=np.array([[[0.0, 0.0, 0.1173],
                          [0.0, 0.7572, -0.4692],
                          [0.0, -0.7572, -0.4692]]]),
    vibfreqs=np.array([1600.0, 3700.0, 3800.0]),
    scfenergies=np.array([-76.4 * _H_TO_EV]),
)


def _fake_ccread(monkeypatch, **data):
    monkeypatch.setattr(cclib.io, "ccread", lambda path: types.SimpleNamespace(**data))


def _fake_ccread_returns(monkeypatch, value):
    monkeypatch.setattr(cclib.io, "ccread", lambda path: value)


def test_read_cclib_atoms_frequencies_energy(monkeypatch, tmp_path):
    _fake_ccread(monkeypatch, **_WATER)
    atoms, freqs, energy = read_cclib(str(tmp_path / "water.log"))

    assert list(atoms.get_chemical_symbols()) == ["O", "H", "H"]
    assert atoms.positions[1, 1] == pytest.approx(0.7572)
    assert list(freqs) == [1600.0, 3700.0, 3800.0]
    assert energy == pytest.approx(-76.4)  # eV -> Hartree


def test_best_energy_prefers_cc_then_mp_then_scf():
    scf = types.SimpleNamespace(scfenergies=np.array([-1.0]))
    assert _best_energy_ev(scf) == -1.0

    mp = types.SimpleNamespace(scfenergies=np.array([-1.0]),
                               mpenergies=np.array([[-2.0, -2.5]]))
    assert _best_energy_ev(mp) == -2.5  # last step, highest MP order

    cc = types.SimpleNamespace(scfenergies=np.array([-1.0]),
                               mpenergies=np.array([[-2.0]]),
                               ccenergies=np.array([-3.0]))
    assert _best_energy_ev(cc) == -3.0

    assert _best_energy_ev(types.SimpleNamespace()) is None


def test_read_cclib_unparseable_raises(monkeypatch, tmp_path):
    _fake_ccread_returns(monkeypatch, None)
    with pytest.raises(TSValueError, match="could not parse"):
        read_cclib(str(tmp_path / "x.log"))


def test_read_cclib_without_frequencies_raises(monkeypatch, tmp_path):
    _fake_ccread(monkeypatch, atomnos=np.array([1]),
                 atomcoords=np.array([[[0.0, 0.0, 0.0]]]), vibfreqs=np.array([]))
    with pytest.raises(TSValueError, match="no vibrational frequencies"):
        read_cclib(str(tmp_path / "x.log"))


def test_read_cclib_missing_dependency(monkeypatch, tmp_path):
    def _raise():
        raise TSValueError("cclib is required to import QM outputs; install ...")

    monkeypatch.setattr(qm, "_import_cclib", _raise)
    with pytest.raises(TSValueError, match="cclib is required"):
        read_cclib(str(tmp_path / "x.log"))


def test_cclib_thermo_uses_file_energy(monkeypatch, tmp_path):
    _fake_ccread(monkeypatch, **_WATER)
    thermo = cclib_thermo(str(tmp_path / "water.log"))
    assert thermo.electronic_energy() == pytest.approx(-76.4)
    assert math.isfinite(thermo.total_EeGtot())


def test_cclib_thermo_transition_state(monkeypatch, tmp_path):
    ts_data = dict(_WATER, vibfreqs=np.array([-300.0, 1600.0, 3800.0]))
    _fake_ccread(monkeypatch, **ts_data)
    thermo = cclib_thermo(str(tmp_path / "ts.log"), transition_state=True)
    assert thermo.imaginary_mode_wavenumber() == pytest.approx(-300.0)


def test_cclib_thermo_energy_override(monkeypatch, tmp_path):
    _fake_ccread(monkeypatch, **_WATER)
    thermo = cclib_thermo(str(tmp_path / "water.log"), energy=-77.0)
    assert thermo.electronic_energy() == pytest.approx(-77.0)


def test_cclib_thermo_requires_energy(monkeypatch, tmp_path):
    data = dict(_WATER)
    del data["scfenergies"]  # no energy anywhere and none passed
    _fake_ccread(monkeypatch, **data)
    with pytest.raises(TSValueError, match="No energy"):
        cclib_thermo(str(tmp_path / "water.log"))


def test_cclib_thermo_requires_energy_multi_file_error_is_readable(monkeypatch):
    data = dict(_WATER)
    del data["scfenergies"]
    _fake_ccread(monkeypatch, **data)
    with pytest.raises(TSValueError, match=r"No energy for 'control, coord, aoforce\.out'"):
        cclib_thermo(["control", "coord", "aoforce.out"])


# --- multi-file (Turbomole) support --- #
#
# Turbomole splits a job's output across many small files instead of one
# logfile, so cclib.io.ccread must receive a LIST of paths for it, not a
# single path. read_cclib previously always did cclib.io.ccread(str(path)),
# which silently could not support this at all.


def test_normalize_source_single_path_is_a_string():
    assert _normalize_source("a.log") == "a.log"
    assert _normalize_source(Path("a.log")) == "a.log"


def test_normalize_source_list_stays_a_list_of_strings():
    assert _normalize_source(["control", Path("coord"), "aoforce.out"]) == [
        "control", "coord", "aoforce.out",
    ]


def test_describe_source_formats_both_forms():
    assert _describe_source("a.log") == "a.log"
    assert _describe_source(["a", "b"]) == "a, b"


def test_read_cclib_passes_a_list_through_to_ccread_unmodified(monkeypatch):
    captured = {}

    def fake_ccread(source):
        captured["source"] = source
        return types.SimpleNamespace(**_WATER)

    monkeypatch.setattr(cclib.io, "ccread", fake_ccread)
    files = ["control", "coord", "aoforce.out"]
    read_cclib(files)

    assert captured["source"] == files  # not str(files) / a single joined string


def test_read_cclib_real_turbomole_output():
    # a genuine Turbomole 7.2 aoforce (frequency) calculation; see
    # tests/data/calculator/turbomole/README.md for provenance
    atoms, freqs, energy = read_cclib(_REAL_TURBOMOLE_FILES)

    assert list(atoms.get_chemical_symbols()) == ["Cl", "Au", "N", "N", "Au", "N", "N"]
    assert len(freqs) == 15  # 3*7 - 6, all real (a genuine minimum)
    assert freqs.min() > 0
    assert freqs.min() == pytest.approx(17.75)
    assert freqs.max() == pytest.approx(2303.92)
    # cclib's scfenergies (eV) converted to Hartree, cross-checked against the
    # raw eV value cclib itself reports for this file (-25880.26134295 eV)
    assert energy == pytest.approx(-25880.26134295 / _H_TO_EV)


def test_cclib_thermo_real_turbomole_output():
    thermo = cclib_thermo(_REAL_TURBOMOLE_FILES)

    assert thermo.electronic_energy() == pytest.approx(-951.0820621320888)
    assert math.isfinite(thermo.total_EeGtot())
    assert thermo.total_entropy("cal/(mol*K)") > 0


# --- real Gaussian 16 output --- #
#
# See tests/data/calculator/gaussian/README.md for provenance.


def test_read_cclib_real_gaussian_output():
    atoms, freqs, energy = read_cclib(_REAL_GAUSSIAN_WATER)

    assert list(atoms.get_chemical_symbols()) == ["O", "H", "H"]
    assert list(freqs) == pytest.approx([2169.7613, 4141.3837, 4392.5759])
    assert energy == pytest.approx(-74.96589788571758)


def test_cclib_thermo_real_gaussian_output():
    thermo = cclib_thermo(_REAL_GAUSSIAN_WATER)

    assert thermo.electronic_energy() == pytest.approx(-74.96589788571758)
    # gas-phase water's experimental standard entropy is ~45 cal/(mol K)
    assert thermo.total_entropy("cal/(mol*K)") == pytest.approx(45.0, abs=2.0)
    assert math.isfinite(thermo.total_EeGtot())


def test_read_cclib_wraps_a_real_cclib_parser_crash():
    # this real Gaussian 16 file's "Leave Link" timing line crashes cclib
    # 1.8.1's own Gaussian parser with an unhandled ValueError (confirmed:
    # still the latest cclib release at the time this test was written) --
    # that must not leak out of read_cclib as a raw, arbitrary exception
    with pytest.raises(TSValueError, match="cclib raised"):
        read_cclib(_REAL_GAUSSIAN_UNPARSEABLE)


def test_read_cclib_wraps_any_ccread_exception(monkeypatch):
    # a cclib-version-independent guard for the same behaviour: whatever
    # cclib.io.ccread raises internally must come out as a clean TSValueError,
    # not the original exception type
    def fake_ccread(source):
        raise ValueError("some internal cclib parsing failure")

    monkeypatch.setattr(cclib.io, "ccread", fake_ccread)
    with pytest.raises(TSValueError, match="cclib raised ValueError.*some internal"):
        read_cclib("whatever.log")
