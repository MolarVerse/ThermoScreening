import math
import types

import numpy as np
import pytest

import cclib.io

from ThermoScreening.calculator import qm
from ThermoScreening.calculator.qm import read_cclib, _best_energy_ev
from ThermoScreening.thermo.api import cclib_thermo
from ThermoScreening.exceptions import TSValueError

_H_TO_EV = 27.211386245988

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
