import math
import sys

import numpy as np
import pytest

import ThermoScreening.thermo.api as api
from ThermoScreening.thermo.api import pyscf_thermo, _pyscf_frequencies_from_hessian
from ThermoScreening.exceptions import TSValueError


class _FakeMol:
    natm = 3

    def atom_coords(self):  # Bohr
        return np.array([[0.0, 0.0, 0.221722],
                         [0.0, 1.430901, -0.886569],
                         [0.0, -1.430901, -0.886569]])

    def atom_pure_symbol(self, index):
        return ["O", "H", "H"][index]


class _FakeMeanField:
    def __init__(self, e_tot=-76.4):
        self.mol = _FakeMol()
        self.e_tot = e_tot


def test_pyscf_thermo_from_frequencies():
    thermo = pyscf_thermo(_FakeMeanField(), frequencies=[1600.0, 3700.0, 3800.0])
    assert thermo.electronic_energy() == pytest.approx(-76.4)  # from mf.e_tot
    assert math.isfinite(thermo.total_EeGtot())


def test_pyscf_thermo_energy_override():
    thermo = pyscf_thermo(_FakeMeanField(), frequencies=[1600.0, 3700.0, 3800.0], energy=-77.0)
    assert thermo.electronic_energy() == pytest.approx(-77.0)


def test_pyscf_thermo_requires_frequencies_or_hessian():
    with pytest.raises(TSValueError, match="Provide frequencies"):
        pyscf_thermo(_FakeMeanField())


def test_pyscf_thermo_hessian_path(monkeypatch):
    monkeypatch.setattr(
        api, "_pyscf_frequencies_from_hessian",
        lambda mol, hessian: np.array([1600.0, 3700.0, 3800.0]),
    )
    thermo = pyscf_thermo(_FakeMeanField(), hessian=object())
    assert math.isfinite(thermo.total_EeGtot())


def test_pyscf_frequencies_from_hessian_missing_dependency(monkeypatch):
    # force the pyscf import inside the helper to fail, regardless of whether
    # pyscf happens to be installed in the test environment (None -> ImportError)
    monkeypatch.setitem(sys.modules, "pyscf.hessian.thermo", None)
    with pytest.raises(TSValueError, match="pyscf is required"):
        _pyscf_frequencies_from_hessian(object(), object())
