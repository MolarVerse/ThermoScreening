import importlib.util
import math

import numpy as np
import pytest

_HAS_PYSCF = importlib.util.find_spec("pyscf") is not None

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
    # force the pyscf import inside the helper to fail, robustly -- independent
    # of whether pyscf is installed or already imported (which would defeat a
    # sys.modules sentinel) -- by intercepting the import itself
    import builtins

    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name.startswith("pyscf"):
            raise ImportError("no pyscf")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    with pytest.raises(TSValueError, match="pyscf is required"):
        _pyscf_frequencies_from_hessian(object(), object())


@pytest.mark.skipif(not _HAS_PYSCF, reason="pyscf is not installed")
def test_pyscf_frequencies_from_hessian_real():
    # real pyscf: harmonic_analysis on a water HF/STO-3G Hessian returns the
    # 3N-6 vibrational modes as a real (cm^-1) array (imaginary_freq=False)
    from pyscf import gto, scf

    mol = gto.M(atom="O 0 0 0.117; H 0 0.757 -0.469; H 0 -0.757 -0.469",
                basis="sto-3g", verbose=0)
    mf = scf.RHF(mol).run()
    hessian = mf.Hessian().kernel()

    frequencies = _pyscf_frequencies_from_hessian(mol, hessian)
    assert frequencies.dtype == float  # real (negative for any imaginary mode)
    assert len(frequencies) == 3       # 3N - 6 for a bent triatomic
