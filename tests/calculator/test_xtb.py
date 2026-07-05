import importlib.util

import numpy as np
import pytest

from ThermoScreening.calculator.xtb import (
    _eV_to_hartree,
    _real_frequencies_cm,
    optimise_and_frequencies,
)


def test_eV_to_hartree_conversion():
    # 1 Hartree is ~27.2114 eV
    assert _eV_to_hartree(27.211386) == pytest.approx(1.0, rel=1e-4)
    assert _eV_to_hartree(0.0) == 0.0


def test_real_frequencies_cm_maps_imaginary_negative_and_sorts():
    # imaginary modes (in the imaginary part) -> negative; then sorted ascending
    freqs = np.array([1000 + 0j, 0 + 50j, 200 + 0j, 0 + 5j])
    out = _real_frequencies_cm(freqs)

    assert list(out) == [-50.0, -5.0, 200.0, 1000.0]


def test_optimise_and_frequencies_with_emt(monkeypatch, tmp_path):
    # exercise the real ASE optimiser + finite-difference vibrations with a
    # dependency-free calculator (EMT), so no tblite is required
    from ase import Atoms
    from ase.calculators.emt import EMT

    monkeypatch.chdir(tmp_path)
    cu2 = Atoms("Cu2", positions=[[0, 0, 0], [0, 0, 2.4]])

    optimized, energy_hartree, frequencies = optimise_and_frequencies(
        cu2, EMT(), fmax=0.05
    )

    assert isinstance(energy_hartree, float)
    assert len(frequencies) == 6  # 3N modes for 2 atoms
    assert np.all(np.diff(frequencies) >= 0)  # sorted ascending
    # the single real stretch (top mode) is a positive vibration
    assert frequencies[-1] > 0


tblite_available = importlib.util.find_spec("tblite") is not None


@pytest.mark.skipif(not tblite_available, reason="tblite (GFN-xTB) is not installed.")
def test_xtb_thermo_runs_real_gfn2(tmp_path):
    # real GFN2-xTB end-to-end: gas-phase water entropy is close to experiment
    from ase.build import molecule
    from ThermoScreening.thermo.api import xtb_thermo

    thermo = xtb_thermo(molecule("H2O"), directory=str(tmp_path / "w"))
    entropy = thermo.total_entropy("cal/(mol*K)")

    # experimental standard molar entropy of gaseous water ~ 45.1 cal/mol/K
    assert entropy == pytest.approx(45.1, abs=2.0)
