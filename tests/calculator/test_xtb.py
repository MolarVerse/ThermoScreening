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


def test_optimise_and_frequencies_cleans_vibration_cache_on_failure(
    monkeypatch, tmp_path
):
    from ase import Atoms
    from ase.calculators.emt import EMT

    clean_calls = []

    class FakeOptimizer:
        def __init__(self, atoms, logfile=None):
            self.atoms = atoms

        def run(self, fmax):
            return None

    class FakeVibrations:
        def __init__(self, atoms, name):
            self.name = name

        def clean(self):
            clean_calls.append(self.name)

        def run(self):
            raise RuntimeError("interrupted")

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr("ase.optimize.BFGS", FakeOptimizer)
    monkeypatch.setattr("ase.vibrations.Vibrations", FakeVibrations)

    atoms = Atoms("Cu", positions=[[0.0, 0.0, 0.0]])
    with pytest.raises(RuntimeError, match="interrupted"):
        optimise_and_frequencies(atoms, EMT())

    assert clean_calls == ["xtb_vib", "xtb_vib"]


def test_xtb_thermo_has_no_solvation_parameter():
    # regression guard: xTB-side implicit solvation is intentionally not wired up
    import inspect
    from ThermoScreening.thermo.api import xtb_thermo

    assert "solvent" not in inspect.signature(xtb_thermo).parameters


def test_xtb_thermo_passes_charge_and_unpaired_electrons_to_tblite(
    monkeypatch, tmp_path
):
    from ase import Atoms
    import ThermoScreening.thermo.api as api

    captured = {}
    sentinel = object()

    def fake_optimise(atoms, calculator, fmax):
        captured["charge"] = atoms.get_initial_charges().sum()
        captured["unpaired"] = atoms.get_initial_magnetic_moments().sum()
        return atoms, -1.0, np.array([100.0, 200.0, 300.0])

    monkeypatch.setattr(api, "optimise_and_frequencies", fake_optimise)
    monkeypatch.setattr(api, "xtb_calculator", lambda method: object())
    monkeypatch.setattr(api, "run_thermo", lambda *args, **kwargs: sentinel)

    result = api.xtb_thermo(
        Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.75]]),
        charge=-1,
        spin=0.5,
        directory=tmp_path,
    )

    assert result is sentinel
    assert captured["charge"] == pytest.approx(-1.0)
    assert captured["unpaired"] == pytest.approx(1.0)


tblite_available = importlib.util.find_spec("tblite") is not None
xtb_skip = pytest.mark.skipif(not tblite_available, reason="tblite (GFN-xTB) is not installed.")


@xtb_skip
@pytest.mark.parametrize("name,exp_S", [("H2O", 45.1), ("CH4", 44.5), ("N2", 45.8)])
def test_xtb_thermo_entropy_vs_experiment(name, exp_S, tmp_path):
    # real GFN2-xTB gas-phase entropy is within ~2 cal/mol/K of experiment
    from ase.build import molecule
    from ThermoScreening.thermo.api import xtb_thermo

    thermo = xtb_thermo(molecule(name), pressure=100000.0, directory=str(tmp_path / name))
    assert thermo.total_entropy("cal/(mol*K)") == pytest.approx(exp_S, abs=2.0)


# n-octane geometry (embedded so the test needs only tblite, not a 3D builder)
_OCTANE = [
    ("C", (4.1592, -0.4455, -0.0648)), ("C", (2.7620, 0.1265, 0.1227)),
    ("C", (1.7204, -0.9840, 0.2660)), ("C", (0.3124, -0.4653, 0.5771)),
    ("C", (-0.3124, 0.3231, -0.5771)), ("C", (-1.7204, 0.8418, -0.2660)),
    ("C", (-2.7620, -0.2686, -0.1227)), ("C", (-4.1592, 0.3033, 0.0648)),
    ("H", (4.8889, 0.3636, -0.1690)), ("H", (4.4517, -1.0575, 0.7944)),
    ("H", (4.2106, -1.0677, -0.9640)), ("H", (2.5272, 0.7621, -0.7373)),
    ("H", (2.7520, 0.7646, 1.0138)), ("H", (1.6959, -1.5906, -0.6475)),
    ("H", (2.0238, -1.6535, 1.0806)), ("H", (0.3407, 0.1568, 1.4797)),
    ("H", (-0.3149, -1.3331, 0.8095)), ("H", (-0.3407, -0.2989, -1.4797)),
    ("H", (0.3149, 1.1909, -0.8095)), ("H", (-1.6959, 1.4484, 0.6475)),
    ("H", (-2.0238, 1.5113, -1.0806)), ("H", (-2.5272, -0.9043, 0.7374)),
    ("H", (-2.7520, -0.9068, -1.0137)), ("H", (-4.8889, -0.5058, 0.1691)),
    ("H", (-4.4517, 0.9153, -0.7944)), ("H", (-4.2106, 0.9255, 0.9640)),
]


@xtb_skip
def test_xtb_quasi_rrho_on_real_floppy_molecule(tmp_path):
    # a real molecule with engine-generated sub-100 cm^-1 modes; quasi-RRHO tames
    # them, and (optimising once) leaves enthalpy and Cv invariant.
    import numpy as np
    from ase import Atoms
    from ThermoScreening.calculator.xtb import optimise_and_frequencies, xtb_calculator
    from ThermoScreening.thermo.api import run_thermo

    octane = Atoms([s for s, _ in _OCTANE], positions=[p for _, p in _OCTANE])
    octane.info["charge"] = 0
    octane.info["spin"] = 0
    opt, energy, freqs = optimise_and_frequencies(
        octane.copy(), xtb_calculator("GFN2-xTB"), fmax=0.02
    )

    kept = np.sort(freqs)[-(3 * len(opt) - 6):]
    assert np.any(kept < 100.0)  # genuine engine-generated low-frequency modes

    common = dict(atoms=opt, energy=energy, engine="xtb", pressure=101325)
    harmonic = run_thermo(freqs, quasi_rrho=False, **common)
    qrrho = run_thermo(freqs, quasi_rrho=True, **common)

    assert qrrho.total_entropy("cal/(mol*K)") < harmonic.total_entropy("cal/(mol*K)") - 2.0
    # invariance (same optimised freqs -> only entropy differs)
    assert qrrho.total_enthalpy("H") == pytest.approx(harmonic.total_enthalpy("H"), abs=1e-12)
    assert qrrho.total_heat_capacity("cal/(mol*K)") == pytest.approx(
        harmonic.total_heat_capacity("cal/(mol*K)"), abs=1e-12
    )


@xtb_skip
def test_xtb_thermo_is_reproducible_single_thread(monkeypatch, tmp_path):
    # same inputs + OMP_NUM_THREADS=1 -> bit-identical energy, thread-invariant S
    from ase import Atoms
    from ThermoScreening.thermo.api import xtb_thermo

    monkeypatch.setenv("OMP_NUM_THREADS", "1")
    oh = lambda: Atoms("OH", positions=[[0, 0, 0], [0, 0, 0.97]])
    a = xtb_thermo(oh(), directory=str(tmp_path / "a"), fmax=1e-3)
    b = xtb_thermo(oh(), directory=str(tmp_path / "b"), fmax=1e-3)

    assert abs(a.electronic_energy() - b.electronic_energy()) < 1e-9
    assert abs(a.total_entropy("cal/(mol*K)") - b.total_entropy("cal/(mol*K)")) < 1e-6
