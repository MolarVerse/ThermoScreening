import os
import shutil

import numpy as np
import pytest
from ase import Atoms

from ThermoScreening.exceptions import TSValueError
from ThermoScreening.thermo import conformers

xtb_available = shutil.which("xtb") is not None or "XTB_COMMAND" in os.environ


class _FakeThermo:
    """A stand-in exposing only what a Thermo consumer needs to see it succeeded."""

    def __init__(self, eegtot=-1.0):
        self._eegtot = eegtot

    def total_EeGtot(self):
        return self._eegtot


def test_generate_returns_ase_conformers():
    # n-butane is flexible enough to embed at least one conformer
    result = conformers.generate("CCCC", max_conformers=5)

    assert 1 <= len(result) <= 5
    assert all(isinstance(atoms, Atoms) for atoms in result)
    assert all(atoms.get_chemical_formula() == "C4H10" for atoms in result)


def test_generate_respects_max_conformers():
    result = conformers.generate("CCCCCCCC", max_conformers=3)
    assert len(result) <= 3


def test_generate_energy_window_filters():
    wide = conformers.generate("CCCCCC", max_conformers=15, energy_window=100.0)
    narrow = conformers.generate("CCCCCC", max_conformers=15, energy_window=0.1)
    assert len(narrow) <= len(wide)


def test_generate_without_optimize_runs():
    result = conformers.generate("CCO", max_conformers=3, optimize=False)
    assert len(result) >= 1


def test_generate_rejects_bad_smiles():
    with pytest.raises(ValueError, match="Could not parse SMILES"):
        conformers.generate("this-is-not-smiles")


def test_write_conformers_writes_readable_xyz(tmp_path):
    import ase.io

    result = conformers.generate("CCO", max_conformers=2)
    paths = conformers.write_conformers(result, tmp_path / "confs", prefix="c")

    assert len(paths) == len(result)
    assert all(p.name.startswith("c_") and p.suffix == ".xyz" for p in paths)
    back = ase.io.read(str(paths[0]))
    assert back.get_chemical_formula() == "C2H6O"


def test_public_api_is_exported():
    from ThermoScreening.thermo import generate_conformers, write_conformers, generate_thermo_ensemble

    assert generate_conformers is conformers.generate
    assert write_conformers is conformers.write_conformers
    assert generate_thermo_ensemble is conformers.generate_thermo_ensemble


def test_generate_thermo_ensemble_returns_up_to_max_conformers():
    def thermo_fn(atoms, charge, **kwargs):
        return _FakeThermo()

    result = conformers.generate_thermo_ensemble("CCCCCC", thermo_fn, max_conformers=3)

    assert 1 <= len(result) <= 3
    assert all(isinstance(t, _FakeThermo) for t in result)


def test_generate_thermo_ensemble_forwards_charge_and_kwargs():
    seen = []

    def thermo_fn(atoms, charge, **kwargs):
        seen.append((charge, kwargs))
        return _FakeThermo()

    conformers.generate_thermo_ensemble("CCO", thermo_fn, charge=-1.0, max_conformers=1, solvent="water")

    assert seen == [(-1.0, {"solvent": "water"})]


def test_generate_thermo_ensemble_default_charge_is_zero():
    seen_charges = []

    def thermo_fn(atoms, charge, **kwargs):
        seen_charges.append(charge)
        return _FakeThermo()

    conformers.generate_thermo_ensemble("CCO", thermo_fn, max_conformers=1)

    assert seen_charges == [0.0]


def test_generate_thermo_ensemble_skips_ts_value_error_and_retries():
    call_count = [0]

    def thermo_fn(atoms, charge, **kwargs):
        call_count[0] += 1
        if call_count[0] <= 2:
            raise TSValueError("Imaginary (non-positive) vibrational frequencies are present")
        return _FakeThermo()

    result = conformers.generate_thermo_ensemble("CCCCCC", thermo_fn, max_conformers=2, max_attempts=5)

    assert len(result) == 2
    assert call_count[0] > 2  # the first two failures were skipped and retried past


def test_generate_thermo_ensemble_skips_runtime_error():
    call_count = [0]

    def thermo_fn(atoms, charge, **kwargs):
        call_count[0] += 1
        if call_count[0] == 1:
            raise RuntimeError("xtb failed (exit 1)")
        return _FakeThermo()

    result = conformers.generate_thermo_ensemble("CCCCCC", thermo_fn, max_conformers=1, max_attempts=3)

    assert len(result) == 1


def test_generate_thermo_ensemble_stops_once_target_reached_mid_attempt(monkeypatch):
    # a single attempt's conformer list can be longer than what's still
    # needed (e.g. after an earlier attempt partially succeeded) -- the
    # inner loop must stop as soon as the target is hit, not exhaust the list
    dummy_atoms = [Atoms("H") for _ in range(3)]

    def fake_generate(smiles, max_conformers, prune_rms_thresh, energy_window, random_seed):
        return dummy_atoms

    monkeypatch.setattr(conformers, "generate", fake_generate)

    call_count = [0]

    def thermo_fn(atoms, charge, **kwargs):
        call_count[0] += 1
        return _FakeThermo()

    result = conformers.generate_thermo_ensemble("CCCC", thermo_fn, max_conformers=2, max_attempts=1)

    assert len(result) == 2
    assert call_count[0] == 2  # the 3rd dummy conformer was never processed


def test_generate_thermo_ensemble_raises_when_nothing_succeeds():
    def thermo_fn(atoms, charge, **kwargs):
        raise TSValueError("Imaginary (non-positive) vibrational frequencies are present")

    with pytest.raises(ValueError, match="No conformer"):
        conformers.generate_thermo_ensemble("CCCC", thermo_fn, max_conformers=2, max_attempts=2)


def test_generate_thermo_ensemble_propagates_other_exceptions():
    def thermo_fn(atoms, charge, **kwargs):
        raise KeyError("not a recognised failure mode")

    with pytest.raises(KeyError):
        conformers.generate_thermo_ensemble("CCCC", thermo_fn, max_conformers=1)


@pytest.mark.skipif(not xtb_available, reason="the native xtb binary is not available.")
def test_generate_thermo_ensemble_end_to_end_with_real_xtb(tmp_path):
    # 4-hydroxybutanoic acid's flexible backbone means naive conformer
    # generation occasionally lands a starting geometry on a saddle point
    # (TSValueError) at the anion charge state -- this is exactly the
    # manual retry-by-hand loop this helper replaces.
    from ThermoScreening.thermo.api import xtb_cli_thermo
    from ThermoScreening.thermo.ensemble import EnsembleThermo

    call_count = [0]

    def thermo_fn(atoms, charge, **kwargs):
        call_count[0] += 1
        return xtb_cli_thermo(atoms, charge=charge, directory=str(tmp_path / str(call_count[0])), **kwargs)

    thermos = conformers.generate_thermo_ensemble(
        "OCCCC(=O)[O-]", thermo_fn, charge=-1, max_conformers=3, solvent="water"
    )

    assert 1 <= len(thermos) <= 3
    ensemble = EnsembleThermo(thermos)
    assert ensemble.total_EeGtot() <= min(t.total_EeGtot() for t in thermos)
