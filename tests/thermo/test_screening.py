import csv
import json
from argparse import Namespace
from pathlib import Path

import pytest

from ThermoScreening.exceptions import TSValueError
from ThermoScreening.thermo import screening


class _FakeThermo:
    def electronic_energy(self):
        return -100.0

    def total_energy(self, unit):
        return -10.0

    def total_enthalpy(self, unit):
        return -9.0

    def total_gibbs_free_energy(self, unit):
        return -11.0

    def total_EeGtot(self):
        return -111.0

    def total_entropy(self, unit):
        return 50.0

    def total_heat_capacity(self, unit):
        return 6.0


def _write_xyz(path):
    path.write_text("1\n\nH 0.0 0.0 0.0\n", encoding="utf-8")


def test_load_jobs_from_directory(tmp_path):
    _write_xyz(tmp_path / "mol_b.xyz")
    _write_xyz(tmp_path / "mol_a.xyz")
    (tmp_path / "notes.txt").write_text("ignore me", encoding="utf-8")

    jobs = screening._load_jobs(tmp_path, charge=2.0)

    assert [job.name for job in jobs] == ["mol_a", "mol_b"]  # sorted, .txt ignored
    assert all(job.charge == 2.0 for job in jobs)


def test_load_jobs_from_manifest(tmp_path):
    _write_xyz(tmp_path / "a.xyz")
    _write_xyz(tmp_path / "b.xyz")
    manifest = tmp_path / "m.csv"
    manifest.write_text(
        "name,path,charge\nfirst,a.xyz,-1\n,b.xyz,\n", encoding="utf-8"
    )

    jobs = screening._load_jobs(manifest, charge=0.0)

    assert jobs[0].name == "first" and jobs[0].charge == -1.0
    assert jobs[0].path == tmp_path / "a.xyz"  # resolved relative to the manifest
    assert jobs[1].name == "b"  # blank name falls back to the file stem
    assert jobs[1].charge == 0.0  # blank charge falls back to the default


def test_load_jobs_reads_spin_from_manifest(tmp_path):
    _write_xyz(tmp_path / "a.xyz")
    _write_xyz(tmp_path / "b.xyz")
    manifest = tmp_path / "m.csv"
    manifest.write_text(
        "name,path,charge,spin\nradical,a.xyz,0,0.5\nclosed,b.xyz,0,\n",
        encoding="utf-8",
    )

    jobs = screening._load_jobs(manifest, charge=0.0, spin=None)

    assert jobs[0].spin == 0.5
    assert jobs[1].spin is None  # blank -> default (electron-count guess)


def test_screen_passes_spin_to_dftbplus_thermo(monkeypatch, tmp_path):
    _write_xyz(tmp_path / "mol.xyz")
    manifest = tmp_path / "m.csv"
    manifest.write_text("name,path,charge,spin\nmol,mol.xyz,0,1.0\n", encoding="utf-8")

    captured = {}

    def fake_thermo(atoms, spin=None, **kwargs):
        captured["spin"] = spin
        return _FakeThermo()

    monkeypatch.setattr(screening, "dftbplus_thermo", fake_thermo)
    screening.screen(str(manifest), out=str(tmp_path / "r"), directory=str(tmp_path / "runs"))

    assert captured["spin"] == 1.0


def test_screen_selects_mio_parameter_set(monkeypatch, tmp_path):
    from ThermoScreening.calculator.dftbplus import SPIN_CONSTANTS_MIO

    _write_xyz(tmp_path / "mol.xyz")

    captured = {}

    def fake_thermo(atoms, spin_constants=None, **kwargs):
        captured["spin_constants"] = spin_constants
        captured["kwargs"] = kwargs
        return _FakeThermo()

    monkeypatch.setattr(screening, "dftbplus_thermo", fake_thermo)
    screening.screen(
        str(tmp_path),
        out=str(tmp_path / "r"),
        directory=str(tmp_path / "runs"),
        parameter_set="mio",
    )

    # mio spin constants and DFTB2 Hamiltonian (no third order) flow through
    assert captured["spin_constants"] is SPIN_CONSTANTS_MIO
    assert "Hamiltonian_ThirdOrderFull" not in captured["kwargs"]


def test_screen_rejects_unknown_parameter_set(tmp_path):
    _write_xyz(tmp_path / "mol.xyz")
    with pytest.raises(ValueError, match="Unknown DFTB parameter set"):
        screening.screen(str(tmp_path), out=str(tmp_path / "r"), parameter_set="nope")


def test_screen_passes_solvent_to_dftbplus_thermo(monkeypatch, tmp_path):
    _write_xyz(tmp_path / "mol.xyz")

    captured = {}

    def fake_thermo(atoms, solvent=None, **kwargs):
        captured["solvent"] = solvent
        return _FakeThermo()

    monkeypatch.setattr(screening, "dftbplus_thermo", fake_thermo)
    screening.screen(
        str(tmp_path),
        out=str(tmp_path / "r"),
        directory=str(tmp_path / "runs"),
        solvent="water",
    )

    assert captured["solvent"] == "water"


def test_load_jobs_rejects_unknown_source(tmp_path):
    bad = tmp_path / "thing.txt"
    bad.write_text("x", encoding="utf-8")
    with pytest.raises(TSValueError):
        screening._load_jobs(bad, charge=0.0)


def test_load_jobs_rejects_empty_directory(tmp_path):
    with pytest.raises(TSValueError):
        screening._load_jobs(tmp_path, charge=0.0)


def test_load_jobs_rejects_manifest_without_path_column(tmp_path):
    manifest = tmp_path / "m.csv"
    manifest.write_text("name,charge\nfoo,0\n", encoding="utf-8")
    with pytest.raises(TSValueError, match="path"):
        screening._load_jobs(manifest, charge=0.0)


def test_load_jobs_rejects_empty_manifest(tmp_path):
    manifest = tmp_path / "m.csv"
    manifest.write_text("name,path,charge\n", encoding="utf-8")  # header only
    with pytest.raises(TSValueError, match="no molecules"):
        screening._load_jobs(manifest, charge=0.0)


def test_screen_collects_results_and_isolates_failures(monkeypatch, tmp_path):
    _write_xyz(tmp_path / "good.xyz")
    _write_xyz(tmp_path / "bad.xyz")
    manifest = tmp_path / "m.csv"
    manifest.write_text(
        "name,path,charge\ngood,good.xyz,0\nbad,bad.xyz,-99\n", encoding="utf-8"
    )

    captured = {"dirs": []}

    def fake_thermo(atoms, charge=0, directory=None, **kwargs):
        captured["dirs"].append(directory)
        if charge == -99:
            raise RuntimeError("kaboom")
        return _FakeThermo()

    monkeypatch.setattr(screening, "dftbplus_thermo", fake_thermo)

    out = tmp_path / "out" / "results"
    results = screening.screen(
        str(manifest), out=str(out), directory=str(tmp_path / "runs")
    )

    # one failure does not abort the run
    assert [record["status"] for record in results] == ["ok", "error"]
    assert results[1]["error"] == "kaboom"
    assert results[0]["G_hartree"] == -11.0
    # electronic energy (carries solvation) and absolute Gibbs are reported too
    assert results[0]["Eelec_hartree"] == -100.0
    assert results[0]["G_total_hartree"] == -111.0

    # each molecule ran in its own directory
    assert str(tmp_path / "runs" / "good") in captured["dirs"]

    csv_rows = list(csv.DictReader(out.with_suffix(".csv").open(encoding="utf-8")))
    assert csv_rows[0]["name"] == "good"
    assert csv_rows[1]["status"] == "error"
    assert set(screening._RESULT_FIELDS) <= set(csv_rows[0].keys())

    data = json.loads(out.with_suffix(".json").read_text(encoding="utf-8"))
    assert data[1]["error"] == "kaboom"


def test_screen_uses_default_charge_for_directory_input(monkeypatch, tmp_path):
    _write_xyz(tmp_path / "m1.xyz")
    charges = []

    def fake_thermo(atoms, charge=0, directory=None, **kwargs):
        charges.append(charge)
        return _FakeThermo()

    monkeypatch.setattr(screening, "dftbplus_thermo", fake_thermo)

    screening.screen(
        str(tmp_path), out=str(tmp_path / "r"), charge=-2.0,
        directory=str(tmp_path / "runs"),
    )

    assert charges == [-2.0]


def test_cli_parse_args_routes_screen():
    import ThermoScreening.cli.thermo as cli

    args = cli.parse_args(
        ["screen", "molecules.csv", "-o", "out", "--charge", "-1", "--temperature", "300"]
    )

    assert args.command == "screen"
    assert args.source == "molecules.csv"
    assert args.out == "out"
    assert args.charge == -1.0
    assert args.temperature == 300.0
    assert args.parameter_set == "3ob"  # default set
    assert args.solvent is None  # gas phase by default

    mio_args = cli.parse_args(["screen", "molecules.csv", "--parameter-set", "mio"])
    assert mio_args.parameter_set == "mio"

    solv_args = cli.parse_args(["screen", "molecules.csv", "--solvent", "water"])
    assert solv_args.solvent == "water"


def test_cli_run_screen_returns_failure_count(monkeypatch):
    import ThermoScreening.cli.thermo as cli

    monkeypatch.setattr(
        cli, "screen", lambda *args, **kwargs: [{"status": "ok"}, {"status": "error"}]
    )

    args = Namespace(
        source="x", out="res", charge=0.0, temperature=298.15,
        pressure=101325.0, directory="screening", parameter_set="3ob",
        solvent=None,
    )

    assert cli.run_screen(args) == 1  # one molecule failed
