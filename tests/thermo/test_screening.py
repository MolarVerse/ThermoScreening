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


def test_screen_passes_dispersion_to_dftbplus_thermo(monkeypatch, tmp_path):
    _write_xyz(tmp_path / "mol.xyz")

    captured = {}

    def fake_thermo(atoms, dispersion=None, **kwargs):
        captured["dispersion"] = dispersion
        return _FakeThermo()

    monkeypatch.setattr(screening, "dftbplus_thermo", fake_thermo)
    screening.screen(
        str(tmp_path),
        out=str(tmp_path / "r"),
        directory=str(tmp_path / "runs"),
        dispersion="d3-bj",
    )

    assert captured["dispersion"] == "d3-bj"


def test_screen_dispatches_to_xtb_engine(monkeypatch, tmp_path):
    _write_xyz(tmp_path / "mol.xyz")

    captured = {}

    def fake_xtb(atoms, method="GFN2-xTB", **kwargs):
        captured["method"] = method
        captured["called"] = "xtb"
        return _FakeThermo()

    def fake_dftb(*args, **kwargs):
        captured["called"] = "dftb+"
        return _FakeThermo()

    monkeypatch.setattr(screening, "xtb_thermo", fake_xtb)
    monkeypatch.setattr(screening, "dftbplus_thermo", fake_dftb)
    screening.screen(
        str(tmp_path),
        out=str(tmp_path / "r"),
        directory=str(tmp_path / "runs"),
        engine="xtb",
        method="GFN1-xTB",
    )

    assert captured["called"] == "xtb"
    assert captured["method"] == "GFN1-xTB"


def test_screen_resume_skips_completed_and_reruns_failed(monkeypatch, tmp_path):
    from pathlib import Path

    _write_xyz(tmp_path / "mol_a.xyz")
    _write_xyz(tmp_path / "mol_b.xyz")
    out = tmp_path / "out"

    calls = []

    def thermo_fail_b(atoms, directory=None, **kwargs):
        name = Path(directory).name
        calls.append(name)
        if name == "mol_b":
            raise RuntimeError("boom")
        return _FakeThermo()

    # first run: mol_a ok, mol_b fails
    monkeypatch.setattr(screening, "dftbplus_thermo", thermo_fail_b)
    first = screening.screen(str(tmp_path), out=str(out), directory=str(tmp_path / "r1"))
    assert {r["name"]: r["status"] for r in first} == {"mol_a": "ok", "mol_b": "error"}

    # resume: mol_a (ok) skipped, only mol_b (was error) re-run
    calls.clear()
    monkeypatch.setattr(
        screening, "dftbplus_thermo",
        lambda atoms, directory=None, **kw: (calls.append(Path(directory).name), _FakeThermo())[1],
    )
    second = screening.screen(
        str(tmp_path), out=str(out), directory=str(tmp_path / "r2"), resume=True
    )
    assert calls == ["mol_b"]  # only the previously-failed molecule re-run
    assert {r["name"]: r["status"] for r in second} == {"mol_a": "ok", "mol_b": "ok"}


def test_load_completed_handles_missing_and_corrupt(tmp_path):
    # no prior file -> nothing to resume
    assert screening._load_completed(str(tmp_path / "nope")) == {}
    # unreadable/corrupt json -> empty, does not crash
    (tmp_path / "bad.json").write_text("{ not valid json", encoding="utf-8")
    assert screening._load_completed(str(tmp_path / "bad")) == {}
    # valid JSON but not a list of records -> empty, does not crash
    (tmp_path / "scalar.json").write_text("42", encoding="utf-8")
    assert screening._load_completed(str(tmp_path / "scalar")) == {}


def test_screen_writes_results_incrementally(monkeypatch, tmp_path):
    for name in ("mol_a", "mol_b", "mol_c"):
        _write_xyz(tmp_path / f"{name}.xyz")

    writes = []
    real_write = screening._write_results

    def spy_write(results, out):
        writes.append(len(results))
        return real_write(results, out)

    monkeypatch.setattr(screening, "_write_results", spy_write)
    monkeypatch.setattr(screening, "dftbplus_thermo", lambda atoms, **kw: _FakeThermo())

    screening.screen(str(tmp_path), out=str(tmp_path / "out"), directory=str(tmp_path / "r"))
    # written after each molecule with a growing result set (durable/resumable)
    assert writes == [1, 2, 3]


def test_screen_dispatches_to_xtb_cli_engine(monkeypatch, tmp_path):
    _write_xyz(tmp_path / "mol.xyz")

    captured = {}

    def fake_xtb_cli(atoms, method="GFN2-xTB", solvent=None, **kwargs):
        captured.update(called="xtb-cli", method=method, solvent=solvent)
        return _FakeThermo()

    monkeypatch.setattr(screening, "xtb_cli_thermo", fake_xtb_cli)
    screening.screen(
        str(tmp_path),
        out=str(tmp_path / "r"),
        directory=str(tmp_path / "runs"),
        engine="xtb-cli",
        solvent="water",
    )

    assert captured["called"] == "xtb-cli"
    assert captured["solvent"] == "water"  # xtb-cli supports solvation


def test_screen_rejects_unknown_engine(tmp_path):
    _write_xyz(tmp_path / "mol.xyz")
    with pytest.raises(TSValueError, match="Unknown engine"):
        screening.screen(str(tmp_path), out=str(tmp_path / "r"), engine="orca")


def test_screen_passes_quasi_rrho_to_dftbplus_thermo(monkeypatch, tmp_path):
    _write_xyz(tmp_path / "mol.xyz")

    captured = {}

    def fake_thermo(atoms, quasi_rrho=False, **kwargs):
        captured["quasi_rrho"] = quasi_rrho
        return _FakeThermo()

    monkeypatch.setattr(screening, "dftbplus_thermo", fake_thermo)
    screening.screen(
        str(tmp_path),
        out=str(tmp_path / "r"),
        directory=str(tmp_path / "runs"),
        quasi_rrho=True,
    )

    assert captured["quasi_rrho"] is True


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
    assert solv_args.quasi_rrho is False  # harmonic by default

    qrrho_args = cli.parse_args(["screen", "molecules.csv", "--quasi-rrho"])
    assert qrrho_args.quasi_rrho is True
    assert args.engine == "dftb+"  # default engine
    assert args.method == "GFN2-xTB"

    xtb_args = cli.parse_args(
        ["screen", "molecules.csv", "--engine", "xtb", "--method", "GFN1-xTB"]
    )
    assert xtb_args.engine == "xtb"
    assert xtb_args.method == "GFN1-xTB"

    cli_args = cli.parse_args(["screen", "molecules.csv", "--engine", "xtb-cli"])
    assert cli_args.engine == "xtb-cli"

    assert args.resume is False  # default
    assert cli.parse_args(["screen", "molecules.csv", "--resume"]).resume is True


def test_cli_run_screen_returns_failure_count(monkeypatch):
    import ThermoScreening.cli.thermo as cli

    monkeypatch.setattr(
        cli, "screen", lambda *args, **kwargs: [{"status": "ok"}, {"status": "error"}]
    )

    args = Namespace(
        source="x", out="res", charge=0.0, temperature=298.15,
        pressure=101325.0, directory="screening", parameter_set="3ob",
        solvent=None, dispersion=None, quasi_rrho=False, engine="dftb+",
        method="GFN2-xTB", resume=False,
    )

    assert cli.run_screen(args) == 1  # one molecule failed


def test_cli_parse_args_routes_conformers():
    import ThermoScreening.cli.thermo as cli

    args = cli.parse_args(
        ["conformers", "CCCC", "-o", "confs", "--max-conformers", "5",
         "--energy-window", "2.0", "--no-optimize"]
    )
    assert args.command == "conformers"
    assert args.smiles == "CCCC"
    assert args.out_dir == "confs"
    assert args.max_conformers == 5
    assert args.energy_window == 2.0
    assert args.no_optimize is True


def test_cli_run_conformers(monkeypatch, tmp_path):
    import ThermoScreening.cli.thermo as cli

    captured = {}

    def fake_generate(smiles, max_conformers=10, optimize=True, energy_window=None):
        captured.update(smiles=smiles, max_conformers=max_conformers,
                        optimize=optimize, energy_window=energy_window)
        return ["conf_a", "conf_b"]

    def fake_write(conformers, out_dir, prefix="conformer"):
        captured["written"] = len(conformers)
        return [tmp_path / "a.xyz", tmp_path / "b.xyz"]

    monkeypatch.setattr(cli, "generate_conformers", fake_generate)
    monkeypatch.setattr(cli, "write_conformers", fake_write)

    args = Namespace(smiles="CCO", out_dir=str(tmp_path / "out"),
                     max_conformers=7, energy_window=1.5, no_optimize=True)

    assert cli.run_conformers(args) == 0
    assert captured["smiles"] == "CCO"
    assert captured["max_conformers"] == 7
    assert captured["optimize"] is False  # --no-optimize
    assert captured["energy_window"] == 1.5
    assert captured["written"] == 2


def test_cli_run_conformers_reports_bad_smiles(monkeypatch, capsys):
    import ThermoScreening.cli.thermo as cli

    def fail(*args, **kwargs):
        raise ValueError("Could not parse SMILES: 'oops'")

    monkeypatch.setattr(cli, "generate_conformers", fail)
    args = Namespace(smiles="oops", out_dir="out", max_conformers=10,
                     energy_window=None, no_optimize=False)

    # clean exit code + message on stderr, not a traceback
    assert cli.run_conformers(args) == 1
    assert "Conformer generation failed" in capsys.readouterr().err
