import csv
import json
from argparse import Namespace
from concurrent.futures import Future

import pytest

from ThermoScreening.exceptions import TSValueError
from ThermoScreening.thermo import screening
from ThermoScreening.thermo._units import HARTREE_TO_EV
from ThermoScreening.thermo.reactions import SHE_ABSOLUTE_POTENTIAL


class _SyncExecutor:
    def __init__(self, max_workers=None):
        self.max_workers = max_workers

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def submit(self, function, *args, **kwargs):
        future = Future()
        future.set_result(function(*args, **kwargs))
        return future


class _ChargeThermo:
    def __init__(self, charge):
        self.charge = charge

    def electronic_energy(self):
        return -100.0 + self.charge

    def total_energy(self, unit):
        return -10.0

    def total_enthalpy(self, unit):
        return -9.0

    def total_gibbs_free_energy(self, unit):
        return -11.0

    def total_EeGtot(self):
        return -100.0 + 0.1 * self.charge

    def total_entropy(self, unit):
        return 50.0

    def total_heat_capacity(self, unit):
        return 6.0


def _write_xyz(path):
    path.write_text("1\n\nH 0.0 0.0 0.0\n", encoding="utf-8")


def _fake_state_screen(energies, captured=None, failed=None):
    failed = failed or {}

    def fake_screen(source, **kwargs):
        if captured is not None:
            captured.update(source=source, kwargs=kwargs)
        with open(source, newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        results = []
        for row in rows:
            molecule, state = row["name"].rsplit("--", 1)
            error = failed.get((molecule, state))
            if error:
                results.append({"name": row["name"], "status": "error", "error": error})
            else:
                results.append(
                    {
                        "name": row["name"],
                        "status": "ok",
                        "G_total_hartree": energies[molecule][state],
                    }
                )
        return results

    return fake_screen


def test_redox_screen_calculates_all_states_in_parallel(monkeypatch, tmp_path):
    _write_xyz(tmp_path / "a.xyz")
    _write_xyz(tmp_path / "b.xyz")
    charges = []

    def fake_thermo(atoms, charge=0.0, **kwargs):
        charges.append(charge)
        return _ChargeThermo(charge)

    monkeypatch.setattr(screening, "dftbplus_thermo", fake_thermo)
    monkeypatch.setattr(screening, "ProcessPoolExecutor", _SyncExecutor)

    results = screening.redox_screen(
        tmp_path,
        out=tmp_path / "redox",
        directory=tmp_path / "work",
        jobs=4,
    )

    assert [result["name"] for result in results] == ["a", "b"]
    assert sorted(charges) == [-2.0, -2.0, -1.0, -1.0, 0.0, 0.0]
    assert all(result["status"] == "ok" for result in results)
    assert all(result["potential_scale"] == "SHE" for result in results)
    assert results[0]["E1_V"] == pytest.approx(0.1 * HARTREE_TO_EV - SHE_ABSOLUTE_POTENTIAL)
    assert (tmp_path / "redox-states.csv").is_file()
    assert (tmp_path / "redox.csv").is_file()


def test_redox_screen_resumes_all_completed_states(monkeypatch, tmp_path):
    _write_xyz(tmp_path / "molecule.xyz")
    calls = []

    def fake_thermo(atoms, charge=0.0, **kwargs):
        calls.append(charge)
        return _ChargeThermo(charge)

    monkeypatch.setattr(screening, "dftbplus_thermo", fake_thermo)
    arguments = {
        "source": tmp_path / "molecule.xyz",
        "out": tmp_path / "redox",
        "directory": tmp_path / "work",
    }

    first = screening.redox_screen(**arguments)
    second = screening.redox_screen(**arguments, resume=True)

    assert calls == [0.0, -1.0, -2.0]
    assert first == second


def test_redox_screen_calibrates_each_reduction_separately(monkeypatch, tmp_path):
    _write_xyz(tmp_path / "target.xyz")
    _write_xyz(tmp_path / "reference.xyz")
    energies = {
        "target": {
            "oxidized": -200.0,
            "reduced_once": -200.11,
            "reduced_twice": -200.20,
        },
        "reference": {
            "oxidized": -100.0,
            "reduced_once": -100.10,
            "reduced_twice": -100.18,
        },
    }
    monkeypatch.setattr(screening, "screen", _fake_state_screen(energies))

    result = screening.redox_screen(
        tmp_path / "target.xyz",
        out=tmp_path / "out",
        directory=tmp_path / "work",
        reference=tmp_path / "reference.xyz",
        reference_e1=-0.75,
        reference_e2=-1.40,
        potential_scale="Fc/Fc+",
    )[0]

    assert result["E1_V"] == pytest.approx(-0.75 + 0.01 * HARTREE_TO_EV)
    assert result["E2_V"] == pytest.approx(-1.40 + 0.01 * HARTREE_TO_EV)
    assert result["E2e_V"] == pytest.approx((result["E1_V"] + result["E2_V"]) / 2)
    assert result["potential_scale"] == "Fc/Fc+"


def test_redox_screen_reuses_reference_from_candidate_set(monkeypatch, tmp_path):
    _write_xyz(tmp_path / "target.xyz")
    _write_xyz(tmp_path / "reference.xyz")
    energies = {
        "target": {
            "oxidized": -200.0,
            "reduced_once": -200.11,
            "reduced_twice": -200.20,
        },
        "reference": {
            "oxidized": -100.0,
            "reduced_once": -100.10,
            "reduced_twice": -100.18,
        },
    }
    captured = {}
    monkeypatch.setattr(
        screening, "screen", _fake_state_screen(energies, captured=captured)
    )

    results = screening.redox_screen(
        tmp_path,
        out=tmp_path / "out",
        directory=tmp_path / "work",
        reference="reference",
        reference_e1=-0.75,
        reference_e2=-1.40,
    )

    with open(captured["source"], newline="", encoding="utf-8") as handle:
        state_rows = list(csv.DictReader(handle))
    assert len(state_rows) == 6
    reference = next(result for result in results if result["name"] == "reference")
    assert reference["E1_V"] == pytest.approx(-0.75)
    assert reference["E2_V"] == pytest.approx(-1.40)


def test_redox_screen_rejects_candidate_reference_charge_mismatch(tmp_path):
    _write_xyz(tmp_path / "reference.xyz")

    with pytest.raises(TSValueError, match="does not match candidate"):
        screening.redox_screen(
            tmp_path,
            out=tmp_path / "out",
            directory=tmp_path / "work",
            reference="reference",
            reference_e1=-0.75,
            reference_e2=-1.40,
            reference_charge=1,
        )


def test_redox_screen_reuses_reference_from_same_structure_path(monkeypatch, tmp_path):
    _write_xyz(tmp_path / "reference.xyz")
    energies = {
        "reference": {
            "oxidized": -100.0,
            "reduced_once": -100.10,
            "reduced_twice": -100.18,
        }
    }
    captured = {}
    monkeypatch.setattr(
        screening, "screen", _fake_state_screen(energies, captured=captured)
    )

    result = screening.redox_screen(
        tmp_path / "reference.xyz",
        out=tmp_path / "out",
        directory=tmp_path / "work",
        reference=tmp_path / "reference.xyz",
        reference_e1=-0.75,
        reference_e2=-1.40,
    )[0]

    with open(captured["source"], newline="", encoding="utf-8") as handle:
        assert len(list(csv.DictReader(handle))) == 3
    assert result["E1_V"] == pytest.approx(-0.75)


def test_redox_screen_renames_a_distinct_reference_that_collides_with_candidate(
    monkeypatch, tmp_path
):
    candidate_directory = tmp_path / "candidates"
    reference_directory = tmp_path / "references"
    candidate_directory.mkdir()
    reference_directory.mkdir()
    _write_xyz(candidate_directory / "same.xyz")
    _write_xyz(reference_directory / "same.xyz")
    energies = {
        "same": {
            "oxidized": -200.0,
            "reduced_once": -200.11,
            "reduced_twice": -200.20,
        },
        "reference": {
            "oxidized": -100.0,
            "reduced_once": -100.10,
            "reduced_twice": -100.18,
        },
    }
    captured = {}
    monkeypatch.setattr(
        screening, "screen", _fake_state_screen(energies, captured=captured)
    )

    result = screening.redox_screen(
        candidate_directory / "same.xyz",
        out=tmp_path / "out",
        directory=tmp_path / "work",
        reference=reference_directory / "same.xyz",
        reference_e1=-0.75,
        reference_e2=-1.40,
    )[0]

    with open(captured["source"], newline="", encoding="utf-8") as handle:
        names = [row["name"] for row in csv.DictReader(handle)]
    assert names == [
        "same--oxidized",
        "same--reduced_once",
        "same--reduced_twice",
        "reference--oxidized",
        "reference--reduced_once",
        "reference--reduced_twice",
    ]
    assert result["status"] == "ok"


def test_redox_screen_propagates_state_and_reference_failures(monkeypatch, tmp_path):
    _write_xyz(tmp_path / "target.xyz")
    _write_xyz(tmp_path / "reference.xyz")
    energies = {
        "target": {
            "oxidized": -200.0,
            "reduced_once": -200.1,
            "reduced_twice": -200.2,
        },
        "reference": {
            "oxidized": -100.0,
            "reduced_once": -100.1,
            "reduced_twice": -100.2,
        },
    }
    monkeypatch.setattr(
        screening,
        "screen",
        _fake_state_screen(
            energies, failed={("reference", "reduced_once"): "no minimum"}
        ),
    )

    result = screening.redox_screen(
        tmp_path / "target.xyz",
        out=tmp_path / "out",
        directory=tmp_path / "work",
        reference=tmp_path / "reference.xyz",
        reference_e1=-0.75,
        reference_e2=-1.40,
    )[0]

    assert result["status"] == "error"
    assert "reference: reduced_once: no minimum" in result["error"]
    assert "E1_V" not in result


def test_redox_screen_reports_missing_state_results_and_energies(monkeypatch, tmp_path):
    _write_xyz(tmp_path / "molecule.xyz")

    def incomplete_screen(source, **kwargs):
        return [
            {
                "name": "molecule--oxidized",
                "status": "ok",
                "G_total_hartree": -100.0,
            },
            {"name": "molecule--reduced_once", "status": "ok"},
        ]

    monkeypatch.setattr(screening, "screen", incomplete_screen)

    result = screening.redox_screen(
        tmp_path / "molecule.xyz",
        out=tmp_path / "out",
        directory=tmp_path / "work",
    )[0]

    assert result["status"] == "error"
    assert "reduced_once: total Gibbs free energy is missing" in result["error"]
    assert "reduced_twice: state calculation produced no result" in result["error"]


@pytest.mark.parametrize("invalid_energy", [float("nan"), float("inf"), "invalid"])
def test_redox_screen_rejects_nonfinite_state_energies(
    monkeypatch, tmp_path, invalid_energy
):
    _write_xyz(tmp_path / "molecule.xyz")
    energies = {
        "molecule": {
            "oxidized": -100.0,
            "reduced_once": invalid_energy,
            "reduced_twice": -100.2,
        }
    }
    monkeypatch.setattr(screening, "screen", _fake_state_screen(energies))

    result = screening.redox_screen(
        tmp_path / "molecule.xyz",
        out=tmp_path / "out",
        directory=tmp_path / "work",
    )[0]

    assert result["status"] == "error"
    assert "reduced_once: total Gibbs free energy is not finite" in result["error"]
    assert "G_reduced_once_hartree" not in result


def test_redox_screen_flags_potential_inversion(monkeypatch, tmp_path):
    _write_xyz(tmp_path / "molecule.xyz")
    energies = {
        "molecule": {
            "oxidized": -100.0,
            "reduced_once": -100.05,
            "reduced_twice": -100.12,
        }
    }
    monkeypatch.setattr(screening, "screen", _fake_state_screen(energies))

    result = screening.redox_screen(
        tmp_path / "molecule.xyz",
        out=tmp_path / "out",
        directory=tmp_path / "work",
    )[0]

    assert result["potential_inversion"] is True
    assert result["potential_gap_V"] < 0


def test_redox_screen_accepts_smiles_and_writes_generated_input(monkeypatch, tmp_path):
    energies = {
        "molecule": {
            "oxidized": -100.0,
            "reduced_once": -100.1,
            "reduced_twice": -100.2,
        }
    }
    monkeypatch.setattr(screening, "screen", _fake_state_screen(energies))

    result = screening.redox_screen(
        "CCO",
        out=tmp_path / "out",
        directory=tmp_path / "work",
        max_conformers=2,
    )[0]

    assert result["status"] == "ok"
    assert (tmp_path / "work" / "inputs" / "molecule_0.xyz").is_file()


def test_redox_screen_infers_formal_charge_from_ionic_smiles(monkeypatch, tmp_path):
    from ase import Atoms

    atoms = Atoms("H", positions=[[0.0, 0.0, 0.0]])
    atoms.info["formal_charge"] = -1
    energies = {
        "molecule": {
            "oxidized": -100.0,
            "reduced_once": -100.1,
            "reduced_twice": -100.2,
        }
    }
    captured = {}
    monkeypatch.setattr(screening, "generate_conformers", lambda *args, **kwargs: [atoms])
    monkeypatch.setattr(
        screening, "screen", _fake_state_screen(energies, captured=captured)
    )

    result = screening.redox_screen(
        "[H-]",
        out=tmp_path / "out",
        directory=tmp_path / "work",
    )[0]

    with open(captured["source"], newline="", encoding="utf-8") as handle:
        charges = [int(row["charge"]) for row in csv.DictReader(handle)]
    assert result["charge"] == -1
    assert charges == [-1, -2, -3]


def test_redox_manifest_supports_paths_smiles_charges_and_spins(monkeypatch, tmp_path):
    _write_xyz(tmp_path / "first.xyz")
    manifest = tmp_path / "molecules.csv"
    manifest.write_text(
        "name,path,smiles,charge,spin,oxidized_spin,reduced_once_spin,"
        "reduced_twice_spin\n"
        "first,first.xyz,,-1,0,,0.5,0\n"
        "second,,CCO,0,,1,1.5,2\n",
        encoding="utf-8",
    )
    energies = {
        "first": {
            "oxidized": -100.0,
            "reduced_once": -100.1,
            "reduced_twice": -100.2,
        },
        "second": {
            "oxidized": -200.0,
            "reduced_once": -200.1,
            "reduced_twice": -200.2,
        },
    }
    captured = {}
    monkeypatch.setattr(
        screening, "screen", _fake_state_screen(energies, captured=captured)
    )

    results = screening.redox_screen(
        manifest,
        out=tmp_path / "out",
        directory=tmp_path / "work",
        max_conformers=2,
        jobs=3,
    )

    with open(captured["source"], newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert [result["charge"] for result in results] == [-1, 0]
    assert [row["charge"] for row in rows[:3]] == ["-1", "-2", "-3"]
    assert [row["spin"] for row in rows[:3]] == ["0.0", "0.5", "0.0"]
    assert captured["kwargs"]["jobs"] == 3


@pytest.mark.parametrize(
    ("manifest_text", "message"),
    [
        ("name,path,smiles\na,,\n", "exactly one"),
        ("name,path,smiles\na,a.xyz,CCO\n", "exactly one"),
        ("name,path,smiles\na,,CCO\na,,CCC\n", "Duplicate molecule"),
        ("name,path,charge\na,a.xyz,bad\n", "charge in row 2 must be numeric"),
    ],
)
def test_redox_manifest_rejects_ambiguous_or_invalid_rows(
    monkeypatch, tmp_path, manifest_text, message
):
    _write_xyz(tmp_path / "a.xyz")
    manifest = tmp_path / "molecules.csv"
    manifest.write_text(manifest_text, encoding="utf-8")
    monkeypatch.setattr(screening, "generate_conformers", lambda *args, **kwargs: [object()])
    monkeypatch.setattr(
        screening,
        "write_conformers",
        lambda conformers, directory, prefix: [tmp_path / f"{prefix}.xyz"],
    )

    with pytest.raises(TSValueError, match=message):
        screening.redox_screen(
            manifest,
            out=tmp_path / "out",
            directory=tmp_path / "work",
        )


def test_redox_screen_rejects_missing_empty_and_unsupported_inputs(tmp_path):
    empty = tmp_path / "empty"
    empty.mkdir()
    unsupported = tmp_path / "molecule.txt"
    unsupported.write_text("not a structure", encoding="utf-8")
    invalid_name = tmp_path / "invalid-name.csv"
    invalid_name.write_text("name,path\n../escape,molecule.xyz\n", encoding="utf-8")
    _write_xyz(tmp_path / "molecule.xyz")

    with pytest.raises(TSValueError, match="input not found"):
        screening.redox_screen(
            tmp_path / "missing.xyz", out=tmp_path / "out", directory=tmp_path / "work"
        )
    with pytest.raises(TSValueError, match="No molecules found"):
        screening.redox_screen(empty, out=tmp_path / "out", directory=tmp_path / "work")
    with pytest.raises(TSValueError, match="must be .xyz or .gen"):
        screening.redox_screen(
            unsupported, out=tmp_path / "out", directory=tmp_path / "work"
        )
    with pytest.raises(TSValueError, match="Invalid molecule name"):
        screening.redox_screen(
            invalid_name, out=tmp_path / "out", directory=tmp_path / "work"
        )


def test_redox_screen_requires_one_reference_molecule(tmp_path):
    _write_xyz(tmp_path / "target.xyz")
    references = tmp_path / "references"
    references.mkdir()
    _write_xyz(references / "one.xyz")
    _write_xyz(references / "two.xyz")

    with pytest.raises(TSValueError, match="exactly one molecule"):
        screening.redox_screen(
            tmp_path / "target.xyz",
            out=tmp_path / "out",
            directory=tmp_path / "work",
            reference=references,
            reference_e1=-0.75,
            reference_e2=-1.40,
        )


def test_redox_screen_reports_smiles_generation_failure(monkeypatch, tmp_path):
    monkeypatch.setattr(
        screening,
        "generate_conformers",
        lambda *args, **kwargs: (_ for _ in ()).throw(ValueError("embedding failed")),
    )

    with pytest.raises(TSValueError, match="Could not generate.*embedding failed"):
        screening.redox_screen(
            "invalid-smiles",
            out=tmp_path / "out",
            directory=tmp_path / "work",
        )


def test_redox_screen_reports_empty_smiles_generation(monkeypatch, tmp_path):
    monkeypatch.setattr(screening, "generate_conformers", lambda *args, **kwargs: [])

    with pytest.raises(TSValueError, match="no conformers returned"):
        screening.redox_screen(
            "CCO",
            out=tmp_path / "out",
            directory=tmp_path / "work",
        )


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"reference": "CCO"}, "supplied together"),
        ({"reference_e1": -0.7}, "supplied together"),
        ({"max_conformers": 0}, "integer >= 1"),
        ({"max_conformers": 1.5}, "integer >= 1"),
        ({"charge": float("inf")}, "charge must be finite"),
        ({"charge": -0.5}, "charge must be an integer"),
        ({"spin": 0.25}, "spin must be a non-negative integer or half-integer"),
        ({"spin": 0.0}, "oxidized_spin=0 is incompatible"),
        ({"potential_scale": "Fc/Fc+"}, "requires reference calibration"),
        ({"engine": "xtb", "solvent": "water"}, "does not support solvent"),
        ({"engine": "xtb", "dispersion": "d3-bj"}, "only by the DFTB"),
        ({"engine": "xtb", "parameter_set": "3ob"}, "only to the DFTB"),
        ({"engine": "dftb+", "method": "GFN2-xTB"}, "only to xTB"),
    ],
)
def test_redox_screen_rejects_invalid_workflow_options(tmp_path, kwargs, message):
    _write_xyz(tmp_path / "molecule.xyz")
    with pytest.raises(TSValueError, match=message):
        screening.redox_screen(
            tmp_path / "molecule.xyz",
            out=tmp_path / "out",
            directory=tmp_path / "work",
            **kwargs,
        )


def test_redox_screen_writes_machine_readable_aggregate(monkeypatch, tmp_path):
    _write_xyz(tmp_path / "molecule.xyz")
    energies = {
        "molecule": {
            "oxidized": -100.0,
            "reduced_once": -100.1,
            "reduced_twice": -100.18,
        }
    }
    monkeypatch.setattr(screening, "screen", _fake_state_screen(energies))

    results = screening.redox_screen(
        tmp_path / "molecule.xyz",
        out=tmp_path / "nested" / "redox",
        directory=tmp_path / "work",
    )

    with open(tmp_path / "nested" / "redox.json", encoding="utf-8") as handle:
        assert json.load(handle) == results
    with open(tmp_path / "nested" / "redox.csv", newline="", encoding="utf-8") as handle:
        row = next(csv.DictReader(handle))
    assert row["name"] == "molecule"
    assert row["potential_inversion"] == "False"
    metadata = json.loads(
        (tmp_path / "nested" / "redox-run.json").read_text(encoding="utf-8")
    )
    assert metadata["single_starting_geometry_approximation"] is True
    assert metadata["candidates"][0]["structure_sha256"]
    assert metadata["run_fingerprint"] == results[0]["run_fingerprint"]


def test_cli_parses_and_runs_redox(monkeypatch, capsys):
    import ThermoScreening.cli.thermo as cli

    args = cli.parse_args(
        [
            "redox",
            "molecules.csv",
            "--reference",
            "AQ",
            "--reference-e1",
            "-0.75",
            "--reference-e2",
            "-1.40",
            "--potential-scale",
            "Fc/Fc+",
            "--jobs",
            "6",
            "--resume",
        ]
    )
    captured = {}

    def fake_redox(source, **kwargs):
        captured.update(source=source, kwargs=kwargs)
        return [
            {
                "name": "candidate",
                "status": "ok",
                "E1_V": -0.6,
                "E2_V": -1.2,
                "E2e_V": -0.9,
                "potential_scale": "Fc/Fc+",
            }
        ]

    monkeypatch.setattr(cli, "redox_screen", fake_redox)

    assert args.command == "redox"
    assert args.charge is None
    assert args.parameter_set is None
    assert args.method is None
    assert cli.run_redox(args) == 0
    assert captured["kwargs"]["jobs"] == 6
    assert captured["kwargs"]["resume"] is True
    assert captured["kwargs"]["reference_e1"] == -0.75
    assert "E1=-0.6000 V" in capsys.readouterr().out


def test_cli_redox_reports_workflow_errors(monkeypatch, capsys):
    import ThermoScreening.cli.thermo as cli

    args = cli.parse_args(["redox", "missing.xyz"])
    monkeypatch.setattr(
        cli,
        "redox_screen",
        lambda *args, **kwargs: (_ for _ in ()).throw(TSValueError("missing")),
    )

    assert cli.run_redox(args) == 1
    assert "Redox screen failed: missing" in capsys.readouterr().err


def test_cli_redox_reports_per_molecule_failures(monkeypatch, capsys):
    import ThermoScreening.cli.thermo as cli

    args = cli.parse_args(["redox", "molecules.csv"])
    monkeypatch.setattr(
        cli,
        "redox_screen",
        lambda *args, **kwargs: [
            {"name": "failed", "status": "error", "error": "anion failed"}
        ],
    )

    assert cli.run_redox(args) == 1
    output = capsys.readouterr().out
    assert "Failed (1)" in output
    assert "failed: anion failed" in output
