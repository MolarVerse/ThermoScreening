"""Scientific invariants for the anthraquinone redox workflow."""

import csv
import json
from pathlib import Path

import ase.io
import pytest

from ThermoScreening.thermo import screening
from ThermoScreening.thermo._units import HARTREE_TO_EV


_DATA = (
    Path(__file__).resolve().parents[1]
    / "data"
    / "regression"
    / "anthraquinones.csv"
)
_ENERGIES = {
    "anthraquinone": {
        "oxidized": -100.0,
        "reduced_once": -100.10,
        "reduced_twice": -100.18,
    },
    "1-hydroxyanthraquinone": {
        "oxidized": -110.0,
        "reduced_once": -110.115,
        "reduced_twice": -110.205,
    },
    "2-hydroxyanthraquinone": {
        "oxidized": -120.0,
        "reduced_once": -120.08,
        "reduced_twice": -120.20,
    },
}


def _state_screen(source, **kwargs):
    with open(source, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    records = []
    for row in rows:
        molecule, state = row["name"].rsplit("--", 1)
        records.append(
            {
                "name": row["name"],
                "status": "ok",
                "G_total_hartree": _ENERGIES[molecule][state],
            }
        )
    screening._write_results(records, kwargs["out"])
    return records


def _run_regression(monkeypatch, tmp_path, name):
    monkeypatch.setattr(screening, "screen", _state_screen)
    out = tmp_path / name
    results = screening.redox_screen(
        _DATA,
        out=out,
        directory=tmp_path / f"{name}-work",
        reference="anthraquinone",
        reference_e1=-0.75,
        reference_e2=-1.40,
        potential_scale="AQ reference",
        max_conformers=3,
    )
    metadata = json.loads(
        out.with_name(f"{out.name}-run.json").read_text(encoding="utf-8")
    )
    return results, metadata


def test_anthraquinone_workflow_regression(monkeypatch, tmp_path):
    results, metadata = _run_regression(monkeypatch, tmp_path, "first")
    by_name = {result["name"]: result for result in results}

    assert list(by_name) == [
        "anthraquinone",
        "1-hydroxyanthraquinone",
        "2-hydroxyanthraquinone",
    ]
    assert all(result["status"] == "ok" for result in results)
    assert by_name["anthraquinone"]["E1_V"] == pytest.approx(-0.75)
    assert by_name["anthraquinone"]["E2_V"] == pytest.approx(-1.40)
    assert by_name["1-hydroxyanthraquinone"]["E1_V"] == pytest.approx(
        -0.75 + 0.015 * HARTREE_TO_EV
    )
    assert by_name["1-hydroxyanthraquinone"]["E2_V"] == pytest.approx(
        -1.40 + 0.010 * HARTREE_TO_EV
    )
    assert by_name["2-hydroxyanthraquinone"]["E1_V"] == pytest.approx(
        -0.75 - 0.020 * HARTREE_TO_EV
    )
    assert by_name["2-hydroxyanthraquinone"]["E2_V"] == pytest.approx(
        -1.40 + 0.040 * HARTREE_TO_EV
    )
    assert by_name["2-hydroxyanthraquinone"]["potential_inversion"] is True

    formulas = {
        candidate["name"]: ase.io.read(candidate["path"]).get_chemical_formula()
        for candidate in metadata["candidates"]
    }
    assert formulas == {
        "anthraquinone": "C14H8O2",
        "1-hydroxyanthraquinone": "C14H8O3",
        "2-hydroxyanthraquinone": "C14H8O3",
    }
    assert metadata["input_set_fingerprint"]

    _results, repeated = _run_regression(monkeypatch, tmp_path, "repeated")
    assert repeated["input_set_fingerprint"] == metadata["input_set_fingerprint"]
