"""Mandatory DFTB+ release-gate workflow."""

import json
import math
import os
from pathlib import Path
import subprocess
import sys

import pytest


@pytest.mark.integration
def test_dftb_screen_end_to_end(tmp_path):
    if os.getenv("THERMOSCREENING_REQUIRE_DFTB") != "1":
        pytest.skip("set THERMOSCREENING_REQUIRE_DFTB=1 for the release gate")

    structure = (
        Path(__file__).resolve().parents[1] / "data" / "regression" / "water.xyz"
    )
    manifest = tmp_path / "molecules.csv"
    manifest.write_text(
        f"name,path,charge\nwater,{structure},0\n",
        encoding="utf-8",
    )
    out = tmp_path / "results"
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "ThermoScreening",
            "screen",
            str(manifest),
            "--out",
            str(out),
            "--directory",
            str(tmp_path / "calculations"),
            "--engine",
            "dftb+",
        ],
        capture_output=True,
        text=True,
        timeout=300,
        check=False,
    )

    if completed.returncode:
        calculation_dir = tmp_path / "calculations" / "water"
        diagnostics = [completed.stdout, completed.stderr]
        for filename in ("dftb_in.hsd", "second_derivative.out", "detailed.out"):
            path = calculation_dir / filename
            if path.is_file():
                diagnostics.append(
                    f"\n--- {filename} ---\n{path.read_text(encoding='utf-8')}"
                )
        pytest.fail("".join(diagnostics))
    records = json.loads(out.with_suffix(".json").read_text(encoding="utf-8"))
    assert len(records) == 1
    record = records[0]
    assert record["name"] == "water"
    assert record["formula"] == "H2O"
    assert record["status"] == "ok"
    for field in (
        "Eelec_hartree",
        "G_total_hartree",
        "S_cal_per_mol_K",
        "Cv_cal_per_mol_K",
    ):
        assert math.isfinite(float(record[field]))

    metadata = json.loads(
        out.with_name(f"{out.name}-run.json").read_text(encoding="utf-8")
    )
    assert metadata["workflow"] == "thermochemistry_screen"
    assert metadata["settings"]["engine"] == "dftb+"
    assert metadata["input_set_fingerprint"]
    assert metadata["jobs"][0]["structure_sha256"]
