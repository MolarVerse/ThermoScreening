import os
import shutil

import numpy as np
import pytest
from ase import Atoms

from ThermoScreening.calculator import xtb_cli


# --- executable resolution --------------------------------------------------- #

def test_resolve_xtb_prefers_explicit_then_env(monkeypatch):
    monkeypatch.setenv("XTB_COMMAND", "/opt/xtb")
    assert xtb_cli.resolve_xtb("/usr/bin/xtb") == "/usr/bin/xtb"  # explicit wins
    assert xtb_cli.resolve_xtb() == "/opt/xtb"                    # then env


def test_resolve_xtb_missing_raises(monkeypatch):
    monkeypatch.delenv("XTB_COMMAND", raising=False)
    monkeypatch.setattr(xtb_cli.shutil, "which", lambda command: None)
    with pytest.raises(FileNotFoundError, match="xtb"):
        xtb_cli.resolve_xtb()


# --- output parsing ---------------------------------------------------------- #

def test_parse_vibspectrum_takes_wavenumber_and_sorts(tmp_path):
    vib = tmp_path / "vibspectrum"
    vib.write_text(
        "$vibrational spectrum\n"
        "#  mode  symmetry  wave number  IR intensity  selection rules\n"
        "     1                      -0.00         0.00000          - \n"
        "     6                       0.00         0.00000          - \n"
        "     7        a           1541.16       133.03481         YES\n"
        "     8        a           3635.99         6.71206         YES\n"
        "$end\n",
        encoding="utf-8",
    )
    out = xtb_cli._parse_vibspectrum(str(vib))

    assert len(out) == 4
    assert np.all(np.diff(out) >= 0)          # sorted ascending
    assert list(out[-2:]) == [1541.16, 3635.99]  # real modes on top


def test_parse_optimised_energy(tmp_path):
    xyz = tmp_path / "xtbopt.xyz"
    xyz.write_text(
        "2\n energy: -4.438237080538 gnorm: 0.0004 xtb: 6.7.1\n"
        "O 0.0 0.0 0.0\nH 0.0 0.0 0.97\n",
        encoding="utf-8",
    )
    assert xtb_cli._parse_optimised_energy(str(xyz)) == pytest.approx(-4.438237080538)


# --- run_xtb (mocked subprocess, no binary needed) --------------------------- #

def test_run_xtb_builds_command_and_parses(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("XTB_COMMAND", "/fake/xtb")
    captured = {}

    def fake_run(argv, **kwargs):
        captured["argv"] = argv
        (tmp_path / "xtbopt.xyz").write_text(
            "2\n energy: -4.4 gnorm: 0 xtb: x\nO 0 0 0\nH 0 0 0.97\n", encoding="utf-8"
        )
        (tmp_path / "vibspectrum").write_text(
            "$vibrational spectrum\n     1   -0.00  0.0  -\n"
            "     6   a  3500.0  1.0  YES\n$end\n", encoding="utf-8"
        )
        return type("R", (), {"returncode": 0, "stdout": "", "stderr": ""})()

    monkeypatch.setattr(xtb_cli.subprocess, "run", fake_run)

    atoms, energy, freqs = xtb_cli.run_xtb(
        Atoms("OH", positions=[[0, 0, 0], [0, 0, 0.97]]),
        charge=-1, unpaired=1, method="GFN2-xTB", solvent="water",
    )

    assert energy == pytest.approx(-4.4)
    assert freqs[-1] == pytest.approx(3500.0)
    argv = captured["argv"]
    assert "--ohess" in argv
    assert argv[argv.index("--gfn") + 1] == "2"
    assert argv[argv.index("--uhf") + 1] == "1"
    assert argv[argv.index("--chrg") + 1] == "-1"
    assert argv[argv.index("--alpb") + 1] == "water"


def test_parse_fukui_extracts_table():
    output = (
        "some preceding SCF chatter\n"
        "\n"
        "Fukui functions:\n"
        "     #        f(+)     f(-)     f(0)\n"
        "     1C       0.024    0.024    0.024\n"
        "     8O       0.131    0.129    0.130\n"
        "\n"
        "trailing text after the table\n"
    )
    rows = xtb_cli._parse_fukui(output)

    assert rows == [("C", 0.024, 0.024, 0.024), ("O", 0.131, 0.129, 0.130)]


def test_parse_fukui_raises_without_table():
    with pytest.raises(ValueError, match="No 'Fukui functions:' table"):
        xtb_cli._parse_fukui("no such table here")


def test_run_xtb_fukui_builds_command_and_parses(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("XTB_COMMAND", "/fake/xtb")
    captured = {}

    def fake_run(argv, **kwargs):
        captured["argv"] = argv
        stdout = "Fukui functions:\n     #        f(+)     f(-)     f(0)\n     1O       0.5      0.5      0.5\n"
        return type("R", (), {"returncode": 0, "stdout": stdout, "stderr": ""})()

    monkeypatch.setattr(xtb_cli.subprocess, "run", fake_run)

    rows = xtb_cli.run_xtb_fukui(
        Atoms("OH", positions=[[0, 0, 0], [0, 0, 0.97]]),
        charge=-1, unpaired=1, method="GFN2-xTB", solvent="water",
    )

    assert rows == [("O", 0.5, 0.5, 0.5)]
    argv = captured["argv"]
    assert "--vfukui" in argv
    assert "--ohess" not in argv
    assert argv[argv.index("--gfn") + 1] == "2"
    assert argv[argv.index("--uhf") + 1] == "1"
    assert argv[argv.index("--chrg") + 1] == "-1"
    assert argv[argv.index("--alpb") + 1] == "water"


def test_run_xtb_fukui_rejects_unknown_method(monkeypatch):
    with pytest.raises(ValueError, match="Unknown xTB method"):
        xtb_cli.run_xtb_fukui(Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.7]]), method="B3LYP")


def test_run_xtb_fukui_raises_on_failure(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("XTB_COMMAND", "/fake/xtb")
    monkeypatch.setattr(
        xtb_cli.subprocess, "run",
        lambda argv, **kw: type("R", (), {"returncode": 1, "stdout": "boom", "stderr": "e"})(),
    )
    with pytest.raises(RuntimeError, match="xtb failed"):
        xtb_cli.run_xtb_fukui(Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.7]]))


def test_run_xtb_rejects_unknown_method(monkeypatch):
    # method check happens before the executable is even resolved
    with pytest.raises(ValueError, match="Unknown xTB method"):
        xtb_cli.run_xtb(Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.7]]), method="B3LYP")


def test_run_xtb_raises_on_failure(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("XTB_COMMAND", "/fake/xtb")
    monkeypatch.setattr(
        xtb_cli.subprocess, "run",
        lambda argv, **kw: type("R", (), {"returncode": 1, "stdout": "boom", "stderr": "e"})(),
    )
    with pytest.raises(RuntimeError, match="xtb failed"):
        xtb_cli.run_xtb(Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.7]]))


# --- real integration (needs the xtb binary) -------------------------------- #

xtb_available = shutil.which("xtb") is not None or "XTB_COMMAND" in os.environ


@pytest.mark.skipif(not xtb_available, reason="the native xtb binary is not available.")
def test_xtb_cli_radical_solvation_end_to_end(tmp_path):
    from ThermoScreening.thermo.api import xtb_cli_thermo

    def oh():
        return Atoms("OH", positions=[[0, 0, 0], [0, 0, 0.97]])

    gas = xtb_cli_thermo(oh(), directory=str(tmp_path / "gas"))
    solvated = xtb_cli_thermo(oh(), directory=str(tmp_path / "sol"), solvent="water")

    # open-shell radical (doublet: S ~ 42-43 cal/mol/K incl. R ln2)
    assert gas.total_entropy("cal/(mol*K)") == pytest.approx(42.5, abs=2.0)
    # solvation stabilises the radical
    assert solvated.total_EeGtot() < gas.total_EeGtot()


@pytest.mark.skipif(not xtb_available, reason="the native xtb binary is not available.")
def test_xtb_fukui_indices_end_to_end(tmp_path):
    import ase.io
    from ThermoScreening.thermo.api import xtb_cli_thermo, xtb_fukui_indices

    water = Atoms("OH2", positions=[[0, 0, 0.119], [0, 0.763, -0.477], [0, -0.763, -0.477]])
    opt_dir = tmp_path / "opt"
    xtb_cli_thermo(water, directory=str(opt_dir))
    optimized = ase.io.read(str(opt_dir / "xtbopt.xyz"))

    rows = xtb_fukui_indices(optimized, directory=str(tmp_path / "fukui"))

    assert len(rows) == 3
    symbols = [symbol for symbol, *_ in rows]
    assert symbols == ["O", "H", "H"]
    # Fukui functions are normalised: each column sums to ~1 over all atoms
    assert sum(f_plus for _, f_plus, _, _ in rows) == pytest.approx(1.0, abs=0.05)
    assert sum(f_minus for _, _, f_minus, _ in rows) == pytest.approx(1.0, abs=0.05)
    # the oxygen lone pairs make O the dominant site for electron loss (f_minus,
    # i.e. oxidation); the two symmetric H atoms dominate electron gain (f_plus)
    oxygen_row = rows[0]
    hydrogen_rows = rows[1:]
    assert oxygen_row[2] > max(row[2] for row in hydrogen_rows)
    assert oxygen_row[1] < min(row[1] for row in hydrogen_rows)
