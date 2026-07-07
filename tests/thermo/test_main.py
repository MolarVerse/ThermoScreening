
import pytest
import runpy
import unittest as unittest
from io import StringIO
from unittest.mock import patch
from filecmp import cmp as filecmp

import sys
import argparse

import ThermoScreening.cli.thermo as thermo
from ThermoScreening.utils.header import print_header

sys.path.append("/home/stk/dev/ThermoScreening/tests/data/thermo/")

path = "/home/stk/dev/ThermoScreening/tests/data/thermo/"
#


class TestMain(unittest.TestCase):
    def test_parse_args(self):
        test_args = ['thermo', 'input_file', '-v']
        with patch.object(sys, 'argv', test_args):
            args = thermo.parse_args()
            self.assertEqual(args.input_file, 'input_file')
            self.assertTrue(args.verbose)

    def test_module_help_entrypoint(self):
        test_args = ['ThermoScreening', '--help']
        with patch.object(sys, 'argv', test_args):
            with self.assertRaises(SystemExit) as e:
                runpy.run_module("ThermoScreening.__main__", run_name="__main__")
        self.assertEqual(e.exception.code, 0)

    def test_print_header_to_file(self):
        output = StringIO()
        print_header(file=output)

        self.assertIn("ThermoScreening - v", output.getvalue())


def test_main_executes_input_file_with_timing(monkeypatch, capsys):
    calls = []
    times = iter([10.0, 12.5])

    monkeypatch.setattr(
        thermo,
        "parse_args",
        lambda: argparse.Namespace(input_file="thermo.in", verbose=True),
    )
    monkeypatch.setattr(thermo, "execute", lambda input_file: calls.append(input_file))
    monkeypatch.setattr(thermo.time, "time", lambda: next(times))

    assert thermo.main() is None

    assert calls == ["thermo.in"]
    output = capsys.readouterr().out
    assert "Input file:  thermo.in" in output
    assert "Verbose:  True" in output
    assert "Time elapsed:  2.5  s" in output


def test_main_executes_input_file_without_verbose(monkeypatch, capsys):
    calls = []

    monkeypatch.setattr(
        thermo,
        "parse_args",
        lambda: argparse.Namespace(input_file="thermo.in", verbose=False),
    )
    monkeypatch.setattr(thermo, "execute", lambda input_file: calls.append(input_file))

    assert thermo.main() is None

    assert calls == ["thermo.in"]
    assert capsys.readouterr().out == ""


def test_parse_args_setup_dftb_command():
    args = thermo.parse_args(
        [
            "setup-dftb",
            "--install-root",
            "/tmp/slakos",
            "--url",
            "file:///tmp/3ob-3-1.tar.xz",
            "--force",
            "--parameter-set",
            "mio",
        ]
    )

    assert args.command == "setup-dftb"
    assert args.install_root == "/tmp/slakos"
    assert args.url == "file:///tmp/3ob-3-1.tar.xz"
    assert args.force is True
    assert args.parameter_set == "mio"


def test_parse_args_doctor_command():
    args = thermo.parse_args(["doctor"])

    assert args.command == "doctor"


def test_main_runs_setup_dftb(monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(
        thermo,
        "parse_args",
        lambda: argparse.Namespace(
            command="setup-dftb",
            install_root=str(tmp_path),
            url="file:///tmp/3ob-3-1.tar.xz",
            force=True,
            parameter_set="3ob",
            solvent=None,
        ),
    )
    monkeypatch.setattr(
        thermo,
        "install_slakos",
        lambda install_root, url, force, parameter_set: tmp_path / "3ob-3-1",
    )

    assert thermo.main() == 0

    output = capsys.readouterr().out
    assert "Slater-Koster files:" in output
    assert "export DFTB_PREFIX=" in output


def test_main_setup_dftb_downloads_solvent(monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(
        thermo,
        "parse_args",
        lambda: argparse.Namespace(
            command="setup-dftb",
            install_root=str(tmp_path),
            url=None,
            force=False,
            parameter_set="3ob",
            solvent="water",
        ),
    )
    monkeypatch.setattr(
        thermo,
        "install_gbsa_param",
        lambda solvent, install_root, url, force: tmp_path / "param_gbsa_h2o.txt",
    )

    assert thermo.main() == 0
    assert "GBSA solvation parameters:" in capsys.readouterr().out


def test_main_runs_doctor(monkeypatch, capsys):
    diagnostic = type("DiagnosticStub", (), {"ok": False, "optional": False})()

    monkeypatch.setattr(
        thermo,
        "parse_args",
        lambda: argparse.Namespace(command="doctor"),
    )
    monkeypatch.setattr(thermo, "check_dftb_setup", lambda: [diagnostic])
    monkeypatch.setattr(thermo, "format_diagnostics", lambda diagnostics: "not ready")

    assert thermo.main() == 1
    assert capsys.readouterr().out == "not ready\n"


def test_run_screen_forwards_dispersion(monkeypatch):
    captured = {}

    def fake_screen(source, **kwargs):
        captured.update(kwargs)
        return []

    monkeypatch.setattr(thermo, "screen", fake_screen)
    args = argparse.Namespace(
        source="mols", out="results", charge=0.0, temperature=298.15,
        pressure=101325, directory="screening", parameter_set="3ob",
        solvent=None, dispersion="d3-bj", quasi_rrho=False, engine="dftb+",
        method="GFN2-xTB", resume=False, jobs=1,
    )
    assert thermo.run_screen(args) == 0
    assert captured["dispersion"] == "d3-bj"


def test_parse_args_screen_dispersion_choice():
    assert thermo.parse_args(["screen", "mols"]).dispersion is None
    assert thermo.parse_args(["screen", "mols", "--dispersion", "d3-bj"]).dispersion == "d3-bj"


def test_parse_args_screen_jobs():
    assert thermo.parse_args(["screen", "mols"]).jobs == 1
    assert thermo.parse_args(["screen", "mols", "-j", "4"]).jobs == 4
    assert thermo.parse_args(["screen", "mols", "--jobs", "8"]).jobs == 8


def test_main_runs_conformers(monkeypatch):
    monkeypatch.setattr(
        thermo, "parse_args", lambda: argparse.Namespace(command="conformers")
    )
    called = {}

    def fake_run_conformers(args):
        called["ran"] = True
        return 0

    monkeypatch.setattr(thermo, "run_conformers", fake_run_conformers)

    assert thermo.main() == 0
    assert called["ran"] is True


def test_main_doctor_ignores_missing_optional_backend(monkeypatch, capsys):
    required_ok = type("D", (), {"ok": True, "optional": False})()
    optional_missing = type("D", (), {"ok": False, "optional": True})()

    monkeypatch.setattr(
        thermo, "parse_args", lambda: argparse.Namespace(command="doctor")
    )
    monkeypatch.setattr(
        thermo, "check_dftb_setup", lambda: [required_ok, optional_missing]
    )
    monkeypatch.setattr(thermo, "format_diagnostics", lambda diagnostics: "ok")

    # a missing optional backend (xtb/tblite) must not fail the doctor
    assert thermo.main() == 0


if __name__ == '__main__':
    unittest.main()
