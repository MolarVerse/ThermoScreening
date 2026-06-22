
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

if __name__ == '__main__':
    unittest.main()
