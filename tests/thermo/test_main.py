
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

if __name__ == '__main__':
    unittest.main()
