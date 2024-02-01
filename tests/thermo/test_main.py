
import pytest
import unittest as unittest
from unittest.mock import patch
from filecmp import cmp as filecmp

import sys
import argparse

import ThermoScreening.cli.thermo as thermo

sys.path.append("/home/stk/dev/ThermoScreening/tests/data/thermo/")

path = "/home/stk/dev/ThermoScreening/tests/data/thermo/"
#


class TestMain(unittest.TestCase):
    def test_parse_args(self):
        test_args = ['thermo', '-i', 'input_file', '-o', 'output_file', '-p', 'plot_file', '-v', '-t']
        with patch.object(sys, 'argv', test_args):
            args = thermo.parse_args()
            self.assertEqual(args.input, 'input_file')
            self.assertEqual(args.output, 'output_file')
            self.assertEqual(args.plot, 'plot_file')
            self.assertTrue(args.verbose)
            self.assertTrue(args.test)


if __name__ == '__main__':
    unittest.main()