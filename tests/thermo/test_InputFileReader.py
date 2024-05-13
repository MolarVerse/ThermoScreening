import unittest
import numpy as np
from unittest.mock import patch, mock_open
from ThermoScreening.thermo.inputFileReader import InputFileReader
from ThermoScreening.exceptions import TSValueError

class TestFileReader(unittest.TestCase):
    @patch("builtins.open",new_callable=mock_open,read_data="coord_file = geo_opt.xyz\ntemperature = 298.15\npressure = 101325\nengine = dftb+\nvibrational_file = frequency.txt\nenergy = -33.6052447996")
    def test_read(self,mock_open):
        """
        Tests if the input file is read correctly.
        """
        input_file_reader = InputFileReader("test.in")
        assert input_file_reader._input_file == "test.in"
        assert input_file_reader._raw_input_file == ["coord_file = geo_opt.xyz\n","temperature = 298.15\n","pressure = 101325\n","engine = dftb+\n","vibrational_file = frequency.txt\n","energy = -33.6052447996"]
        assert input_file_reader._dictionary == {"coord_file":"geo_opt.xyz","temperature":"298.15","pressure":"101325","engine":"dftb+","vibrational_file":"frequency.txt","energy":"-33.6052447996"}
        
    def test_not_input_file(self):
        """
        Tests if the input file is not given.
        """
        with self.assertRaises(TSValueError) as e:
            InputFileReader()
        self.assertEqual(str(e.exception),"The input file has to be given to initialize the InputFileReader.")
        
    @patch("builtins.open",new_callable=mock_open,read_data="temperature = 298.15\npressure = 101325\nengine = dftb+\nvibrational_file = frequency.txt\nenergy = -33.6052447996")
    def test_not_all_required_keys(self, mock_open):
        """
        Tests if not all required keys are given.
        """
        with self.assertRaises(TSValueError) as e:
            InputFileReader("test.in")
        self.assertEqual(str(e.exception),"The key coord_file is not set in the input file.")
        
    @patch("builtins.open",new_callable=mock_open,read_data="_unknown_key_ = 21\ncoord_file = geo_opt.xyz\ntemperature = 298.15\npressure = 101325\nengine = dftb+\nvibrational_file = frequency.txt\nenergy = -33.6052447996")
    def test_not_unknown_keys(self, mock_open):
        """
        Tests if not all required keys are given.
        """
        with self.assertRaises(TSValueError) as e:
            InputFileReader("test.in")
        self.assertEqual(str(e.exception),"The key _unknown_key_ is not known.")

if __name__ == '__main__':
    unittest.main()
    
