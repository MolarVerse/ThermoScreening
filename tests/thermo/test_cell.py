import unittest
import numpy as np
from unittest.mock import patch, mock_open
from  ThermoScreening.thermo.cell import Cell 
import pytest

class TestCell(unittest.TestCase):
   
    def test_cell(self):
        cell = Cell(a=1, b=2, c=3, alpha=90, beta=90, gamma=90)
        expected_vectors = np.array([[1, 0, 0], [0, 2, 0], [0, 0, 3]])
        np.testing.assert_array_equal(cell.cell_vectors, expected_vectors)

        expected_parameters = np.array([1, 2, 3, 90, 90, 90])
        np.testing.assert_array_equal(cell.cell_parameters, expected_parameters)


        expected_volume = 1 * 2 * 3  # a * b * c
        self.assertEqual(cell.volume, expected_volume)


        expected_inverse = np.linalg.inv(cell.cell_vectors)
        np.testing.assert_array_equal(cell.cell_inverse, expected_inverse)

if __name__ == '__main__':
    unittest.main()