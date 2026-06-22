import unittest
import numpy as np
from unittest.mock import patch, mock_open
from  ThermoScreening.thermo.cell import Cell, cell_parameters_calc
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

    def test_cell_parameters_calc(self):
        cell_vector = np.array([[2, 0, 0], [0, 3, 0], [0, 0, 4]])

        np.testing.assert_allclose(
            cell_parameters_calc(cell_vector),
            np.array([2, 3, 4, 90, 90, 90]),
        )

    def test_cell_scalar_properties_and_repr(self):
        cell = Cell(a=1, b=2, c=3, alpha=80, beta=90, gamma=100)

        assert cell.a == 1
        assert cell.b == 2
        assert cell.c == 3
        assert cell.alpha == 80
        assert cell.beta == 90
        assert cell.gamma == 100
        assert repr(cell) == "Cell(a=1, b=2, c=3, alpha=80, beta=90, gamma=100)"

if __name__ == '__main__':
    unittest.main()
