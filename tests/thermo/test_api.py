import pytest

import unittest
import numpy as np
import pytest
from pathlib import Path
from unittest.mock import patch, mock_open
from ThermoScreening.calculator.dftbplus import dftb_3ob_parameters
from ThermoScreening.thermo.api import (
    dftbplus_thermo,
    read_coord,
    read_vibrational,
    read_gen,
    read_xyz, 
    run_thermo,
    unit_length,
    unit_mass,
    unit_energy,
    unit_frequency,
)
from ase.build import molecule
from ThermoScreening.exceptions import TSNotImplementedError, TSValueError

class TestApi(unittest.TestCase):
    
    def test_read_coord_not_implemented(self):
        with pytest.raises(TSNotImplementedError) as e:
            read_coord("test.xyz", engine="not_implemented")
        assert str(e.value) == "The engine is not supported."
        
        with pytest.raises(TSNotImplementedError) as e:
            read_coord("test.com", engine="dftb+")
        assert str(e.value) == "The input file is not supported."
        
    def test_read_vibrational_not_implemented(self):
        with pytest.raises(TSNotImplementedError) as e:
            read_vibrational("test.vib", engine="not_implemented")
        assert str(e.value) == "The engine is not supported."
        
    def test_unit_length(self):
        with pytest.raises(TSNotImplementedError) as e:
            unit_length(engine="not_implemented")
        assert str(e.value) == "The engine is not supported."
        
        assert unit_length(engine="dftb+") == "Angstrom"
        
    def test_unit_mass(self):
        with pytest.raises(TSNotImplementedError) as e:
            unit_mass(engine="not_implemented")
        assert str(e.value) == "The engine is not supported."
        
        assert unit_mass(engine="dftb+") == "amu"
        
    def test_unit_energy(self):
        with pytest.raises(TSNotImplementedError) as e:
            unit_energy(engine="not_implemented")
        assert str(e.value) == "The engine is not supported."
        
        assert unit_energy(engine="dftb+") == "Hartree"
        
    def test_unit_frequency(self):
        with pytest.raises(TSNotImplementedError) as e:
            unit_frequency(engine="not_implemented")
        assert str(e.value) == "The engine is not supported."
        
        assert unit_frequency(engine="dftb+") == "cm^-1"
            
   
    @patch(
        "builtins.open",
        new_callable=mock_open,
        read_data="24\n\nO     0.00000003     -0.00000060     -2.14906255      6.45170211\n    O     -0.00000001      0.00000032     -7.56375986      6.45170117\n    C     -3.70943130     -0.00000041     -5.55541507      4.06573254\n    C     -2.50358224      0.00000055     -6.25344004      4.05635404\n    C     -1.28152402      0.00000064     -5.56329596      4.05319169\n    C     -1.28152226     -0.00000048     -4.14952765      4.05319299\n    C     -2.50358591     -0.00000218     -3.45938372      4.05635329\n    C     -3.70942822     -0.00000196     -4.15740715      4.06573042\n    C      0.00000004      0.00000235     -6.33321909      3.57768382\n    C     -0.00000003      0.00000103     -3.37960436      3.57768340\n    C      1.28152286      0.00000564     -4.14952745      4.05319251\n    C      1.28152346      0.00000677     -5.56329581      4.05319217\n    C      2.50358346      0.00001200     -6.25343979      4.05635377\n    H      2.48494730      0.00001319     -7.34030311      0.89492446\n    C      3.70943026      0.00001569     -5.55541540      4.06573185\n    C      3.70942928      0.00001410     -4.15740743      4.06573112\n    C      2.50358473      0.00000924     -3.45938351      4.05635352\n    H     -4.65151464      0.00000015     -6.09787561      0.91510467\n    H     -2.48494736      0.00000157     -7.34030301      0.89492435\n    H     -2.48494723     -0.00000334     -2.37251966      0.89492456\n    H     -4.65151524     -0.00000279     -3.61494721      0.91510621\n    H      4.65151483      0.00001973     -6.09787570      0.91510518\n    H      4.65151507      0.00001672     -3.61494731      0.91510569\n    H      2.48494732      0.00000821     -2.37251974      0.89492448",
    )
    def test_read_coord(self,mock_open):
        
        coord = np.array([[ 0.00000003  ,   -0.00000060  ,   -2.14906255],
                            [-0.00000001 ,    0.00000032  ,   -7.56375986],
                            [-3.70943130 ,   -0.00000041  ,   -5.55541507],
                            [-2.50358224 ,    0.00000055  ,   -6.25344004],
                            [-1.28152402 ,    0.00000064  ,   -5.56329596],
                            [-1.28152226 ,   -0.00000048  ,   -4.14952765],
                            [-2.50358591 ,   -0.00000218  ,   -3.45938372],
                            [-3.70942822 ,   -0.00000196  ,   -4.15740715],
                            [ 0.00000004 ,    0.00000235  ,   -6.33321909],
                            [-0.00000003 ,    0.00000103  ,   -3.37960436],
                            [ 1.28152286 ,    0.00000564  ,   -4.14952745],
                            [ 1.28152346 ,    0.00000677  ,   -5.56329581],
                            [ 2.50358346 ,    0.00001200  ,   -6.25343979],
                            [ 2.48494730 ,    0.00001319  ,   -7.34030311],
                            [ 3.70943026 ,    0.00001569  ,   -5.55541540],
                            [ 3.70942928 ,    0.00001410  ,   -4.15740743],
                            [ 2.50358473 ,    0.00000924  ,   -3.45938351],
                            [-4.65151464 ,    0.00000015  ,   -6.09787561],
                            [-2.48494736 ,    0.00000157  ,   -7.34030301],
                            [-2.48494723 ,   -0.00000334  ,   -2.37251966],
                            [-4.65151524 ,   -0.00000279  ,   -3.61494721],
                            [ 4.65151483 ,    0.00001973  ,   -6.09787570],
                            [ 4.65151507 ,    0.00001672  ,   -3.61494731],
                            [ 2.48494732 ,    0.00000821  ,   -2.37251974]])
        atom_number = 24
        atoms = np.array(
            [
                "O",
                "O",
                "C",
                "C",
                "C",
                "C",
                "C",
                "C",
                "C",
                "C",
                "C",
                "C",
                "C",
                "H",
                "C",
                "C",
                "C",
                "H",
                "H",
                "H",
                "H",
                "H",
                "H",
                "H",
            ]
        )

        cell = None

        data_N, data_atoms, data_xyz, cell, pbc = read_xyz("test.xyz")

        assert data_N == atom_number
        np.testing.assert_array_equal(data_atoms, atoms)
        np.testing.assert_array_almost_equal(data_xyz[:, 0], coord[:, 0], decimal=3)
        assert cell == None
        assert pbc == False

    @patch(
        "builtins.open",
        new_callable=mock_open,
        read_data=(
            "2 10.0 11.0 12.0 90.0 90.0 90.0\n"
            "\n"
            "Cl 0.0 0.0 0.0\n"
            "Na 1.0 2.0 3.0\n"
        ),
    )
    def test_read_xyz_keeps_multichar_symbols(self, mock_open):
        data_N, data_atoms, data_xyz, cell, pbc = read_xyz("test.xyz")

        assert data_N == 2
        np.testing.assert_array_equal(data_atoms, np.array(["Cl", "Na"], dtype=object))
        np.testing.assert_allclose(data_xyz, np.array([[0.0, 0.0, 0.0], [1.0, 2.0, 3.0]]))
        np.testing.assert_allclose(cell, np.array([10.0, 11.0, 12.0, 90.0, 90.0, 90.0]))
        assert pbc is True

    @patch(
        "builtins.open",
        new_callable=mock_open,
        read_data=(
            "1 bad 11.0 12.0 90.0 90.0 90.0\n"
            "\n"
            "Cl 0.0 0.0 0.0\n"
        ),
    )
    def test_read_xyz_ignores_invalid_cell_header(self, mock_open):
        data_N, data_atoms, data_xyz, cell, pbc = read_xyz("test.xyz")

        assert data_N == 1
        np.testing.assert_array_equal(data_atoms, np.array(["Cl"], dtype=object))
        np.testing.assert_allclose(data_xyz, np.array([[0.0, 0.0, 0.0]]))
        assert cell is None
        assert pbc is False

    @patch(
        "builtins.open",
        new_callable=mock_open,
        read_data="1.0  2.0\n2.0  3.0\n3.0  4.0\n4.0  5.0\n5.0  6.0\n6.0  7.0\n",
    )
    def test_read_vib_file(self, mock_open):
        vibrational_frequencies = np.array([2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
        data_vib = read_vibrational("test.vib", "dftb+")
        np.testing.assert_array_almost_equal(
            data_vib, vibrational_frequencies, decimal=8
        )

    def test_read_gen(self):
        gen_file = Path(__file__).resolve().parents[1] / "data/thermo/geo_opt.gen"

        data_N, data_atoms, data_xyz, cell, pbc = read_gen(str(gen_file))

        assert data_N == 24
        assert pbc is True
        np.testing.assert_array_equal(data_atoms[:4], np.array(["O", "O", "C", "C"]))
        np.testing.assert_allclose(
            data_xyz[0],
            np.array([3.0e-08, -6.0e-07, -2.14906255]),
            atol=1.0e-12,
        )
        np.testing.assert_allclose(
            cell,
            np.array(
                [
                    [10000.0, 0.0, 0.0],
                    [0.0, 10000.0, 0.0],
                    [0.0, 0.0, 10000.0],
                ]
            ),
        )

        coord_data = read_coord(str(gen_file), engine="dftb+")
        assert coord_data[0] == data_N
        np.testing.assert_array_equal(coord_data[1], data_atoms)
        np.testing.assert_allclose(coord_data[2], data_xyz)
        np.testing.assert_allclose(coord_data[3], cell)
        assert coord_data[4] is True

    @patch(
        "builtins.open",
        new_callable=mock_open,
        read_data=(
            "2 C\n"
            "Cl Na\n"
            "1 1 0.0 0.0 0.0\n"
            "2 2 1.0 2.0 3.0\n"
        ),
    )
    def test_read_gen_keeps_multichar_symbols(self, mock_open):
        data_N, data_atoms, data_xyz, cell, pbc = read_gen("test.gen")

        assert data_N == 2
        np.testing.assert_array_equal(data_atoms, np.array(["Cl", "Na"], dtype=object))
        np.testing.assert_allclose(data_xyz, np.array([[0.0, 0.0, 0.0], [1.0, 2.0, 3.0]]))
        assert cell is None
        assert pbc is False

    def test_run_thermo_rejects_wrong_frequency_count(self):
        coord_file = Path(__file__).resolve().parents[1] / "data/thermo/geo_opt.xyz"

        with pytest.raises(TSValueError, match="number of vibrational frequencies"):
            run_thermo(
                np.array([1.0]),
                coord_file=str(coord_file),
                engine="dftb+",
            )

if __name__ == "__main__":
    unittest.main()
