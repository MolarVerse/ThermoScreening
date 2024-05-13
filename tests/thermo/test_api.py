import pytest

import unittest
import numpy as np
import pytest
from unittest.mock import patch, mock_open
from ThermoScreening.calculator.dftbplus import dftb_3ob_parameters
from ThermoScreening.thermo.api import (
     read_coord,
     read_vibrational,
     read_gen,
     read_xyz, 
     dftbplus_thermo
)
from ase.build import (
    molecule,
    unit_length,
    unit_mass,
    unit_energy,
    unit_frequency,
)
from ThermoScreening.exceptions import TSNotImplementedError

class TestApi(unittest.TestCase):
    
    def test_read_coord_not_implemented(self):
        with pytest.raises(TSNotImplementedError) as e:
            read_coord("test.xyz", engine="not_implemented")
        assert str(e.value) == "The engine is not supported."
        
        with pytest.raises(TSNotImplementedError) as e:
            read_coord("test.gen", engine="dftb+")
        assert str(e.value) == "The gen file is not tested yet."
        
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
        read_data="1.0  2.0\n2.0  3.0\n3.0  4.0\n4.0  5.0\n5.0  6.0\n6.0  7.0\n",
    )
    def test_read_vib_file(self, mock_open):
        vibrational_frequencies = np.array([2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
        data_vib = read_vibrational("test.vib", "dftb+")
        np.testing.assert_array_almost_equal(
            data_vib, vibrational_frequencies, decimal=8
        )

if __name__ == "__main__":
    unittest.main()
