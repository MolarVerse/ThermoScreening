import unittest
import numpy as np
import pytest
from unittest.mock import patch, mock_open
from ThermoScreening.thermo.api import (
    read_coord,
    read_vibrational,
    read_gen,
    read_xyz,
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
            
   
    @patch("builtins.open",new_callable=mock_open,read_data="24\n\nO     0.00000003     -0.00000060     -2.14906255      6.45170211\n    O     -0.00000001      0.00000032     -7.56375986      6.45170117\n    C     -3.70943130     -0.00000041     -5.55541507      4.06573254\n    C     -2.50358224      0.00000055     -6.25344004      4.05635404\n    C     -1.28152402      0.00000064     -5.56329596      4.05319169\n    C     -1.28152226     -0.00000048     -4.14952765      4.05319299\n    C     -2.50358591     -0.00000218     -3.45938372      4.05635329\n    C     -3.70942822     -0.00000196     -4.15740715      4.06573042\n    C      0.00000004      0.00000235     -6.33321909      3.57768382\n    C     -0.00000003      0.00000103     -3.37960436      3.57768340\n    C      1.28152286      0.00000564     -4.14952745      4.05319251\n    C      1.28152346      0.00000677     -5.56329581      4.05319217\n    C      2.50358346      0.00001200     -6.25343979      4.05635377\n    H      2.48494730      0.00001319     -7.34030311      0.89492446\n    C      3.70943026      0.00001569     -5.55541540      4.06573185\n    C      3.70942928      0.00001410     -4.15740743      4.06573112\n    C      2.50358473      0.00000924     -3.45938351      4.05635352\n    H     -4.65151464      0.00000015     -6.09787561      0.91510467\n    H     -2.48494736      0.00000157     -7.34030301      0.89492435\n    H     -2.48494723     -0.00000334     -2.37251966      0.89492456\n    H     -4.65151524     -0.00000279     -3.61494721      0.91510621\n    H      4.65151483      0.00001973     -6.09787570      0.91510518\n    H      4.65151507      0.00001672     -3.61494731      0.91510569\n    H      2.48494732      0.00000821     -2.37251974      0.89492448")
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
        atoms = np.array(['O', 'O', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'H', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H'])
    
        cell = None


        data_N,data_atoms,data_xyz,cell,pbc = read_xyz("test.xyz")
       
      
        assert data_N == atom_number
        np.testing.assert_array_equal(data_atoms, atoms)
        np.testing.assert_array_almost_equal(data_xyz[:,0],coord[:,0],decimal=3)
        assert cell == None
        assert pbc == False


    # @patch("builtins.open",new_callable=mock_open,read_data="24 S\nO C H\n 1   1   3.0000000000E-08  -6.0000000000E-07  -2.1490625500E+00\n 2   1  -1.0000000000E-08   3.2000000000E-07  -7.5637598600E+00\n 3   2  -3.7094313000E+00  -4.1000000000E-07  -5.5554150700E+00\n 4   2  -2.5035822400E+00   5.5000000000E-07  -6.2534400400E+00\n 5   2  -1.2815240200E+00   6.4000000000E-07  -5.5632959600E+00\n 6   2  -1.2815222600E+00  -4.8000000000E-07  -4.1495276500E+00\n 7   2  -2.5035859100E+00  -2.1800000000E-06  -3.4593837200E+00\n 8   2  -3.7094282200E+00  -1.9600000000E-06  -4.1574071500E+00\n 9   2   4.0000000000E-08   2.3500000000E-06  -6.3332190900E+00\n10   2  -3.0000000000E-08   1.0300000000E-06  -3.3796043600E+00\n11   2   1.2815228600E+00   5.6400000000E-06  -4.1495274500E+00\n12   2   1.2815234600E+00   6.7700000000E-06  -5.5632958100E+00\n13   2   2.5035834600E+00   1.2000000000E-05  -6.2534397900E+00\n14   3   2.4849473000E+00   1.3190000000E-05  -7.3403031100E+00\n15   2   3.7094302600E+00   1.5690000000E-05  -5.5554154000E+00\n16   2   3.7094292800E+00   1.4100000000E-05  -4.1574074300E+00\n17   2   2.5035847300E+00   9.2400000000E-06  -3.4593835100E+00\n18   3  -4.6515146400E+00   1.5000000000E-07  -6.0978756100E+00\n19   3  -2.4849473600E+00   1.5700000000E-06  -7.3403030100E+00\n20   3  -2.4849472300E+00  -3.3400000000E-06  -2.3725196600E+00\n21   3  -4.6515152400E+00  -2.7900000000E-06  -3.6149472100E+00\n22   3   4.6515148300E+00   1.9730000000E-05  -6.0978757000E+00\n23   3   4.6515150700E+00   1.6720000000E-05  -3.6149473100E+00\n24   3   2.4849473200E+00   8.2100000000E-06  -2.3725197400E+00\n0 0 0\n10000 0 0\n0 10000 0\n0 0 10000\n")

    # def test_read_coord_pbc(self,mock_open):
    #     coord = np.array([[ 0.00000003  ,   -0.00000060  ,   -2.14906255],
    #                         [-0.00000001 ,    0.00000032  ,   -7.56375986],
    #                         [-3.70943130 ,   -0.00000041  ,   -5.55541507],
    #                         [-2.50358224 ,    0.00000055  ,   -6.25344004],
    #                         [-1.28152402 ,    0.00000064  ,   -5.56329596],
    #                         [-1.28152226 ,   -0.00000048  ,   -4.14952765],
    #                         [-2.50358591 ,   -0.00000218  ,   -3.45938372],
    #                         [-3.70942822 ,   -0.00000196  ,   -4.15740715],
    #                         [ 0.00000004 ,    0.00000235  ,   -6.33321909],
    #                         [-0.00000003 ,    0.00000103  ,   -3.37960436],
    #                         [ 1.28152286 ,    0.00000564  ,   -4.14952745],
    #                         [ 1.28152346 ,    0.00000677  ,   -5.56329581],
    #                         [ 2.50358346 ,    0.00001200  ,   -6.25343979],
    #                         [ 2.48494730 ,    0.00001319  ,   -7.34030311],
    #                         [ 3.70943026 ,    0.00001569  ,   -5.55541540],
    #                         [ 3.70942928 ,    0.00001410  ,   -4.15740743],
    #                         [ 2.50358473 ,    0.00000924  ,   -3.45938351],
    #                         [-4.65151464 ,    0.00000015  ,   -6.09787561],
    #                         [-2.48494736 ,    0.00000157  ,   -7.34030301],
    #                         [-2.48494723 ,   -0.00000334  ,   -2.37251966],
    #                         [-4.65151524 ,   -0.00000279  ,   -3.61494721],
    #                         [ 4.65151483 ,    0.00001973  ,   -6.09787570],
    #                         [ 4.65151507 ,    0.00001672  ,   -3.61494731],
    #                         [ 2.48494732 ,    0.00000821  ,   -2.37251974]])
    #     atom_number = 24
    #     atoms = np.array(['O', 'O', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'H', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H'])
    #     cell = np.array([[10000,0,0],[0,10000,0],[0,0,10000]])
    #     data_N,data_atoms,data_xyz,cell_vectors,pbc= read_gen("test.gen")

    #     assert data_N == atom_number
    #     np.testing.assert_array_equal(data_atoms, atoms)
    #     np.testing.assert_array_almost_equal(data_xyz,coord,decimal=8)
    #     np.testing.assert_array_almost_equal(cell,cell_vectors,decimal=8)
    #     assert pbc == True
        

    @patch('builtins.open',new_callable=mock_open,read_data="1.0  2.0\n2.0  3.0\n3.0  4.0\n4.0  5.0\n5.0  6.0\n6.0  7.0\n")
    def test_read_vib_file(self,mock_open):
        vibrational_frequencies = np.array([2.0,3.0,4.0,5.0,6.0,7.0])
        data_vib = read_vibrational("test.vib", "dftb+")
        np.testing.assert_array_almost_equal(data_vib,vibrational_frequencies,decimal=8)
         
       


if __name__ == '__main__':
    unittest.main()
    
