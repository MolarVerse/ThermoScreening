import unittest
import numpy as np
from ThermoScreening.thermo.system import System, dim, dof, linearity
from ThermoScreening.thermo.atoms import Atom
from ThermoScreening.thermo.cell import Cell
from ThermoScreening.exceptions import TSValueError
import pytest

def test_dim():
    atoms = []
    with pytest.raises(TSValueError) as e:
        dim(atoms)
    assert str(e.value) == "The number of atoms must be greater than 0."
    
def test_dof():
    atoms = []
    with pytest.raises(TSValueError) as e:
        dof(atoms)
    assert str(e.value) == "The number of atoms must be greater than 0."
    
def test_linearity():
    atoms = [Atom(symbol='H', position=np.array([0, 0, 0]))]
    
    with pytest.raises(TSValueError) as e:
        linearity(atoms)
    assert str(e.value) == "Number of atoms must be greater than 1. The system is monoatomic."
             

class TestSystem(unittest.TestCase):
    atoms = [
                Atom(symbol='H', position=np.array([0, 0, 0])),
                Atom(symbol='H', position=np.array([1, 0, 0])),
                Atom(symbol='H', position=np.array([2, 0, 0]))
            ]
    
    def test_invalid_init(self):
        with pytest.raises(TSValueError) as e:
            System(atoms=None,periodicity=False,cell=None,solvation=None,solvent=None,charge=0,electronic_energy=-33.6052447996,vibrational_frequencies=np.array([-1,-2,0.1,1,2,3,4,5,6]))
        assert str(e.value) == "Atoms must be provided."
        
        with pytest.raises(TSValueError) as e:
            System(atoms=self.atoms,periodicity=False,cell=None,solvation=None,solvent=None,charge=0,electronic_energy=-33.6052447996,vibrational_frequencies=None)
        assert str(e.value) == "Vibrational frequencies must be provided."
        
        with pytest.raises(TSValueError) as e:
            System(atoms=[],periodicity=False,cell=None,solvation=None,solvent=None,charge=0,electronic_energy=-33.6052447996,vibrational_frequencies=np.array([-1,-2,0.1,1,2,3,4,5,6]))
        assert str(e.value) == "The number of atoms must be greater than 0."
        
        
    def test_system(self):
        
            
        system = System(self.atoms,periodicity=False,cell=None,solvation=None,solvent=None,charge=0,electronic_energy=-33.6052447996,vibrational_frequencies=np.array([-1,-2,0.1,1,2,3,4,5,6]))

        np.testing.assert_array_equal(system.atom_names(), ['H', 'H', 'H'])
        
        np.testing.assert_array_equal(system.coord(),np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0]]))

        np.testing.assert_array_equal(system.atomic_masses()[0], [1.00794, 1.00794, 1.00794])

        assert len(system.atoms) == 3

        assert system.charge == 0

        assert system.spin == 0

        assert system.number_of_atoms == 3

        assert system.mass == 3.02382

        assert system.dim == 1

        assert system.dof == 3 * 3 - 5

        assert system.rotational_symmetry_number() == 1

        assert system.rotational_group == 'D*h'

        assert system._spacegroup_number is None

        assert system._spacegroup is None

        assert system.periodicity is False

        np.testing.assert_array_equal(system.center_of_mass,[1, 0, 0])

        assert system.pbc == [False, False, False]

        assert system.cell ==  None

        assert system.solvation == False

        assert system.solvent == ''

        assert system.electronic_energy == -33.6052447996

        np.testing.assert_array_equal(system.vibrational_frequencies,np.array([-1,-2,0.1,1,2,3,4,5,6]))

        np.testing.assert_array_equal(system.imaginary_frequencies,np.array([-1,-2]))

        np.testing.assert_array_equal(system._real_vibrational_frequencies,np.array([3,4,5,6]))

        assert system._has_imaginary_frequencies == True

        assert system._check_frequency_length == True

        np.testing.assert_array_equal(system.x(),np.array([0,1,2]))

        np.testing.assert_array_equal(system.y(),np.array([0,0,0]))

        np.testing.assert_array_equal(system.z(),np.array([0,0,0]))







if __name__ == '__main__':
    unittest.main()