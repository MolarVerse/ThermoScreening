import unittest
import numpy as np
from ase.atoms import Atoms
from ThermoScreening.thermo.system import System
from ThermoScreening.thermo.atoms import Atom
from ThermoScreening.thermo.cell import Cell
import pytest




class TestSystem(unittest.TestCase):
    def test_system(self):
        atoms = [
            Atom(symbol='H', position=np.array([0, 0, 0])),
            Atom(symbol='H', position=np.array([1, 0, 0])),
            Atom(symbol='H', position=np.array([2, 0, 0]))
        ]

        system = System(atoms, periodicity=False, cell=None, solvation=None, solvent=None, charge=0,
                        electronic_energy=-33.6052447996, vibrational_frequencies=np.array([-1, -2, 0.1, 1, 2, 3, 4, 5, 6]))

        np.testing.assert_array_equal(system.atom_names(), ['H', 'H', 'H'])

        np.testing.assert_array_equal(
            system.coord(), np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0]]))

        np.testing.assert_array_equal(system.atomic_masses()[0], [
                                      1.00794, 1.00794, 1.00794])

        assert len(system.atoms) == 3

        assert system.charge == 0

        assert system.spin == 0

        assert system.number_of_atoms == 3

        assert system.mass == 3.02382

        assert system.dim == 1

        assert system.dof == 3 * 3 - 5

        assert system.rotational_symmetry_number() == 1

        assert system.rotational_group == 'D*h'

        assert system._spacegroup_number == None

        assert system._spacegroup == None

        assert system.periodicity == False

        np.testing.assert_array_equal(system.center_of_mass, [1, 0, 0])

        assert system.pbc == [False, False, False]

        assert system.cell == None

        assert system.solvation == False

        assert system.solvent == ''

        assert system.electronic_energy == -33.6052447996

        np.testing.assert_array_equal(
            system.vibrational_frequencies, np.array([-1, -2, 0.1, 1, 2, 3, 4, 5, 6]))

        np.testing.assert_array_equal(
            system.imaginary_frequencies, np.array([-1, -2]))

        np.testing.assert_array_equal(
            system._real_vibrational_frequencies, np.array([3, 4, 5, 6]))

        assert system._has_imaginary_frequencies == True

        assert system._check_frequency_length == True

        np.testing.assert_array_equal(system.x(), np.array([0, 1, 2]))

        np.testing.assert_array_equal(system.y(), np.array([0, 0, 0]))

        np.testing.assert_array_equal(system.z(), np.array([0, 0, 0]))

    def test_system_aq(self):

        atoms = [
            Atom(symbol='O', position=np.array(
                [0.00008112, 0.00000188, -2.05402286])),
            Atom(symbol='O', position=np.array(
                [-0.00009254, 0.00000502, -7.65886167])),
            Atom(symbol='C', position=np.array(
                [-3.68652488, -0.00000123, -5.56937006])),
            Atom(symbol='C', position=np.array(
                [-2.49045658, 0.00000104, -6.25271675])),
            Atom(symbol='C', position=np.array(
                [-1.22985769, 0.00000200, -5.58346285])),
            Atom(symbol='C', position=np.array(
                [-1.22981936, 0.00000095, -4.12932133])),
            Atom(symbol='C', position=np.array(
                [-2.49034390, -0.00000124, -3.45992559])),
            Atom(symbol='C', position=np.array(
                [-3.68647387, -0.00000247, -4.14315959])),
            Atom(symbol='C', position=np.array(
                [-0.00004268, 0.00000431, -6.32159897])),
            Atom(symbol='C', position=np.array(
                [0.00004092, 0.00000241, -3.39127855])),
            Atom(symbol='C', position=np.array(
                [1.22985982, 0.00000519, -4.12939595])),
            Atom(symbol='C', position=np.array(
                [1.22981336, 0.00000624, -5.58353992])),
            Atom(symbol='C', position=np.array(
                [2.49035536, 0.00001035, -6.25289256])),
            Atom(symbol='H', position=np.array(
                [2.48994553, 0.00001106, -7.34077798])),
            Atom(symbol='C', position=np.array(
                [3.68646897, 0.00001431, -5.56962708])),
            Atom(symbol='C', position=np.array(
                [3.68653147, 0.00001307, -4.14341821])),
            Atom(symbol='C', position=np.array(
                [2.49044502, 0.00000805, -3.46009695])),
            Atom(symbol='H', position=np.array(
                [-4.63230122, -0.00000198, -6.10807371])),
            Atom(symbol='H', position=np.array(
                [-2.49018198, 0.00000208, -7.34059729])),
            Atom(symbol='H', position=np.array(
                [-2.48990606, -0.00000205, -2.37203685])),
            Atom(symbol='H', position=np.array(
                [-4.63218409, -0.00000429, -3.60434535])),
            Atom(symbol='H', position=np.array(
                [4.63217171, 0.00001845, -6.10844479])),
            Atom(symbol='H', position=np.array(
                [4.63230600, 0.00001612, -3.60469800])),
            Atom(symbol='H', position=np.array(
                [2.49016578, 0.00000689, -2.37221336]))
        ]

        charge = -2

        vibrational_frequencies = np.array([-13.800,-2.54,6.36,25.33,46.98,53.790,90.3, 143.96, 207.59, 222.69, 229.17, 233.44, 281.5, 364.2, 388.97, 393.46, 418.36, 439.41, 449.72, 469.27, 513.78, 597.7, 606.64, 624.86, 629.79, 678.06, 686.75, 701.51, 717.21, 729.74, 734.03, 802.32, 807.81, 825.31, 896.02, 896.94, 898.48, 898.59, 959.91,
                                           960.79, 1043.79, 1047.42, 1056.91, 1120, 1127.16, 1130.62, 1131.3, 1192.85, 1202.44, 1275.69, 1321.66, 1323.16, 1333.16, 1397.04, 1433.18, 1504.85, 1514.35, 1518.2, 1543.43, 1596.63, 1625.67, 1677.6, 1699.25, 1701.29, 2979.33, 2980.88, 2986.11, 2986.64, 2995.14, 2995.94, 3000.44, 3001.33])

        elecronic_energy = -33.8387075830

        system = System(atoms, periodicity=False, cell=None, solvation=True, solvent="DMF", charge=charge,
                        electronic_energy=elecronic_energy, vibrational_frequencies=vibrational_frequencies)

        np.testing.assert_array_equal(system.atom_names(), [
                                      'O', 'O', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'H', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H'])

        np.testing.assert_array_equal(
            system.coord(), np.array([[0.00008112, 0.00000188, -2.05402286], [-0.00009254, 0.00000502, -7.65886167], [-3.68652488, -0.00000123, -5.56937006], [-2.49045658, 0.00000104, -6.25271675], [-1.22985769, 0.00000200, -5.58346285], [-1.22981936, 0.00000095, -4.12932133], [-2.49034390, -0.00000124, -3.45992559], [-3.68647387, -0.00000247, -4.14315959], [-0.00004268, 0.00000431, -6.32159897], [0.00004092, 0.00000241, -3.39127855], [1.22985982, 0.00000519, -4.12939595], [1.22981336, 0.00000624, -5.58353992], [2.49035536, 0.00001035, -6.25289256], [2.48994553, 0.00001106, -7.34077798], [3.68646897, 0.00001431, -5.56962708], [3.68653147, 0.00001307, -4.14341821], [2.49044502, 0.00000805, -3.46009695], [-4.63230122, -0.00000198, -6.10807371], [-2.49018198, 0.00000208, -7.34059729], [-2.48990606, -0.00000205, -2.37203685], [-4.63218409, -0.00000429, -3.60434535], [4.63217171, 0.00001845, -6.10844479], [4.63230600, 0.00001612, -3.60469800], [2.49016578, 0.00000689, -2.37221336]]))

        x = np.array([8.112000e-05 ,-9.254000e-05 ,-3.686525e+00 ,-2.490457e+00 ,-1.229858e+00 ,-1.229819e+00 ,-2.490344e+00 ,-3.686474e+00 ,-4.268000e-05 ,4.092000e-05 ,1.229860e+00 ,1.229813e+00 ,2.490355e+00 ,2.489946e+00 ,3.686469e+00 ,3.686531e+00 ,2.490445e+00 ,-4.632301e+00 ,-2.490182e+00 ,-2.489906e+00 ,-4.632184e+00 ,4.632172e+00 ,4.632306e+00 ,2.490166e+00])
        y= np.array([1.880000e-06 ,5.020000e-06 ,-1.230000e-06 ,1.040000e-06 ,2.000000e-06 ,9.500000e-07 ,-1.240000e-06 ,-2.470000e-06 ,4.310000e-06 ,2.410000e-06 ,5.190000e-06 ,6.240000e-06 ,1.035000e-05 ,1.106000e-05 ,1.431000e-05 ,1.307000e-05 ,8.050000e-06 ,-1.980000e-06 ,2.080000e-06 ,-2.050000e-06 ,-4.290000e-06 ,1.845000e-05 ,1.612000e-05 ,6.890000e-06])
        z = np.array([-2.054023e+00 ,-7.658862e+00 ,-5.569370e+00 ,-6.252717e+00 ,-5.583463e+00 ,-4.129321e+00 ,-3.459926e+00 ,-4.143160e+00 ,-6.321599e+00 ,-3.391279e+00 ,-4.129396e+00 ,-5.583540e+00 ,-6.252893e+00 ,-7.340778e+00 ,-5.569627e+00 ,-4.143418e+00 ,-3.460097e+00 ,-6.108074e+00 ,-7.340597e+00 ,-2.372037e+00 ,-3.604345e+00 ,-6.108445e+00 ,-3.604698e+00 ,-2.372213e+00])
        np.testing.assert_allclose(system.coord(), np.array([x, y, z]).T, atol=1e-6)
        np.testing.assert_allclose(system.x(), x, atol=1e-6)
        np.testing.assert_allclose(system.y(), y, atol=1e-6)
        np.testing.assert_allclose(system.z(), z, atol=1e-6)
        np.testing.assert_array_equal(system.atomic_masses()[0], [15.9994])
        np.testing.assert_array_equal(system.atomic_masses()[1], [15.9994])
        np.testing.assert_array_equal(system.atomic_masses()[2], [12.0107])
        np.testing.assert_array_equal(system.atomic_masses()[23], [1.00794])

        assert len(system.atoms) == 24

        assert system.charge == -2

        assert system.spin == 0

        assert system.number_of_atoms == 24

        assert system.mass == 208.21212

        assert system.dim == 3

        assert system.dof == 3 * 24 - 6

        assert system.rotational_symmetry_number() == 4

        assert system.rotational_group == 'D2h'

        assert system._spacegroup_number == None

        assert system._spacegroup == None

        assert system.periodicity == False

        # use ase to get the center of mass
        ase_atoms = Atoms(symbols=system.atom_names(), positions=system.coord())
        print(ase_atoms.get_center_of_mass())
        def center_of_mass(system):
                total_mass = np.sum(system.atomic_masses())
                center_of_mass = np.zeros(3)
                for i in range(system.number_of_atoms):
                        center_of_mass += system.coord()[i] * system.atomic_masses()[i]
                center_of_mass /= total_mass
                return center_of_mass
        np.testing.assert_allclose(ase_atoms.get_center_of_mass(), center_of_mass(system), atol=1e-10)
        np.testing.assert_array_equal(system.center_of_mass, center_of_mass(system))



        assert system.pbc == [False, False, False]

        assert system.cell == None

        assert system.solvation == True

        assert system.solvent == 'DMF'

        assert system.electronic_energy == -33.8387075830

        np.testing.assert_array_equal(
            system.vibrational_frequencies, np.array([-13.800,-2.54,6.36,25.33,46.98,53.790,90.3, 143.96, 207.59, 222.69, 229.17, 233.44, 281.5, 364.2, 388.97, 393.46, 418.36, 439.41, 449.72, 469.27, 513.78, 597.7, 606.64, 624.86, 629.79, 678.06, 686.75, 701.51, 717.21, 729.74, 734.03, 802.32, 807.81, 825.31, 896.02, 896.94, 898.48, 898.59, 959.91, 960.79, 1043.79, 1047.42, 1056.91, 1120, 1127.16, 1130.62, 1131.3, 1192.85, 1202.44, 1275.69, 1321.66, 1323.16, 1333.16, 1397.04, 1433.18, 1504.85, 1514.35, 1518.2, 1543.43, 1596.63, 1625.67, 1677.6, 1699.25, 1701.29, 2979.33, 2980.88, 2986.11, 2986.64, 2995.14, 2995.94, 3000.44, 3001.33]))

        np.testing.assert_array_equal(
            system.imaginary_frequencies, np.array([-13.800,-2.54]))

        np.testing.assert_array_equal(
            system._real_vibrational_frequencies, np.array([90.3, 143.96, 207.59, 222.69, 229.17, 233.44, 281.5, 364.2, 388.97, 393.46, 418.36, 439.41, 449.72, 469.27, 513.78, 597.7, 606.64, 624.86, 629.79, 678.06, 686.75, 701.51, 717.21, 729.74, 734.03, 802.32, 807.81, 825.31, 896.02, 896.94, 898.48, 898.59, 959.91, 960.79, 1043.79, 1047.42, 1056.91, 1120, 1127.16, 1130.62, 1131.3, 1192.85, 1202.44, 1275.69, 1321.66, 1323.16, 1333.16, 1397.04, 1433.18, 1504.85, 1514.35, 1518.2, 1543.43, 1596.63, 1625.67, 1677.6, 1699.25, 1701.29, 2979.33, 2980.88, 2986.11, 2986.64, 2995.14, 2995.94, 3000.44, 3001.33]))

        assert system._has_imaginary_frequencies == True

        assert system._check_frequency_length == True



if __name__ == '__main__':
    unittest.main()
