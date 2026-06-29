import pytest
import os

import unittest
import numpy as np

from ThermoScreening.thermo.system import System
from ThermoScreening.thermo.atoms import Atom
from ThermoScreening.thermo.cell import Cell
from ThermoScreening.thermo.thermo import Thermo, _real_scalar
from ThermoScreening.exceptions import TSValueError
from ThermoScreening.utils.physicalConstants import PhysicalConstants


def _thermo_with_scalar_totals():
    thermo = object.__new__(Thermo)
    thermo._total_energy_Hartree_per_mol = np.complex128(1.0 + 0.0j)
    thermo._total_energy_kcal = np.complex128(2.0 + 0.0j)
    thermo._total_energy = np.complex128(3.0 + 0.0j)
    thermo._total_enthalpy_Hartree_per_mol = np.complex128(4.0 + 0.0j)
    thermo._total_enthalpy_kcal = np.complex128(5.0 + 0.0j)
    thermo._total_enthalpy = np.complex128(6.0 + 0.0j)
    thermo._total_entropy_Hartree_per_mol = np.complex128(7.0 + 0.0j)
    thermo._total_entropy = np.complex128(8.0 + 0.0j)
    thermo._total_gibbs_free_energy_Hartree_per_mol = np.complex128(9.0 + 0.0j)
    thermo._total_gibbs_free_energy_kcal = np.complex128(10.0 + 0.0j)
    thermo._total_gibbs_free_energy = np.complex128(11.0 + 0.0j)
    thermo._total_heatcapacity = np.complex128(12.0 + 0.0j)
    thermo._EeZPE = np.complex128(13.0 + 0.0j)
    thermo._EeEtot = np.complex128(14.0 + 0.0j)
    thermo._EeHtot = np.complex128(15.0 + 0.0j)
    thermo._EeGtot = np.complex128(16.0 + 0.0j)
    thermo._system = type("SystemStub", (), {
        "electronic_energy": np.complex128(17.0 + 0.0j)
    })()
    return thermo


def test_real_scalar_accepts_zero_imaginary_complex():
    assert _real_scalar(np.complex128(1.25 + 0.0j)) == 1.25
    assert _real_scalar(np.complex128(1.25 + 1e-12j)) == pytest.approx(1.25)


def test_real_scalar_rejects_nonzero_imaginary_complex():
    with pytest.raises(TypeError):
        _real_scalar(np.complex128(1.25 + 0.5j))


def test_scalar_total_accessors_return_real_floats_for_all_units():
    thermo = _thermo_with_scalar_totals()

    assert thermo.total_energy("H") == 1.0
    assert thermo.total_energy("kcal") == 2.0
    assert thermo.total_energy("cal") == 3.0
    assert thermo.total_enthalpy("H") == 4.0
    assert thermo.total_enthalpy("kcal") == 5.0
    assert thermo.total_enthalpy("cal") == 6.0
    assert thermo.total_entropy("H/T") == 7.0
    assert thermo.total_entropy("cal/(mol*K)") == 8.0
    assert thermo.total_gibbs_free_energy("H") == 9.0
    assert thermo.total_gibbs_free_energy("kcal") == 10.0
    assert thermo.total_gibbs_free_energy("cal") == 11.0
    assert thermo.total_heat_capacity("cal/(mol*K)") == 12.0
    assert thermo.total_heat_capacity("H/T") == pytest.approx(
        12.0 * PhysicalConstants["cal"] / PhysicalConstants["H"] / PhysicalConstants["N_A"]
    )
    assert thermo.total_EeZPE() == 13.0
    assert thermo.total_EeEtot() == 14.0
    assert thermo.total_EeHtot() == 15.0
    assert thermo.total_EeGtot() == 16.0
    assert thermo.electronic_energy() == 17.0


@pytest.mark.parametrize(
    "method_name",
    [
        "total_energy",
        "total_enthalpy",
        "total_entropy",
        "total_gibbs_free_energy",
        "total_heat_capacity",
    ],
)
def test_scalar_total_accessors_reject_unknown_units(method_name):
    thermo = _thermo_with_scalar_totals()

    with pytest.raises(TSValueError):
        getattr(thermo, method_name)("unknown")


def _valid_system():
    atoms = [Atom(symbol="H", position=np.array([i, 0.0, 0.0])) for i in range(3)]
    return System(
        atoms,
        periodicity=False,
        cell=None,
        charge=0,
        electronic_energy=-1.0,
        vibrational_frequencies=np.array([-1, -2, 0.1, 1, 2, 3, 4, 5, 6]),
    )


def test_thermo_rejects_unsupported_engine():
    with pytest.raises(TSValueError, match="engine is not supported"):
        Thermo(temperature=298.15, pressure=101325, system=_valid_system(), engine="xtb")


def test_thermo_rejects_negative_temperature():
    with pytest.raises(TSValueError, match="temperature is negative"):
        Thermo(temperature=-1.0, pressure=101325, system=_valid_system(), engine="dftb+")


def test_thermo_rejects_negative_pressure():
    with pytest.raises(TSValueError, match="pressure is negative"):
        Thermo(temperature=298.15, pressure=-1.0, system=_valid_system(), engine="dftb+")


# --- geometry-dependent thermochemistry, validated against ASE IdealGasThermo --- #

from ase import Atoms  # noqa: E402
from ase.thermochemistry import IdealGasThermo  # noqa: E402

_CM_TO_EV = 1.23984198e-4
_EVK_TO_CALMOLK = 1.602176634e-19 * 6.02214076e23 / 4.184
_T, _P = 298.15, 101325.0
_R_CALMOLK = 8.314462618 / 4.184


def _ts_thermo(symbols, positions, real_freqs, dof):
    atoms = [Atom(symbol=s, position=np.array(p, float)) for s, p in zip(symbols, positions)]
    pad = np.concatenate([np.zeros(3 * len(atoms) - dof), np.asarray(real_freqs, float)])
    system = System(
        atoms, periodicity=False, cell=None, charge=0,
        electronic_energy=0.0, vibrational_frequencies=pad,
    )
    thermo = Thermo(temperature=_T, pressure=_P, system=system, engine="dftb+")
    thermo.run()
    return thermo


def _ase_total_entropy(symbols, positions, real_freqs, geometry, sigma):
    ase_thermo = IdealGasThermo(
        vib_energies=[f * _CM_TO_EV for f in real_freqs],
        geometry=geometry,
        atoms=Atoms(symbols=symbols, positions=positions),
        symmetrynumber=sigma, spin=0, potentialenergy=0.0,
    )
    return ase_thermo.get_entropy(_T, _P, verbose=False) * _EVK_TO_CALMOLK


@pytest.mark.parametrize(
    "symbols,positions,freqs,dof,geometry,sigma",
    [
        (["O", "H", "H"], [[0, 0, 0.119], [0, 0.763, -0.477], [0, -0.763, -0.477]],
         [1595.0, 3657.0, 3756.0], 3, "nonlinear", 2),
        (["C", "O", "O"], [[0, 0, 0], [0, 0, 1.16], [0, 0, -1.16]],
         [667.0, 667.0, 1333.0, 2349.0], 4, "linear", 2),
        (["N", "N"], [[0, 0, 0], [0, 0, 1.10]], [2359.0], 1, "linear", 2),
        (["Ar"], [[0, 0, 0]], [], 0, "monatomic", 1),
    ],
)
def test_total_entropy_matches_ase(symbols, positions, freqs, dof, geometry, sigma):
    thermo = _ts_thermo(symbols, positions, freqs, dof)
    ts_total = thermo.total_entropy("cal/(mol*K)")
    ase_total = _ase_total_entropy(symbols, positions, freqs, geometry, sigma)

    assert np.isfinite(ts_total)
    assert ts_total == pytest.approx(ase_total, abs=0.05)


def test_thermo_rejects_imaginary_vibrational_mode():
    # an imaginary (negative) mode in the kept dof set used to give NaN; it now
    # raises (the geometry is not a minimum)
    with pytest.raises(TSValueError, match="[Ii]maginary"):
        _ts_thermo(
            ["O", "H", "H"],
            [[0, 0, 0.119], [0, 0.763, -0.477], [0, -0.763, -0.477]],
            [-500.0, 3657.0, 3756.0],
            3,
        )


def test_rotational_contribution_handles_linear_and_monatomic():
    # linear and monatomic species used to give -inf rotational entropy
    co2 = _ts_thermo(["C", "O", "O"], [[0, 0, 0], [0, 0, 1.16], [0, 0, -1.16]],
                     [667.0, 667.0, 1333.0, 2349.0], 4)
    assert co2._n_rot == 2
    assert np.isfinite(co2._rotational_entropy)
    assert co2._rotational_heat_capacity == pytest.approx(_R_CALMOLK)  # Cv_rot = R (linear)

    argon = _ts_thermo(["Ar"], [[0, 0, 0]], [], 0)
    assert argon._n_rot == 0
    assert argon._rotational_entropy == 0.0
    assert argon._rotational_heat_capacity == 0.0


class TestThermo(unittest.TestCase):

    def test_thermo_aq_2(self):
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

        vibrational_frequencies = np.array([-13.800, -2.54, 6.36, 25.33, 46.98, 53.790, 90.3, 143.96, 207.59, 222.69, 229.17, 233.44, 281.5, 364.2, 388.97, 393.46, 418.36, 439.41, 449.72, 469.27, 513.78, 597.7, 606.64, 624.86, 629.79, 678.06, 686.75, 701.51, 717.21, 729.74, 734.03, 802.32, 807.81, 825.31, 896.02, 896.94, 898.48, 898.59, 959.91,
                                            960.79, 1043.79, 1047.42, 1056.91, 1120, 1127.16, 1130.62, 1131.3, 1192.85, 1202.44, 1275.69, 1321.66, 1323.16, 1333.16, 1397.04, 1433.18, 1504.85, 1514.35, 1518.2, 1543.43, 1596.63, 1625.67, 1677.6, 1699.25, 1701.29, 2979.33, 2980.88, 2986.11, 2986.64, 2995.14, 2995.94, 3000.44, 3001.33])

        elecronic_energy = -33.8387075830

        system = System(atoms, periodicity=False, cell=None, solvation=True, solvent="DMF", charge=charge,
                        electronic_energy=elecronic_energy, vibrational_frequencies=vibrational_frequencies)
        thermo = Thermo(system=system, temperature=298.15,
                        pressure=101325, engine="dftb+")

        assert thermo._temperature == 298.15

        assert thermo._pressure == 101325

        assert thermo._engine == "dftb+"

        assert thermo.electronic_energy() == -33.8387075830

        thermo.run()


        x_coord_relocated = np.array([8.215472e-05 ,-9.150528e-05 ,-3.686524e+00 ,-2.490456e+00 ,-1.229857e+00 ,-1.229818e+00 ,-2.490343e+00 ,-3.686473e+00 ,-4.164528e-05 ,4.195472e-05 ,1.229861e+00 ,1.229814e+00 ,2.490356e+00 ,2.489947e+00 ,3.686470e+00 ,3.686533e+00 ,2.490446e+00 ,-4.632300e+00 ,-2.490181e+00 ,-2.489905e+00 ,-4.632183e+00 ,4.632173e+00 ,4.632307e+00 ,2.490167e+0])

        y_coord_relocated = np.array([-2.507243e-06 ,6.327568e-07 ,-5.617243e-06 ,-3.347243e-06 ,-2.387243e-06 ,-3.437243e-06 ,-5.627243e-06 ,-6.857243e-06 ,-7.724321e-08 ,-1.977243e-06 ,8.027568e-07 ,1.852757e-06 ,5.962757e-06 ,6.672757e-06 ,9.922757e-06 ,8.682757e-06 ,3.662757e-06 ,-6.367243e-06 ,-2.307243e-06 ,-6.437243e-06 ,-8.677243e-06 ,1.406276e-05 ,1.173276e-05 ,2.502757e-06])

        z_coord_relocated = np.array([2.802395e+00 ,-2.802443e+00 ,-7.129518e-01 ,-1.396299e+00 ,-7.270446e-01 ,7.270969e-01 ,1.396493e+00 ,7.132586e-01 ,-1.465181e+00 ,1.465140e+00 ,7.270223e-01 ,-7.271217e-01 ,-1.396474e+00 ,-2.484360e+00 ,-7.132089e-01 ,7.130000e-01 ,1.396321e+00 ,-1.251655e+00 ,-2.484179e+00 ,2.484381e+00 ,1.252073e+00 ,-1.252027e+00 ,1.251720e+00 ,2.484205e+00]) 

       
        np.testing.assert_allclose(
            thermo._reloc_coord[:, 0], x_coord_relocated, atol=1e-6)

        np.testing.assert_allclose(
            thermo._reloc_coord[:, 1], y_coord_relocated, atol=1e-6)

        np.testing.assert_allclose(
            thermo._reloc_coord[:, 2], z_coord_relocated, atol=1e-6)

        #### ROTATIONAL PARTITION FUNCTION ####

        np.testing.assert_allclose(thermo._inertia_tensor, np.array([[477.579476, -0.002293, 0.023808], [-0.002293, 1612.635725, 0.000317],[0.023808, 0.000317, 1135.056249]]), atol=1e-4)

        # somewhere in the array the values should be the same
        np.testing.assert_allclose(np.sort(thermo._eigenvalues_I), np.sort(
            np.array([1612.635724886585, 477.579474789699, 1135.056250097402])), atol=1e-10)

        np.testing.assert_allclose(np.sort(thermo._rotational_temperature), np.sort(
            np.array([0.015040202022, 0.050786033259, 0.021368427414])), atol=1e-10)

        np.testing.assert_allclose(np.sort(thermo._rotational_constant), np.sort(
            np.array([0.313386961074, 1.058209231802, 0.445245783280])), atol=1e-10)

        np.testing.assert_allclose(
            thermo._rotational_temperature_xyz, 0.000016321893, atol=1e-10)

        np.testing.assert_allclose(
            thermo._rotational_partition_function, 564653.342499404098, atol=1e-10)

        np.testing.assert_allclose(
            thermo._rotational_entropy, 29.299274545140, atol=1e-10)

        np.testing.assert_allclose(
            thermo._rotational_heat_capacity, 2.980806387906, atol=1e-10)

        E = (3/2)*8.314462618 * 298.15 / 4.184 * \
            (4.184/(4.3597447222071 * 10**(-18))/(6.02214076 * 10**(23)))
        # print(E)
        np.testing.assert_allclose(thermo._rotational_energy*4.184/(
            4.3597447222071 * 10**(-18))/(6.02214076 * 10**(23)),  0.001416277301, atol=1e-10)

        #### VIBRATIONAL PARTITION FUNCTION ####

        np.testing.assert_allclose(
            thermo._vibrational_partition_function,  0.000000000000, atol=1e-10)

        np.testing.assert_allclose(
            thermo._vibrational_entropy, 28.264590256379, atol=1e-10)

        np.testing.assert_allclose(
            thermo._vibrational_heat_capacity, 41.364292152340, atol=1e-10)

        vib_real = np.array([90.3, 143.96, 207.59, 222.69, 229.17, 233.44, 281.5, 364.2, 388.97, 393.46, 418.36, 439.41, 449.72, 469.27, 513.78, 597.7, 606.64, 624.86, 629.79, 678.06, 686.75, 701.51, 717.21, 729.74, 734.03, 802.32, 807.81, 825.31, 896.02, 896.94, 898.48, 898.59, 959.91, 960.79, 1043.79, 1047.42, 1056.91, 1120, 1127.16, 1130.62, 1131.3, 1192.85, 1202.44, 1275.69, 1321.66, 1323.16, 1333.16, 1397.04, 1433.18, 1504.85, 1514.35, 1518.2, 1543.43, 1596.63, 1625.67, 1677.6, 1699.25, 1701.29, 2979.33, 2980.88, 2986.11, 2986.64, 2995.14, 2995.94, 3000.44, 3001.33])
        vib_temp = 6.62607015 * 10**(-34) * vib_real * 299792458 * 10**2 /( 1.380649 * 10**(-23))
        vib_temp2 = np.array([129.921552038605, 207.126319285466, 298.675692001042, 320.401222851351, 329.724497017576, 335.868074284518, 405.015691017357, 524.002538786933, 559.641042042705, 566.101150222698, 601.926694472546, 632.212947744003, 647.046737351069, 675.174825306271, 739.214784123971, 859.956939684101, 872.819604968986, 899.034119677108, 906.127289683202, 975.577049560317, 988.080020625826, 1009.316367337784, 1031.905164314596, 1049.933038589721, 1056.105391394212, 1154.359464358956, 1162.258349416452, 1187.436944772772, 1289.172857781075, 1290.496532508378, 1292.712248899734, 1292.870514356260, 1381.096312484801, 1382.362436137004, 1501.780916969831, 1507.003677035170, 1520.657669602683, 1611.430102804406, 1621.731745247334, 1626.709913243498, 1627.688281520200, 1716.244998330567, 1730.042868585830, 1835.433274862993, 1901.573847921849, 1903.732013238105, 1918.119782013144, 2010.028848948095, 2062.026245301088, 2165.143384111795, 2178.811764448082, 2184.351055426472, 2220.651396045896, 2297.194325929106, 2338.976406451820, 2413.692089700599, 2444.841609098559, 2447.776713928668, 4286.591114453795, 4288.821218613926, 4296.346021683272, 4297.108573428348, 4309.338176887132, 4310.489198389136, 4316.963694337903, 4318.244205758881])
        np.testing.assert_allclose(vib_temp, vib_temp2, atol=1e-10)
        np.testing.assert_allclose(thermo._vib_temp_K, vib_temp2, atol=1e-10)
        energy_vib =  (8.314462618 * np.sum(vib_temp * ( 0.5 + (1/( np.exp(vib_temp/298.15)-1)))))  
        energy_vib2 =  (8.314462618 * np.sum(vib_temp2 * ( 0.5 + (1/( np.exp(vib_temp2/298.15)-1)))))  


        np.testing.assert_allclose(energy_vib2,energy_vib, atol=1e-10)

        np.testing.assert_allclose(energy_vib2/(
            4.3597447222071 * 10**(-18))/(6.02214076 * 10**(23)),  0.181332151330, atol=1e-10)
      
        np.testing.assert_allclose(thermo._vibrational_energy*4.184,energy_vib, atol=1e-10)

        # # np.testing.assert_allclose(
        # #     thermo._vibrational_temperature, np.prod(vib_temp), atol=1e-6)

        np.testing.assert_allclose(thermo._vibrational_energy*4.184/(
            4.3597447222071 * 10**(-18))/(6.02214076 * 10**(23)),  0.181332151330, atol=1e-10)

        np.testing.assert_allclose(thermo._EZP*4.184/(4.3597447222071 * 10**(-18))/(
            6.02214076 * 10**(23)), 0.173008013562, atol=1e-10)

        # np.testing.assert_allclose(thermo._zpecorr, x, atol=1e-6)

        #### TRANSLATIONAL PARTITION FUNCTION ####

        np.testing.assert_allclose(
            thermo._translational_partition_function,  118089966.152834132314, atol=1e-10)

        np.testing.assert_allclose(
            thermo._translational_entropy, 41.904091381259, atol=1e-10)

        np.testing.assert_allclose(
            thermo._translational_heat_capacity, 2.980806387906, atol=1e-10)

        np.testing.assert_allclose(thermo._translational_energy*4.184/(
            4.3597447222071 * 10**(-18))/(6.02214076 * 10**(23)), 0.001416277301, atol=1e-10)

        #### ELECTRONIC PARTITION FUNCTION ####

        np.testing.assert_allclose(
            thermo._electronic_partition_function, 1.000000000000, atol=1e-10)

        np.testing.assert_allclose(
            thermo._electronic_entropy, 0.000000000000, atol=1e-10)

        np.testing.assert_allclose(
            thermo._electronic_heat_capacity, 0.000000000000, atol=1e-10)

        np.testing.assert_allclose(thermo._electronic_energy*4.184/(
            4.3597447222071 * 10**(-18))/(6.02214076 * 10**(23)), 0.000000000000, atol=1e-10)

        #### SUMMARY ###

        np.testing.assert_allclose(
            thermo.total_energy("H"), 0.184164705933, atol=1e-10)

        np.testing.assert_allclose(thermo.total_enthalpy(
            "H"),  0.185108890800, atol=1e-10)

        np.testing.assert_allclose(thermo.total_entropy(
            "cal/(mol*K)"), 99.467956182779, atol=1e-10)

        np.testing.assert_allclose(thermo.total_heat_capacity(
            "cal/(mol*K)"), 47.325904928153, atol=1e-10)

        np.testing.assert_allclose(
            thermo.total_gibbs_free_energy("H"), 0.137848455123, atol=1e-10)

       
        np.testing.assert_allclose(
            thermo.total_EeZPE(), -33.665699569438, atol=1e-10)

        np.testing.assert_allclose(
            thermo.total_EeEtot(), -33.654542877067, atol=1e-10)

        np.testing.assert_allclose(
            thermo.total_EeHtot(), -33.653598692200, atol=1e-10)

        np.testing.assert_allclose(
            thermo.total_EeGtot(), -33.700859127877, atol=1e-10)

        # difference by different masses and vibrations


if __name__ == '__main__':
    unittest.main()
