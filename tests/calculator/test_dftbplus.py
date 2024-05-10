import pytest

import os
import numpy as np
import ase.io as ase_io
from ase.build import molecule
from ThermoScreening import BASE_PATH
from ThermoScreening.calculator import Geoopt, Hessian, Modes
from ThermoScreening.calculator.dftbplus import dftb_3ob_parameters
from ThermoScreening.thermo.api import dftbplus_thermo

# --------------------------------------------------------------------------- #

@pytest.mark.skip(reason="Issues with the DFTB+ calculator.")
class TestDftbplus:

    def test_dftb(self):
        # Test the DFTB+ calculator
        assert os.system("which dftb+ > /dev/null") == 0, "DFTB+ is not installed."
        assert os.system("which modes > /dev/null") == 0, "Modes is not installed."


    # Test the Geoopt class
    @pytest.mark.usefixtures("tmpdir")
    def test_geoopt(self):
        """
        Test the Geoopt class.
        """

        # read the atoms object
        atoms = molecule("H2O")

        # create an instance of the Geoopt class
        geoopt = Geoopt(atoms, charge=0, **dftb_3ob_parameters)

        atoms.calc = geoopt

        assert geoopt.atoms == atoms
        assert np.allclose(geoopt.potential_energy(), -4.06231229, atol=1e-5)
        assert geoopt.label == "geo_opt"
        assert geoopt.slako_dir == BASE_PATH + "../external/slakos/3ob-3-1/"

        # read the optimized geometry
        geoopt.read()

        # check the optimized geometry
        assert geoopt.atoms == ase_io.read("geo_opt.gen", format="gen")
        assert atoms != geoopt.atoms
        assert atoms != ase_io.read("geo_opt.gen", format="gen")


    @pytest.mark.usefixtures("tmpdir")
    def test_hessian(self):
        """
        Test the Hessian class.
        """

        # read the atoms object
        atoms = molecule("H2O")

        # create an instance of the Hessian class
        optimizer = Geoopt(atoms, charge=0, **dftb_3ob_parameters)
        hessian = Hessian(optimizer.read(), charge=0, **dftb_3ob_parameters, delta=0.001)

        assert hessian.atoms == optimizer.read()
        assert hessian.label == "second_derivative"
        assert hessian.slako_dir == BASE_PATH + "../external/slakos/3ob-3-1/"

        assert os.path.exists("hessian.out")
        assert np.allclose(hessian.read(), hessian.hessian, atol=1e-5)


    @pytest.mark.usefixtures("tmpdir")
    def test_modes(self):
        """
        Test the Modes class.
        """

        # read the atoms object
        atoms = molecule("H2O")

        # create an instance of the Modes class
        optimizer = Geoopt(atoms, charge=0, **dftb_3ob_parameters)
        hessian = Hessian(optimizer.read(), charge=0, **dftb_3ob_parameters, delta=0.001)

        assert os.path.exists("hessian.out")
        assert os.path.exists("geo_opt.gen")

        modes = Modes()

        assert modes.geometry == "geo_opt.gen"
        assert modes.hessian == "hessian.out"

        assert os.path.exists("modes_in.hsd")
        assert os.path.exists("modes.out")
        assert os.path.exists("vibrations.tag")

        wave_numbers = modes.read()

        assert np.allclose(
            wave_numbers,
            [
                -9.59971762e00,
                -8.58386111e00,
                -2.50130203e00,
                1.37799719e-01,
                4.48566805e00,
                8.42985145e00,
                1.46181960e03,
                3.60482906e03,
                3.87719263e03,
            ],
            atol=1e-5,
        )

    @pytest.mark.usefixtures("tmpdir")
    def test_dftbplus_thermo(self):
        atoms = molecule("CH4")
        thermo = dftbplus_thermo(atoms, **dftb_3ob_parameters)
        assert np.allclose(thermo.total_EeGtot(), -3.20470041120975, atol=1e-5)

        atoms = molecule("H2O")
        thermo = dftbplus_thermo(atoms, **dftb_3ob_parameters, delta=0.001)
        assert np.allclose(thermo.total_EeGtot(), -4.063547287431202, atol=1e-5)
