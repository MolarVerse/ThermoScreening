import pytest

import os
import numpy as np
import ase.io as ase_io
from ase.build import molecule
from ThermoScreening import BASE_PATH
from ThermoScreening.calculator import Geoopt, Hessian, Modes
from ThermoScreening.calculator.dftbplus import dftb_3ob_parameters

# --------------------------------------------------------------------------- #


def test_dftb():
    # Test the DFTB+ calculator
    assert os.system("which dftb+ > /dev/null") == 0, "DFTB+ is not installed."
    assert os.system("which modes > /dev/null") == 0, "Modes is not installed."


# Test the Geoopt class
@pytest.mark.usefixtures("tmpdir")
def test_geoopt():
    """
    Test the Geoopt class.
    """

    # read the atoms object
    atoms = molecule("CH4")

    # create an instance of the Geoopt class
    geoopt = Geoopt(atoms, charge=0, **dftb_3ob_parameters)

    atoms.calc = geoopt

    assert geoopt.atoms == atoms
    assert np.allclose(geoopt.potential_energy(), -3.2295255, atol=1e-5)
    assert geoopt.label == "geo_opt"
    assert geoopt.slako_dir == BASE_PATH + "../external/slakos/3ob-3-1/"

    # read the optimized geometry
    geoopt.read()

    # check the optimized geometry
    assert geoopt.atoms == ase_io.read("geo_opt.gen", format="gen")
    assert atoms != geoopt.atoms
    assert atoms != ase_io.read("geo_opt.gen", format="gen")


@pytest.mark.usefixtures("tmpdir")
def test_hessian():
    """
    Test the Hessian class.
    """

    # read the atoms object
    atoms = molecule("CH4")

    # create an instance of the Hessian class
    optimizer = Geoopt(atoms, charge=0, **dftb_3ob_parameters)
    hessian = Hessian(optimizer.read(), charge=0, **dftb_3ob_parameters)

    assert hessian.atoms == optimizer.read()
    assert hessian.label == "second_derivative"
    assert hessian.slako_dir == BASE_PATH + "../external/slakos/3ob-3-1/"

    assert os.path.exists("hessian.out")
    assert np.allclose(hessian.read(), hessian.hessian, atol=1e-5)


@pytest.mark.usefixtures("tmpdir")
def test_modes():
    """
    Test the Modes class.
    """

    # read the atoms object
    atoms = molecule("CH4")

    # create an instance of the Modes class
    optimizer = Geoopt(atoms, charge=0, **dftb_3ob_parameters)
    hessian = Hessian(optimizer.read(), charge=0, **dftb_3ob_parameters)

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
        wave_numbers[8:],
        [
            1308.3501094,
            1485.00285228,
            1485.06847189,
            2858.13153344,
            3049.63204897,
            3050.31867544,
            3051.47541999,
        ],
        atol=1e-5,
    )
