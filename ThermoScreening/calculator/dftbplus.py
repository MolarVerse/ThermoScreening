import os
from ase.calculators.dftb import Dftb
from ase.io import read
import numpy as np

from ..utils.physicalConstants import PhysicalConstants
from ThermoScreening import BASE_PATH

# --------------------------------------------------------------------------- #


class Geoopt(Dftb):
    """
    Custom DFTB+ calculator to optimize the system with LBFGS solver.
    It is a subclass of ase.calculators.dftb.Dftb. It uses the
    'LBFGS' driver.

    Parameters:
    -----------
    atoms : ase.Atoms
        Atoms object.
    label : str
        Label for the calculation. Default is "geo_opt".
    charge : int
        Charge of the system. Default is 0.
    slako_dir : str
        Path to the Slater-Koster files. If None, it will look for
        the DFTB_PREFIX environment variable.
    max_force : float
        Maximum force component. Default is 1.0e-6.

    Other Parameters:
    -----------------
    **kwargs : dict
        Additional keyword arguments to pass to the Dftb class.
    """

    def __init__(
        self,
        atoms,
        label="geo_opt",
        charge=0,
        slako_dir=None,
        max_force=1.0e-6,
        **kwargs,
    ):
        """
        Constructor of the Geoopt class.

        Parameters:
        -----------
        atoms : ase.Atoms
            Atoms object.
        label : str
            Label for the calculation. Default is 'geo_opt'.
        charge : int
            Charge of the system. Default is 0.
        slako_dir : str
            Path to the Slater-Koster files. If None, it will look
            for the DFTB_PREFIX environment variable. Default is 3ob-3-1.
        max_force : float
            Maximum force component. Default is 1.0e-6.

        Other Parameters:
        -----------------
        **kwargs : dict
            Additional keyword arguments to pass to the Dftb
        """

        super().__init__(
            atoms=atoms,
            label=label,
            slako_dir=BASE_PATH + "../external/slakos/3ob-3-1/",
            Hamiltonian_Charge=charge,
            Driver_="LBFGS",
            Driver_MaxForceComponent=max_force,
            Driver_OutputPrefix="geo_opt",
            **kwargs,
        )

        self.calculate(atoms)

        return None

    def potential_energy(self):
        """
        Get the potential energy of the optimized geometry. in Hartree.

        Returns:
        --------
        float
            Potential energy of the optimized geometry.
        """
        self.atoms.calc = self
        return (
            self.atoms.get_potential_energy()
            * PhysicalConstants["eV"]
            / PhysicalConstants["H"]
        )

    def read(self):
        """
        Read the optimized geometry from the "geo_opt.gen" file.

        Returns:
        --------
        ase.Atoms
            Atoms object of the optimized geometry.
        """

        self.atoms = read("geo_opt.gen", format="gen")
        return self.atoms


# --------------------------------------------------------------------------- #


class Hessian(Dftb):
    """
    Custom DFTB+ calculator to calculate the Hessian matrix.
    It is a subclass of ase.calculators.dftb.Dftb. It uses the
    'SecondDerivatives' driver to calculate the Hessian matrix.

    Parameters:
    -----------
    atoms : ase.Atoms
        Atoms object.
    label : str
        Label for the calculation. Default is "hessian".
    charge : int
        Charge of the system. Default is 0.
    slako_dir : str
        Path to the Slater-Koster files. If None, it will look for
        the DFTB_PREFIX environment variable.

    Other Parameters:
    -----------------
    **kwargs : dict
        Additional keyword arguments to pass to the Dftb class.
    """

    def __init__(
        self,
        atoms,
        label="second_derivative",
        charge=0,
        delta=1.0e-4,
        slako_dir=BASE_PATH + "../external/slakos/3ob-3-1/",
        **kwargs,
    ):
        """
        Constructor of the Hessian class.

        Parameters:
        -----------
        atoms : ase.Atoms
            Atoms object.
        label : str
            Label for the calculation. Default is 'hessian'.
        charge : int
            Charge of the system. Default is 0.
        delta : float
            Finite difference step. Default is 1.0e-4.
        slako_dir : str
            Path to the Slater-Koster files. If None, it will look
            for the DFTB_PREFIX environment variable. Default is 3ob-3-1.

        Other Parameters:
        -----------------
        **kwargs : dict
            Additional keyword arguments to pass to the Dftb
        """

        super().__init__(
            atoms=atoms,
            label=label,
            slako_dir=slako_dir,
            Hamiltonian_Charge=charge,
            Driver_="SecondDerivatives",
            Driver_Delta=delta,
            **kwargs,
        )

        self.calculate(atoms)

        return None

    def read(self):
        """
        Read the Hessian matrix from the 'hessian.out' file.

        Returns:
        --------
        numpy.ndarray
            Hessian matrix.
        """

        with open("hessian.out") as f:
            lines = [line.split() for line in f]

        # matrix to array
        hessian = []
        for line in lines:
            hessian += line
        hessian = np.array(hessian, dtype=float)

        self.hessian = hessian.reshape(
            self.atoms.get_global_number_of_atoms() * 3,
            self.atoms.get_global_number_of_atoms() * 3,
        )

        return self.hessian


# --------------------------------------------------------------------------- #


class Modes:
    """
    Class to calculate the vibrational modes from the Hessian matrix.

    Parameters:
    -----------
    geometry : str
        Path to the geometry file. Default is 'geo_opt.gen'.
    hessian : str
        Path to the Hessian matrix file. Default is 'hessian.out'.
    wave_numbers : numpy.ndarray
        Array of wave numbers.
    """

    def __init__(self, geometry="geo_opt.gen", hessian="hessian.out"):
        """
        Constructor of the Modes class.

        Parameters:
        -----------
        geometry : str
            Path to the geometry file. Default is 'geo_opt.gen'.
        hessian : str
            Path to the Hessian matrix file. Default is 'hessian.out'.
        """

        self.geometry = geometry
        self.hessian = hessian

        # write the modes_in.hsd file
        self.write()
        # calculate the vibrational modes
        self.calculate()
        # read the vibrational modes
        self.wave_numbers = self.read()

        return None

    def write(self):
        """
        Writes the modes_in.hsd file.

        Returns:
        --------
        None
        """

        # use f-strings to write the modes_in.hsd file
        string = (
            "Geometry = GenFormat {\n"
            f"    <<< {self.geometry}\n"
            "}\n"
            "\n"
            "Hessian = {\n"
            f"    <<< {self.hessian}\n"
            "}\n"
            "\n"
            "Atoms = 1:-1\n"
            "\n"
        )

        with open("modes_in.hsd", "w") as f:
            f.write(string)

        return None

    def calculate(self):
        """
        Calculate the vibrational modes using the modes executable from DFTB+.

        Throws:
        -------
        FileNotFoundError
            If the modes executable is not found.
        FileNotFoundError
            If the modes_in.hsd file is not found.

        Returns:
        --------
        None
        """

        # run the modes executable
        os.system("modes > modes.out")

        return None

    def read(self):
        """
        Read the vibrational modes from the vibrations.tag file.

        Returns:
        --------
        numpy.ndarray
            Vibrational modes.
        """

        with open("vibrations.tag") as f:
            f.readline()  # skip the first line
            lines = [line.split() for line in f]

        # matrix to array
        modes = []
        for line in lines:
            modes += line
        modes = np.array(modes, dtype=float)

        # Hartree to cm^-1 - 1 Hartree = 219474.63 cm^-1
        self.wave_numbers = modes * 219474.63

        return self.wave_numbers


# --------------------------------------------------------------------------- #

"""
DFTB+ parameters for the 3ob-3-1 Slater-Koster files.
"""

dftb_3ob_parameters = dict(
    # SCC 
    Hamiltonian_SCC="Yes",
    Hamiltonian_MaxSCCIterations=250,
    Hamiltonian_SCCTolerance='1.0e-7',
    Hamiltonian_ReadInitialCharges="No",

    # Fermi smearing
    Hamiltonian_Filling="Fermi {",
    Hamiltonian_Filling_empty="Temperature [Kelvin] = 300",

    # Convergence helper
    Hamiltonian_Mixer="DIIS{}",

    # Are guessed by ase
    Hamiltonian_MaxAngularMomentum_="",

    # Hubbard derivatives
    Hamiltonian_ThirdOrderFull="Yes",
    Hamiltonian_hubbardderivs_="",
    Hamiltonian_hubbardderivs_C=-0.1492,
    Hamiltonian_hubbardderivs_N=-0.1535,
    Hamiltonian_hubbardderivs_O=-0.1575,
    Hamiltonian_hubbardderivs_H=-0.1857,
    Hamiltonian_hubbardderivs_S=-0.11,
    Hamiltonian_hubbardderivs_P=-0.14,
    Hamiltonian_hubbardderivs_F=-0.1623,
    Hamiltonian_hubbardderivs_Cl=-0.0697,
    Hamiltonian_hubbardderivs_Br=-0.0573,
    Hamiltonian_hubbardderivs_I=-0.0433,
    Hamiltonian_hubbardderivs_Zn=-0.03,
    Hamiltonian_hubbardderivs_Mg=-0.02,
    Hamiltonian_hubbardderivs_Ca=-0.0340,
    Hamiltonian_hubbardderivs_K=-0.0339,
    Hamiltonian_hubbardderivs_Na=-0.0454,

    # Analysis
    Analysis_="",
    Analysis_CalculateForces="Yes",
    Analysis_MullikenAnalysis="Yes",
)