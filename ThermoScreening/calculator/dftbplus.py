"""DFTB+ calculators: geometry optimisation, Hessian, and vibrational modes."""

import os
import shutil
import subprocess

import numpy as np
from ase.calculators.dftb import Dftb
from ase.io import read

from ..utils.physicalConstants import PhysicalConstants
from ..exceptions import TSValueError

# --------------------------------------------------------------------------- #


def _read_hessian_matrix(filename, size):
    """
    Read a DFTB+ Hessian matrix from ``filename``.

    DFTB+ writes ``hessian.out`` as a flat stream of second derivatives wrapped
    across a fixed number of values per line (with a ragged final line per
    matrix row), so the matrix is reconstructed by reading every value and
    reshaping to ``(size, size)`` rather than treating each physical line as a
    matrix row.
    """
    with open(filename, "r", encoding="utf-8") as handle:
        values = np.array(handle.read().split(), dtype=float)

    if values.size != size * size:
        raise ValueError(
            "Hessian matrix size does not match the number of atoms."
        )

    return values.reshape(size, size)


def _slako_dir(slako_dir=None):
    selected_dir = slako_dir or os.getenv("DFTB_PREFIX")
    if not selected_dir:
        raise FileNotFoundError(
            "Slater-Koster files are not bundled with ThermoScreening. "
            "Set DFTB_PREFIX or pass slako_dir to the DFTB+ calculator."
        )

    selected_dir = os.path.abspath(os.path.expanduser(selected_dir))
    if not os.path.isdir(selected_dir):
        raise FileNotFoundError(
            f"Slater-Koster directory does not exist: {selected_dir}"
        )

    return selected_dir + os.sep


# Atomic spin constants (Hartree): the spin constant of the highest occupied
# shell per element (Wss for H, Wpp for the p-block, valence Wss for the s-block
# and Zn), used with ShellResolvedSpin = No to match the atom-resolved SCC. These
# are parameters tied to the Slater-Koster set + functional, so each parameter set
# has its own.
#
# 3ob-3-1 (PBE): taken from the authoritative ``spinw.hsd`` shipped with the set
# (calculated with PBE/slateratom) and match it exactly. Verified end-to-end: an
# OH radical runs spin-polarised (0.40 eV below restricted, S_elec = R ln 2).
SPIN_CONSTANTS_3OB = {
    "H": "{ -0.07174 }",
    "C": "{ -0.02265 }",
    "N": "{ -0.02545 }",
    "O": "{ -0.02785 }",
    "F": "{ -0.02990 }",
    "Na": "{ -0.01528 }",
    "Mg": "{ -0.01667 }",
    "P": "{ -0.01490 }",
    "S": "{ -0.01549 }",
    "Cl": "{ -0.01606 }",
    "K": "{ -0.01075 }",
    "Ca": "{ -0.01196 }",
    "Zn": "{ -0.01680 }",
    "Br": "{ -0.01377 }",
    "I": "{ -0.01144 }",
}

# mio-1-1 reuses the 3ob spin constants. mio is the LDA-based Elstner-1998 set
# (PRB 58, 7260) and ships no spin constants of its own; the only well-documented
# organic spin constants (3ob's ``spinw.hsd`` and the DFTB+ manual's own H2O
# example, H = -0.072 / O Wpp = -0.028) are PBE values. The atomic spin constant
# is only weakly functional-dependent for H/C/N/O/S, so the authoritative 3ob
# values are the best available choice for mio too. Verified end-to-end: a mio OH
# radical runs spin-polarised and is stabilised relative to the restricted run.
SPIN_CONSTANTS_MIO = SPIN_CONSTANTS_3OB

# Default (3ob) spin constants.
SPIN_CONSTANTS = SPIN_CONSTANTS_3OB


def _spin_kwargs(atoms, spin, spin_constants=SPIN_CONSTANTS_3OB):
    """
    ASE ``Dftb`` keyword arguments enabling colinear spin polarisation.

    Returns an empty dict for ``spin`` in (None, 0), so the closed-shell
    (restricted) calculation is left exactly as before. For spin S > 0 it enables
    ``SpinPolarisation = Colinear`` with ``UnpairedElectrons = round(2*S)`` and
    injects the ``spin_constants`` for the elements present.

    Parameters
    ----------
    atoms : ase.Atoms
        The atoms whose elements need spin constants.
    spin : float or None
        Spin quantum number S.
    spin_constants : dict
        Element -> spin-constant brace string, matching the Slater-Koster set.

    Raises
    ------
    ValueError
        If an element has no tabulated spin constant.
    """
    if spin is None or float(spin) <= 0.0:
        return {}

    unpaired = int(round(2.0 * float(spin)))
    if unpaired <= 0:
        return {}

    elements = sorted(set(atoms.get_chemical_symbols()))
    missing = [element for element in elements if element not in spin_constants]
    if missing:
        raise ValueError(
            "Spin-polarised DFTB+ is not available for element(s) "
            f"{', '.join(missing)}: no spin constant is tabulated for this "
            "parameter set."
        )

    kwargs = {
        "Hamiltonian_SpinPolarisation": "Colinear {",
        "Hamiltonian_SpinPolarisation_UnpairedElectrons": unpaired,
        "Hamiltonian_SpinConstants_": "",
        "Hamiltonian_SpinConstants_ShellResolvedSpin": "No",
    }
    for element in elements:
        kwargs[f"Hamiltonian_SpinConstants_{element}"] = spin_constants[element]
    return kwargs


def _solvation_kwargs(solvent=None, param_file=None, install_root=None):
    """
    ASE ``Dftb`` keyword arguments enabling GBSA/ALPB implicit solvation.

    Returns an empty dict when neither ``solvent`` nor ``param_file`` is given,
    so the gas-phase calculation is left exactly as before. Otherwise it enables
    ``Solvation = GeneralizedBorn`` with the solvent's parameter file. DFTB+
    includes the solvation term in the energy, gradient, and Hessian, so the
    geometry, energy, and frequencies are all computed in solution consistently.

    The parameter files (grimme-lab/gbsa-parameters) are fit for GFN-xTB; used
    with the DFTB (3ob/mio) Hamiltonians they are an approximation. Pass an
    explicit ``param_file`` to use a method-consistent set instead.

    Parameters
    ----------
    solvent : str, optional
        Solvent name (e.g. ``"water"``). Resolved to a downloaded parameter file.
    param_file : str, optional
        Explicit path to a GBSA parameter file (overrides ``solvent``).
    install_root : str or Path, optional
        Root used to locate downloaded solvent parameters.

    Raises
    ------
    FileNotFoundError
        If the resolved parameter file does not exist.
    """
    if solvent is None and param_file is None:
        return {}

    if param_file is not None:
        path = os.path.abspath(os.path.expanduser(str(param_file)))
        if not os.path.isfile(path):
            raise FileNotFoundError(
                f"GBSA solvation parameter file does not exist: {path}"
            )
    else:
        # Lazy import: this is an environment/paths lookup, not a core dependency.
        from ..cli.dftb_setup import gbsa_param_path

        path = os.path.abspath(str(gbsa_param_path(solvent, install_root)))
        if not os.path.isfile(path):
            raise FileNotFoundError(
                f"No GBSA parameter file for solvent {solvent!r} at {path}. "
                f"Download it first, e.g. `thermo setup-dftb --solvent {solvent}`."
            )

    # DFTB+ resolves ParamFile relative to the run directory, so pass an absolute
    # path (jobs run in per-molecule working directories).
    return {
        "Hamiltonian_Solvation": "GeneralizedBorn {",
        "Hamiltonian_Solvation_ParamFile": path,
    }


# Grimme D3 with Becke-Johnson damping, using the parameters the DFTB+ manual
# recommends for the DFTB3/3ob Hamiltonian (Brandenburg et al., J. Chem. Phys.
# 143, 054110 (2015)).
_D3_BJ_3OB = {"a1": 0.5719, "a2": 3.6017, "s6": 1.0, "s8": 0.5883}


def _dispersion_kwargs(dispersion=None):
    """
    ASE ``Dftb`` keyword arguments adding a Grimme dispersion correction.

    Parameters
    ----------
    dispersion : str, optional
        Dispersion model. ``"d3-bj"`` adds Grimme D3 with Becke-Johnson damping
        using the 3ob-recommended parameters. ``None`` (default) adds nothing,
        so existing gas-phase results are unchanged.

    Returns
    -------
    dict
        The ``Hamiltonian_Dispersion*`` kwargs (empty when ``dispersion`` is None).

    Raises
    ------
    TSValueError
        If ``dispersion`` is not a supported model.
    """
    if dispersion is None:
        return {}

    if dispersion.lower() != "d3-bj":
        raise TSValueError(
            f"Unknown dispersion model {dispersion!r}; supported: 'd3-bj'."
        )

    p = _D3_BJ_3OB
    return {
        "Hamiltonian_Dispersion": "DftD3 {",
        "Hamiltonian_Dispersion_Damping": "BeckeJohnson {",
        "Hamiltonian_Dispersion_Damping_a1": p["a1"],
        "Hamiltonian_Dispersion_Damping_a2": p["a2"],
        "Hamiltonian_Dispersion_s6": p["s6"],
        "Hamiltonian_Dispersion_s8": p["s8"],
    }


class Geoopt(Dftb):
    """
    Custom DFTB+ calculator to optimize the system with the 'GeometryOptimisation' driver (Rational).
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
            for the DFTB_PREFIX environment variable.
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
            slako_dir=_slako_dir(slako_dir),
            Hamiltonian_Charge=charge,
            Driver_="GeometryOptimisation",
            Driver_Optimiser="Rational {}",
            Driver_MaxForceComponent=max_force,
            Driver_OutputPrefix="geo_opt",
            **kwargs,
        )

        self.calculate(atoms)


    def potential_energy(self):
        """
        Get the potential energy of the optimized geometry. in Hartree.

        Returns:
        --------
        float
            Potential energy of the optimized geometry.
        """
        self.atoms.calc = self
        return (self.atoms.get_potential_energy() * PhysicalConstants["eV"] /
                PhysicalConstants["H"])

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
        slako_dir=None,
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
            for the DFTB_PREFIX environment variable.

        Other Parameters:
        -----------------
        **kwargs : dict
            Additional keyword arguments to pass to the Dftb
        """

        super().__init__(
            atoms=atoms,
            label=label,
            slako_dir=_slako_dir(slako_dir),
            Hamiltonian_Charge=charge,
            Driver_="SecondDerivatives",
            Driver_Delta=delta,
            **kwargs,
        )

        self.calculate(atoms)


    def read(self):
        """
        Read the Hessian matrix from the 'hessian.out' file.

        Returns:
        --------
        numpy.ndarray
            Hessian matrix.
        """

        hessian_size = self.atoms.get_global_number_of_atoms() * 3
        self.hessian = _read_hessian_matrix("hessian.out", hessian_size)

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


    def write(self):
        """
        Writes the modes_in.hsd file.

        Returns:
        --------
        None
        """

        # use f-strings to write the modes_in.hsd file
        string = ("Geometry = GenFormat {\n"
                  f"    <<< {self.geometry}\n"
                  "}\n"
                  "\n"
                  "Hessian = {\n"
                  f"    <<< {self.hessian}\n"
                  "}\n"
                  "\n"
                  "Atoms = 1:-1\n"
                  "\n")

        with open("modes_in.hsd", "w") as f:
            f.write(string)


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

        if shutil.which("modes") is None:
            raise FileNotFoundError("The modes executable was not found.")

        with open("modes.out", "w") as output:
            subprocess.run(["modes"], stdout=output, check=True)


    def read(self):
        """
        Read the vibrational modes from the vibrations.tag file.

        Returns:
        --------
        numpy.ndarray
            Vibrational modes.
        """

        modes = []
        with open("vibrations.tag", "r", encoding="utf-8") as f:
            f.readline()  # skip the 'frequencies' tag header
            for line in f:
                try:
                    values = [float(field) for field in line.split()]
                except ValueError:
                    # stop at the next tag section (e.g. 'saved_modes :integer:..')
                    break
                modes.extend(values)

        # Hartree to cm^-1 - 1 Hartree = 219474.63 cm^-1
        self.wave_numbers = np.array(modes, dtype=float) * 219474.63

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

    # Parser options
    ParserOptions_ParserVersion=12,
)


# mio-1-1 (DFTB2): the original mio set is a second-order model, so it has no
# ThirdOrderFull / Hubbard derivatives. Same SCC + Fermi machinery as 3ob;
# spin-polarised runs use SPIN_CONSTANTS_MIO. Verified end-to-end against real
# DFTB+ on an OH radical.
dftb_mio_parameters = dict(
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

    # Analysis
    Analysis_="",
    Analysis_CalculateForces="Yes",
    Analysis_MullikenAnalysis="Yes",

    # Parser options
    ParserOptions_ParserVersion=12,
)


# Selectable DFTB parameter sets: name -> (Hamiltonian parameters, spin constants).
# ``dftbplus_thermo`` / ``screen`` pick a set by name so a run needs only a
# structure + charge (+ optional spin); everything else follows from the set.
DFTB_PARAMETER_SETS = {
    "3ob": (dftb_3ob_parameters, SPIN_CONSTANTS_3OB),
    "3ob-3-1": (dftb_3ob_parameters, SPIN_CONSTANTS_3OB),
    "mio": (dftb_mio_parameters, SPIN_CONSTANTS_MIO),
    "mio-1-1": (dftb_mio_parameters, SPIN_CONSTANTS_MIO),
}


def resolve_parameter_set(parameter_set):
    """
    Resolve a parameter-set name to ``(hamiltonian_parameters, spin_constants)``.

    Parameters
    ----------
    parameter_set : str
        One of the keys of :data:`DFTB_PARAMETER_SETS` (e.g. ``"3ob"``, ``"mio"``).

    Returns
    -------
    tuple(dict, dict)
        A fresh copy of the Hamiltonian parameters and the matching spin
        constants for the set.

    Raises
    ------
    ValueError
        If ``parameter_set`` is not a known set.
    """
    try:
        parameters, spin_constants = DFTB_PARAMETER_SETS[parameter_set]
    except KeyError:
        known = ", ".join(sorted(DFTB_PARAMETER_SETS))
        raise ValueError(
            f"Unknown DFTB parameter set {parameter_set!r}; choose one of: {known}."
        )
    return dict(parameters), spin_constants
