"""Molecular system model with symmetry and degree-of-freedom analysis."""

import logging
import warnings
import numpy as np

from beartype.typing import List
from numpy.exceptions import ComplexWarning
from PQAnalysis.atomic_system import AtomicSystem as PQAtomicSystem
from PQAnalysis.core.atom import Atom as PQAtom

from pymatgen.core import Molecule, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.analyzer import PointGroupAnalyzer

from ThermoScreening.utils.custom_logging import setup_logger
from ThermoScreening.exceptions import TSValueError
from ThermoScreening import __package_name__

from .atoms import Atom
from .cell import Cell


def _real_position(position: np.ndarray) -> np.ndarray:
    """
    Return a real floating-point atom position for symmetry analysis.
    """

    real_position = np.real_if_close(np.asarray(position))
    if np.iscomplexobj(real_position):
        if not np.allclose(np.imag(real_position), 0.0):
            raise TSValueError("Atom positions must be real-valued.")
        real_position = np.real(real_position)

    return np.asarray(real_position, dtype=float)


def _molecule_from_atoms(atoms: List[Atom]) -> Molecule:
    """
    Build a pymatgen molecule with real-valued coordinates.
    """

    names = []
    coordinates = []
    for atom in atoms:
        names.append(atom.symbol)
        coordinates.append(_real_position(atom.position))
    _validate_unique_positions(coordinates)

    return Molecule(names, coordinates)


def _validate_unique_positions(coordinates: List[np.ndarray]) -> None:
    """
    Reject overlapping atoms before symmetry analysis.
    """

    for index, position in enumerate(coordinates):
        for previous in coordinates[:index]:
            if np.allclose(position, previous, rtol=0.0, atol=1e-12):
                raise TSValueError(
                    "Atom positions must not overlap for symmetry analysis."
                )


def _point_group_analyzer(molecule):
    """
    Build a pymatgen point-group analyzer without surfacing benign cast warnings.
    """

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            category=ComplexWarning,
            module=r"pymatgen\.core\.operations",
        )
        return PointGroupAnalyzer(molecule)


def _pq_atomic_system(atoms: List[Atom]) -> PQAtomicSystem:
    return PQAtomicSystem(
        atoms=[PQAtom(atom.symbol) for atom in atoms],
        pos=np.asarray([_real_position(atom.position) for atom in atoms], dtype=float),
    )


def linearity(atoms: List[Atom]) -> bool:
    """
    Checks if the system is linear or non-linear.

    Parameters
    ----------
    atoms : List[Atom]
        A list of Atom objects.
        
    Return
    ------
    bool
        True if the system is linear, False otherwise.
        
    Raises
    ------
    TSValueError
        If the number of atoms is 1.
    """
    if len(atoms) == 1:
        raise TSValueError("Number of atoms must be greater than 1. The system is monoatomic.")
    if len(atoms) == 2:
        return True

    masses = np.array([atom.mass for atom in atoms], dtype=float)
    positions = np.array([_real_position(atom.position) for atom in atoms])
    # relocate to the center of mass so linearity is translation-invariant
    positions = positions - np.average(positions, axis=0, weights=masses)

    x, y, z = positions[:, 0], positions[:, 1], positions[:, 2]
    inertia_tensor = np.zeros((3, 3))
    inertia_tensor[0, 0] = np.sum(masses * (y**2 + z**2))
    inertia_tensor[1, 1] = np.sum(masses * (x**2 + z**2))
    inertia_tensor[2, 2] = np.sum(masses * (x**2 + y**2))
    inertia_tensor[0, 1] = inertia_tensor[1, 0] = -np.sum(masses * x * y)
    inertia_tensor[0, 2] = inertia_tensor[2, 0] = -np.sum(masses * x * z)
    inertia_tensor[1, 2] = inertia_tensor[2, 1] = -np.sum(masses * y * z)

    eigenvalues = np.linalg.eigvalsh(inertia_tensor)
    # a linear molecule has exactly one vanishing principal moment of inertia
    # (the molecular axis); use a relative tolerance rather than an exact zero.
    return bool(eigenvalues[0] <= 1e-6 * eigenvalues[-1])


def dimensionality(atoms: List[Atom]) -> int:
    """
    Check if the system is 1D, 2D or 3D for systems with more than 1 atom.

    Parameters
    ----------
    atoms : List[Atom]
        A list of Atom objects.

    return
    ------
    int
        The dimensionality of the system.
    """
    x = np.array([])
    y = np.array([])
    z = np.array([])
    for atom in atoms:
        pos = atom.position
        x = np.append(x, pos[0])
        y = np.append(y, pos[1])
        z = np.append(z, pos[2])
    zero_array = np.zeros(len(x))

    if (
        (np.array_equal(y, zero_array) and np.array_equal(z, zero_array))
        or (np.array_equal(x, zero_array) and np.array_equal(z, zero_array))
        or (np.array_equal(x, zero_array) and np.array_equal(y, zero_array))
    ):
        return 1
    if (
        np.array_equal(z, zero_array)
        or np.array_equal(y, zero_array)
        or np.array_equal(x, zero_array)
    ):
        return 2
    if (
        np.array_equal(z, zero_array)
        and np.array_equal(y, zero_array)
        and np.array_equal(x, zero_array)
    ):
        raise TSValueError("The system is 0D!")  # pragma: no cover
    return 3


def dim(atoms: List[Atom]) -> int:
    """
    Calculates the dimension of the system.

    Parameters
    ----------
    atoms : List[Atom]
        A list of Atom objects.

    Return
    ------
    int
        The dimension of the system.
        
    Raises
    ------
    TSValueError
        If the number of atoms is less than 1.
    """
    number_of_atoms = len(atoms)
    if number_of_atoms == 1:
        return 1
    if number_of_atoms > 1:
        return dimensionality(atoms)
    raise TSValueError("The number of atoms must be greater than 0.")


def dof(atoms: List[Atom]) -> int:
    """
    Calculates the degree of freedom of the system.
    Check if the system is monoatomic, linear or non-linear.
    
    Parameters
    ----------
    atoms : List[Atom]
        A list of Atom objects.
        
    Return
    ------
    int
        The degree of freedom of the system.
        
    Raises
    ------
    TSValueError
        If the number of atoms is less than 1.
    """
    number_of_atoms = len(atoms)
    if number_of_atoms == 1:
        # a monatomic species has no rotational or vibrational degrees of freedom
        return 0
    if number_of_atoms > 1:
        return (
            (3 * number_of_atoms - 5) if linearity(atoms) else (3 * number_of_atoms - 6)
        )
    raise TSValueError("The number of atoms must be greater than 0.")


def rotational_symmetry_number(atoms: List[Atom]) -> int:
    """
    Calculates the symmetry number of the system.

    Parameters
    ----------
    atoms : List[Atom]

    Returns
    -------
    int
        The symmetry number of the system.
    """
    if len(atoms) == 1:
        _real_position(atoms[0].position)
        return 1

    mol = _molecule_from_atoms(atoms)
    symmetry_number = _point_group_analyzer(mol).get_rotational_symmetry_number
    if callable(symmetry_number):
        return symmetry_number()
    return symmetry_number


def spacegroup_number(atoms: List[Atom], cell: Cell) -> int:
    """
    Calculates the spacegroup of the system.

    Parameters
    ----------
    atoms : List[Atom]

    Returns
    -------
    int
        The spacegroup of the system.
    """
    name = []
    coord = []
    for atom in atoms:
        name.append(atom.symbol)
        coord.append(atom.position)

    structure = Structure(cell.cell_vectors, name, coord)
    return SpacegroupAnalyzer(structure).get_space_group_number


def spacegroup(atoms: List[Atom], cell: Cell) -> str:
    """
    Calculates the spacegroup of the system.

    Parameters
    ----------
    atoms : List[Atom]

    Returns
    -------
    str
        The spacegroup of the system.
    """
    name = []
    coord = []
    for atom in atoms:
        name.append(atom.symbol)
        coord.append(atom.position)

    structure = Structure(cell.cell_vectors, name, coord)
    return SpacegroupAnalyzer(structure).get_space_group_symbol


def mass(atoms: List[Atom]) -> float:
    """
    Calculates the mass of the system.

    Parameters
    ----------
    atoms : List[Atom]

    Returns
    -------
    float
        The mass of the system in u.
    """
    return _pq_atomic_system(atoms).mass


def rotational_group_calc(atoms: List[Atom]) -> str:
    """
    Calculates the rotational group of the system.

    Parameters
    ----------
    atoms : List[Atom]

    Returns
    -------
    str
        The rotational group of the system.
    """
    if len(atoms) == 1:
        _real_position(atoms[0].position)
        return "Kh"

    mol = _molecule_from_atoms(atoms)
    symb = _point_group_analyzer(mol).sch_symbol
    return symb


def center_of_mass(atoms: List[Atom], _mass: float) -> np.ndarray:
    """
    Calculates the center of mass of the system.

    Parameters
    ----------
    atoms : List[Atom]

    Returns
    -------
    np.ndarray
        The center of mass of the system.
    """
    return _pq_atomic_system(atoms).center_of_mass


def imaginary_frequencies(vibrational_frequencies: np.ndarray) -> np.ndarray:
    """
    The imaginary vibrational frequencies of the system.

    Parameters
    ----------
    vibrational_frequencies : np.ndarray
        The vibrational frequencies of the system.

    Returns
    -------
    np.ndarray
        The imaginary vibrational frequencies of the system.
    """
    return vibrational_frequencies[vibrational_frequencies < 0]


def check_imaginary_frequencies(imaginary_frequencies: np.ndarray) -> bool:
    """
    Checks if the system has imaginary vibrational frequencies.

    Parameters
    ----------
    imaginary_frequencies : np.ndarray

    Returns
    -------
    bool
        True if the system has imaginary vibrational frequencies, False otherwise.
    """
    return len(imaginary_frequencies) > 0


def cleaned_frequency(frequency: np.ndarray):
    """
    Removes the imaginary part of the frequency.

    Parameters
    ----------
    frequency : np.ndarray

    Returns
    -------
    np.ndarray
        The cleaned frequency.
    """
    return np.delete(frequency, np.where(frequency < 0))


def frequency_dof(frequency: np.ndarray, dof: int) -> np.ndarray:
    """
    Delete all frequencies that are more than the degree of freedom.

    Parameters
    ----------
    frequency : np.ndarray
        The vibrational frequencies of the system.
    dof : int
        The degree of freedom of the system.

    Returns
    -------
    frequency : np.ndarray
        The vibrational frequencies of the system.
    """
    N = len(frequency)
    cleaned_freq = np.empty(dof)
    i = 0

    while (i + 1) <= dof:
        cleaned_freq[i] = frequency[N - 1 - i]
        i += 1

    cleaned_freq = cleaned_freq[::-1]
    return cleaned_freq


def check_frequency_length(frequency: np.ndarray, dof: int) -> bool:
    """
    Checks if the system has the correct number of vibrational frequencies.

    Parameters
    ----------
    frequency : np.ndarray
        The vibrational frequencies of the system.
    dof : int
        The degree of freedom of the system.

    Returns
    -------
    bool
        True if the system has the correct number of vibrational frequencies, False otherwise.
    """
    return len(frequency) == dof


class System:
    """
    A class for storing the system information of a chemical system.

    Attributes
    ----------
    atoms : List[Atom]
        A list of Atom objects.
    charge : float
        The charge of the system.
    periodicity : bool
        Whether the system has periodic boundary conditions.
    pbc : List[bool]
        The periodic boundary conditions of the system.
    cell : Cell
        The cell of the system.
    solvation : bool
        Whether the system is solvated.
    solvent : str
        The solvent of the system.
    electronic_energy : float
        The electronic energy of the system.
    vibrational_frequencies : np.ndarray
        The vibrational frequencies of the system.

    Methods
    -------
    atom_names()
        The atom names of the system.
    atomic_masses()
        A list of the atomic masses of the system.
    coord()
        The atom positions of the system.
    x()
        The x coordinates of the system.
    y()
        The y coordinates of the system.
    z()
        The z coordinates of the system.

    Properties
    ----------
    number_of_atoms : int
        The number of atoms in the system.
    dim : int
        The dimension of the system.
    charge : float
        The charge of the system.
    dof : int
        The degree of freedom of the system.
    spin : float
        The spin of the system.
    rotational_symmetry_number : int
        The symmetry number of the system.
    rotational_group : str
        The rotational group of the system.
    spacegroup_number : int
        The spacegroup of the system.
    spacegroup : str
        The spacegroup of the system.
    center_of_mass : np.ndarray
        The center of mass of the system.
    periodicity : bool
        The periodicity of the system.
    pbc : List[bool]
        The periodic boundary conditions of the system.
    cell : np.ndarray
        The cell of the system.
    solvation : str
        The solvation of the system.
    solvent : str
        The solvent of the system.
    electronic_energy : float
        The electronic energy of the system.
    vibrational_frequencies : np.ndarray
        The vibrational frequencies of the system.
    imaginary_frequencies : np.ndarray
        The imaginary vibrational frequencies of the system.
    has_imaginary_frequencies : bool
        Whether the system has imaginary vibrational frequencies.
    mass : float
        The mass of the system.
    real_vibrational_frequencies : np.ndarray
        The real frequency of the system.
    check_frequency_length : bool
        Whether the system has the correct number of vibrational frequencies.

    Examples
    --------
    >>> from ThermoScreening.thermo.system import System
    >>> from ThermoScreening.thermo.atoms import Atom
    >>> from ThermoScreening.thermo.cell import Cell
    >>> import numpy as np
    >>> atom1 = Atom("H", [0.0, 0.0, 0.0], 1.0)
    >>> atom2 = Atom("H", [0.0, 0.0, 1.0], 1.0)
    >>> atoms = [atom1, atom2]
    >>> cell = Cell([10.0, 10.0, 10.0])
    >>> system = System(atoms=atoms, cell=cell, vibrational_frequencies=np.array([1.0, 2.0, 3.0]))
    >>> system.atoms
    [Atom(symbol='H', position=[0.0, 0.0, 0.0], mass=1.0), Atom(symbol='H', position=[0.0, 0.0, 1.0], mass=1.0)]
    """

    logger = logging.getLogger(__package_name__).getChild(__name__)
    logger = setup_logger(logger)

    def __init__(
        self,
        atoms: List[Atom] | None = None,
        charge: float | None = None,
        periodicity: bool | None = None,
        pbc: List[bool] | bool | None = None,
        cell: Cell | np.ndarray | None = None,
        solvation: bool | None = None,
        solvent: str | None = None,
        electronic_energy: float | None = None,
        vibrational_frequencies: np.ndarray | None = None,
        spin: float | None = None,
    ) -> None:
        """
        Initializes the System with the given parameters.

        Parameters
        ----------
        atoms : List[Atom], optional, default=[]
            A list of Atom objects.
        charge : Float, optional, default=0
            The charge of the system.
        periodicity : bool, optional, default=False
            Whether the system has periodic boundary conditions.
        pbc : List[bool], optional, default=[False, False, False]
            The periodic boundary conditions of the system.
        cell : Cell, optional, default=Cell()
            The periodic boundary conditions of the system.
        solvation : bool, optional, default=False
            Whether the system is solvated.
        solvent : str, optional, default=''
            The solvent of the system.
        electronic_energy : float, optional, default=0.0
            The electronic energy of the system.
        vibrational_frequencies : np.ndarray, optional, default=None
            The vibrational frequencies of the system.

        Raises
        ------
        TSValueError
            If the number of atoms does not match the length of the atoms list.
            If the number of atoms is less than 1.
            If the vibrational frequencies are not provided.
        """

        if atoms is None:
            raise TSValueError("Atoms must be provided.")

        if len(atoms) == 0:
            raise TSValueError("The number of atoms must be greater than 0.")

        if vibrational_frequencies is None:
            raise TSValueError("Vibrational frequencies must be provided.")

        self._atoms = atoms
        self._charge = charge
        self._periodicity = periodicity
        self._pbc = pbc
        self._cell = cell
        if self._cell is None:
            self._cell = None
        self._solvation = solvation
        self._solvent = solvent
        self._electronic_energy = electronic_energy
        self._vibrational_frequencies = vibrational_frequencies
        self._dof = dof(atoms)
        self._dim = dim(atoms)
        self._mass = mass(atoms)
        self._center_of_mass = center_of_mass(atoms, self._mass)
        self._symmetry_number = rotational_symmetry_number(atoms)
        self._rotational_group = rotational_group_calc(atoms)
        if isinstance(self._cell, Cell):
            self._spacegroup_number = spacegroup_number(atoms, self._cell)
            self._spacegroup = spacegroup(atoms, self._cell)
        else:
            # Periodic input arrives as a raw cell array (from read_xyz/read_gen)
            # and the spacegroup is unused by the gas-phase RRHO calculation, so
            # degrade to None instead of feeding a non-Cell into spacegroup().
            self._spacegroup_number = None
            self._spacegroup = None
        self._imaginary_frequencies = imaginary_frequencies(vibrational_frequencies)
        self._has_imaginary_frequencies = check_imaginary_frequencies(
            self._imaginary_frequencies
        )
        self._real_vibrational_frequencies = frequency_dof(
            self._vibrational_frequencies, self._dof
        )
        self._check_frequency_length = check_frequency_length(
            self._real_vibrational_frequencies, self._dof
        )

        if charge is None:
            self._charge = 0.0

        if periodicity is None:
            self._periodicity = False

        if pbc is None:
            self._pbc = [False, False, False]

        if solvation is None:
            self._solvation = False

        if solvent is None:
            self._solvent = ""

        if electronic_energy is None:
            self._electronic_energy = 0.0

        # Default to a closed-shell singlet (consistent with the spin-restricted
        # DFTB+ calculation); open-shell species (radicals, triplet O2, ...) must
        # set the spin explicitly.
        if spin is None:
            self._spin = 0.0
        else:
            if spin < 0:
                raise TSValueError("The spin must be non-negative.")
            self._spin = float(spin)

    @property
    def atoms(self) -> List[Atom]:
        """
        The atoms of the system.

        Returns
        -------
        List[Atom]
            The atoms of the system.
        """
        return self._atoms

    def atom_names(self) -> List[str]:
        """
        The atom names of the system.

        Returns
        -------
        List[str]
            The atom names of the system.
        """
        return [atom.symbol for atom in self._atoms]

    def atomic_masses(self) -> List[float]:
        """
        A list of the atomic masses of the system.

        Returns
        -------
        List[float]
            A list of the atomic masses of the system.
        """
        return [atom.mass for atom in self._atoms]

    def coord(self) -> np.ndarray:
        """
        The atom positions of the system.

        Returns
        -------
        np.ndarray
            The atom positions of the system.
        """
        return np.array([atom.position for atom in self._atoms])

    @property
    def number_of_atoms(self) -> int:
        """
        The number of atoms in the system.

        Returns
        -------
        int
            The number of atoms in the system.
        """
        return len(self.atoms)

    @property
    def dim(self) -> int:
        """
        The dimension of the system.

        Returns
        -------
        int
            The dimension of the system.
        """
        return self._dim

    @property
    def charge(self) -> float:
        """
        The charge of the system.

        Returns
        -------
        float
            The charge of the system.
        """
        return self._charge

    @property
    def dof(self) -> int:
        """
        The degree of freedom of the system.

        Returns
        -------
        int
            The degree of freedom of the system.
        """
        return self._dof

    @property
    def spin(self) -> float:
        """
        The spin of the system.

        Returns
        -------
        float
            The spin of the system.
        """
        return self._spin

    @property
    def rotational_symmetry_number(self) -> int:
        """
        The symmetry number of the system.

        Returns
        -------
        int
            The symmetry number of the system.
        """
        return self._symmetry_number

    @property
    def rotational_group(self) -> str:
        """
        The rotational group of the system.

        Returns
        -------
        str
            The rotational group of the system.
        """

        return self._rotational_group

    @property
    def spacegroup_number(self) -> None | int:
        """
        The spacegroup number of the system.

        Returns
        -------
        int, None
            The spacegroup number of the system.
        """

        return self._spacegroup_number

    @property
    def spacegroup(self) -> None | str:
        """
        The spacegroup of the system.

        Returns
        -------
        str,None
            The spacegroup of the system.
        """

        return self._spacegroup

    @property
    def center_of_mass(self) -> np.ndarray:
        """
        The center of mass of the system.

        Returns
        -------
        np.ndarray
            The center of mass of the system.
        """
        return self._center_of_mass

    @property
    def periodicity(self) -> bool:
        """
        The periodicity of the system.

        Returns
        -------
        bool
            The periodicity of the system.
        """
        return self._periodicity

    @property
    def pbc(self) -> List[bool] | bool:
        """
        The periodic boundary conditions of the system.

        Returns
        -------
        List[bool]
            The periodic boundary conditions of the system e.g. [True, True, True] for a 3D system, 
            [True, False, False] for a 1D system or [True, True, False] for a 2D system.
        """
        return self._pbc

    @property
    def cell(self) -> Cell | np.ndarray | None:
        """
        The cell of the system.

        Returns
        -------
        np.ndarray
            The cell of the system.
        """
        return self._cell

    @property
    def solvation(self) -> bool:
        """
        The solvation of the system.

        Returns
        -------
        str
            The solvation of the system.
        """
        return self._solvation

    @property
    def solvent(self) -> str:
        """
        The solvent of the system.

        Returns
        -------
        str
            The solvent of the system.
        """
        return self._solvent

    @property
    def electronic_energy(self) -> float:
        """
        The electronic energy of the system.

        Returns
        -------
        float
            The electronic energy of the system.
        """
        return self._electronic_energy

    @property
    def vibrational_frequencies(self) -> np.ndarray:
        """
        The vibrational frequencies of the system.

        Returns
        -------
        np.ndarray
            The vibrational frequencies of the system.
        """
        return self._vibrational_frequencies

    @property
    def imaginary_frequencies(self) -> np.ndarray:
        """
        The imaginary vibrational frequencies of the system.

        Returns
        -------
        np.ndarray
            The imaginary vibrational frequencies of the system.
        """
        return self._imaginary_frequencies

    @property
    def has_imaginary_frequencies(self) -> bool:
        """
        Whether the system has imaginary vibrational frequencies.

        Returns
        -------
        bool
            True if the system has imaginary vibrational frequencies, False otherwise.
        """
        return self._has_imaginary_frequencies

    @property
    def mass(self) -> float:
        """
        The mass of the system.

        Returns
        -------
        float
            The mass of the system.
        """
        return self._mass

    @property
    def real_vibrational_frequencies(self) -> np.ndarray:
        """
        The real frequency of the system.

        Returns
        -------
        np.ndarray
            The real frequency of the system.
        """
        return self._real_vibrational_frequencies

    @property
    def check_frequency_length(self) -> bool:
        """
        Whether the system has the correct number of vibrational frequencies.

        Returns
        -------
        bool
            True if the system has the correct number of vibrational frequencies, 
            False otherwise.
        """
        return self._check_frequency_length

    def x(self) -> np.ndarray:
        """
        The x coordinates of the system.

        Returns
        -------
        np.ndarray
            The x coordinates of the system.
        """
        return np.array([atom.position[0] for atom in self._atoms])

    def y(self) -> np.ndarray:
        """
        The y coordinates of the system.

        Returns
        -------
        np.ndarray
            The y coordinates of the system.
        """
        return np.array([atom.position[1] for atom in self._atoms])

    def z(self) -> np.ndarray:
        """
        The z coordinates of the system.

        Returns
        -------
        np.ndarray
            The z coordinates of the system.
        """
        return np.array([atom.position[2] for atom in self._atoms])
