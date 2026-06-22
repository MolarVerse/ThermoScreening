import logging

import numpy as np
from PQAnalysis.io import read_gen_file

from ThermoScreening.exceptions import TSNotImplementedError, TSValueError
from ThermoScreening.utils.custom_logging import setup_logger
from ThermoScreening import __package_name__

from .inputFileReader import InputFileReader
from .system import System, dof
from .thermo import Thermo
from .atoms import Atom
from ..calculator import Geoopt, Hessian, Modes


logger = logging.getLogger(__package_name__).getChild("api")
logger = setup_logger(logger)


def read_xyz(coord_file: str):
    """
    Reads xyz-file and returns number of atoms, chemical symbols of atoms, coordinates of atoms, and cell parameters.

    Parameter:
    ----------
    coord_file : str
        The name of the coordinate file.

    Returns:
    --------
    List
        number of atoms, chemical symbols of atoms, coordinates of atoms, cell parameters
    """

    cell = None
    pbc = False
    data_N = None

    with open(coord_file, "r") as f:

        line = f.readline()

        line = line.strip()

        # split line and ignore spaces
        line = line.split()
        data_N = int(line[0])

        if len(line) >= 7:
            try:
                cell = np.array([float(value) for value in line[1:7]])
            except ValueError:
                cell = None
            else:
                pbc = True

        line = f.readline().strip()

        i = 0
        data_atoms = np.empty(data_N, dtype=object)
        data_xyz = np.zeros((data_N, 3), dtype=float)
        while True:
            line = f.readline().strip()
            line = line.split()
            if len(line) == 0:
                break
            data_atoms[i] = str(line[0])
            data_xyz[i, 0] = float(line[1])
            data_xyz[i, 1] = float(line[2])
            data_xyz[i, 2] = float(line[3])
            i += 1

    return [data_N, data_atoms, data_xyz, cell, pbc]


def read_gen(coord_file: str):
    """
    Reads gen-file and returns number of atoms, chemical symbols of atoms, coordinates of atoms, and cell vector.

    Parameter:
    ----------
    coord_file : str
        The name of the coordinate file.

    Returns:
    --------
    List
        number of atoms, chemical symbols of atoms, coordinates of atoms, cell vector
    """
    system = read_gen_file(coord_file)

    data_N = system.n_atoms
    data_atoms = np.array([atom.name for atom in system.atoms], dtype=object)
    data_xyz = np.asarray(system.pos, dtype=float)

    if system.cell.is_vacuum:
        cell_vector = None
        pbc = False
    else:
        cell_vector = np.asarray(system.cell.box_matrix, dtype=float)
        cell_vector[np.isclose(cell_vector, 0.0, atol=1e-12)] = 0.0
        pbc = True

    return [data_N, data_atoms, data_xyz, cell_vector, pbc]


def read_vib_file(vibrational_file: str):
    """
    Reads the vibrational file and returns the vibrational frequencies as a numpy array.

    Parameter:
    ----------
    vibrational_file : str
        The name of the vibrational file.

    Returns:
    --------

    np.ndarray
        The vibrational frequencies as a numpy array.
    """
    vibrational_frequencies = np.array([])
    with open(vibrational_file, "r") as f:
        while True:
            line = f.readline().strip()
            line = line.split()
            if len(line) == 0:
                break
            vibrational_frequencies = np.append(vibrational_frequencies, float(line[1]))

    return vibrational_frequencies


def read_coord(coord_file: str, engine: str):
    """
    Reads the coordinate file and returns the coordinates as a numpy array.

    Parameters
    ----------
    coord_file : str
        The name of the coordinate file.
    engine : str
        The name of the engine.

    Returns
    -------
    List
        number of atoms, chemical symbols of atoms, coordinates of atoms
        
    Raises
    ------
    TSNotImplementedError
        If the engine is not supported.
        If the input file is not supported.
        If the gen file is not tested yet.
    """

    data_N, data_atoms, data_xyz, cell = None, None, None, None
    if engine == "dftb+":
        if coord_file.endswith(".xyz"):
            data_N, data_atoms, data_xyz, cell, pbc = read_xyz(coord_file)
            return data_N, data_atoms, data_xyz, cell, pbc
        elif coord_file.endswith(".gen"):
            data_N, data_atoms, data_xyz, cell_vectors, pbc = read_gen(coord_file)
            return data_N, data_atoms, data_xyz, cell_vectors, pbc
        else:
            logger.error(
                "The input file is not supported.",
                exception=TSNotImplementedError
            )

    else:
        logger.error(
            "The engine is not supported.",
            exception=TSNotImplementedError
        )


def read_vibrational(vibrational_file: str, engine: str):
    """
    Reads the vibrational file and returns the vibrational frequencies as a numpy array.

    Parameters
    ----------
    vibrational_file : str
        The name of the vibrational file.
    engine : str

    Returns
    -------
    np.ndarray
        The vibrational frequencies as a numpy array.
        
    Raises
    ------
    TSNotImplementedError
        If the engine is not supported.
    """
    vibrational_frequencies = None
    if engine == "dftb+":
        vibrational_frequencies = read_vib_file(vibrational_file)
    else:
        # throw error
        logger.error(
            "The engine is not supported.",
            exception=TSNotImplementedError
        )

    return vibrational_frequencies


def unit_length(engine: str):
    """
    Returns the unit of length.

    Parameter:
    ----------
    engine : str
        The name of the engine.

    Returns:
    --------
    str
        The unit of length.
        
    Raises:
    -------
    TSNotImplementedError
        If the engine is not supported.
    """
    if engine == "dftb+":
        return "Angstrom"
    else:
        logger.error(
            "The engine is not supported.",
            exception=TSNotImplementedError
        )


def unit_energy(engine: str):
    """
    Returns the unit of energy.

    Parameter:
    ----------
    engine : str
        The name of the engine.

    Returns:
    --------
    str
        The unit of energy.
        
    Raises:
    -------
    TSNotImplementedError
        If the engine is not supported.
    """
    if engine == "dftb+":
        return "Hartree"
    else:
        logger.error(
            "The engine is not supported.",
            exception=TSNotImplementedError
        )


def unit_mass(engine: str):
    """
    Returns the unit of mass.

    Parameter:
    ----------
    engine : str
        The name of the engine.

    Returns:
    --------
    str
        The unit of mass.
        
    Raises:
    -------
    TSNotImplementedError
        If the engine is not supported.
    """
    if engine == "dftb+":
        return "amu"
    else:
        logger.error(
            "The engine is not supported.",
            exception=TSNotImplementedError
        )


def unit_frequency(engine: str):
    """
    Returns the unit of frequency.

    Parameter:
    ----------
    engine : str
        The name of the engine.

    Returns:
    --------
    str
        The unit of frequency.
        
    Raises:
    -------
    TSNotImplementedError
        If the engine is not supported.
    """
    if engine == "dftb+":
        return "cm^-1"
    else:
        logger.error(
            "The engine is not supported.",
            exception=TSNotImplementedError
        )


def run_thermo(
    vibrational_frequencies,
    coord_file="geo_opt.xyz",
    temperature=298.15,
    pressure=101325,
    energy=0.0,
    engine="dftb+",
    charge=0.0,
):
    """
    Run the thermo calculation. Returns thermo calculation object.

    Parameters
    ----------
    vibrational_frequencies : np.ndarray
        Vibrational frequencies in cm^-1.
    coord_file : str
        The name of the coordinate file. Default is 'geo_opt.xyz'.
    temperature : float
        The temperature in K. Default is 298.15.
    pressure : float
        The pressure in Pa. Default is 101325.
    energy : float
        The energy in Hartree. Default is 0.0.
    engine : str
        The name of the engine. Default is 'dftb+'.
    charge : float
        The system charge.

    The coordinate file should be in xyz format.


    Returns
    -------
    Thermo
        The thermo calculation object.
        
    Raises
    ------
    TSValueError
        If the number of vibrational frequencies does not match with the degree of freedom.
    """
    # Read coordinate file
    _, atom_symbol, xyz, cell, pbc = read_coord(coord_file, engine)

    # Initialize atoms
    atom_list = []
    for i, symbol in enumerate(atom_symbol):
        atom_list.append(Atom(symbol=symbol, position=xyz[i, :]))

    expected_dof = dof(atom_list)
    if len(vibrational_frequencies) < expected_dof:
        logger.error(
            "The number of vibrational frequencies does not match with the degree of freedom.",
            exception=TSValueError
        )

    system_info = System(
        atoms=atom_list,
        electronic_energy=energy,
        cell=cell,
        vibrational_frequencies=vibrational_frequencies,
        pbc=pbc,
        charge=charge,
    )

    thermo_setup = Thermo(
        system=system_info, temperature=temperature, pressure=pressure, engine=engine
    )

    thermo_setup.run()

    return thermo_setup


def execute(input_file: str) -> Thermo:
    """
    Execute the thermo calculation. Returns thermo calculation object.
    Input file should be in the following format:

    ```
    coord_file=geo_opt.xyz
    vibrational_file=vibrational.out
    temperature=298.15
    pressure=101325
    energy=0.0
    engine=dftb+
    ```

    Parameters
    ----------
    input_file : str
        The name of the input file.

    Returns
    -------
    Thermo
        The thermo calculation object.
    """
    # Read input file
    input_file_reader = InputFileReader(input_file)
    input_dictionary = input_file_reader._dictionary

    # Read coordinate file
    coordinate_file = input_dictionary["coord_file"]
    engine = input_dictionary["engine"]

    # Read vibrational file
    vibrational_file = input_dictionary["vibrational_file"]
    vibrational_frequencies = read_vibrational(vibrational_file, engine)

    # Read temperature
    temperature = float(input_dictionary["temperature"])

    # Read pressure
    pressure = float(input_dictionary["pressure"])

    # Read energy
    energy = float(input_dictionary["energy"])

    return run_thermo(
        vibrational_frequencies,
        coord_file=coordinate_file,
        temperature=temperature,
        pressure=pressure,
        energy=energy,
        engine=engine,
    )

def dftbplus_thermo(
        atoms, 
        temperature=298.15,
        pressure=101325,
        charge=0.0,
        **kwargs
    ):
    """
    Run the thermo calculation with geometry optimization. Returns thermo calculation object.

    Parameters
    ----------
    atoms : ase.Atoms
        Atoms object. The initial geometry.
    temperature : float
        The temperature in K. Default is 298.15.
    pressure : float
        The pressure in Pa. Default is 101325.
    charge : float
        The system charge. Default is 0.0.

    Other Parameters
    ----------------
    **kwargs : dict
        Additional keyword arguments to pass to the Geoopt and Hessian classes.

    Returns
    -------
    Thermo
        The thermo calculation object.
    """
    
    # run geometry optimization
    geoopt = Geoopt(atoms=atoms, charge=charge, **kwargs)
    potential_energy = geoopt.potential_energy()
    optimized_atoms = geoopt.read()
    
    # run hessian calculation
    Hessian(atoms=optimized_atoms, charge=charge, **kwargs)

    # run normal mode calculation
    modes = Modes()
    frequencies = modes.wave_numbers

    # run thermo calculation
    thermo = run_thermo(
        frequencies,
        coord_file="geo_opt.xyz",
        temperature=temperature,
        pressure=pressure,
        energy=potential_energy,
        engine='dftb+',
        charge=charge,
    )

    return thermo
