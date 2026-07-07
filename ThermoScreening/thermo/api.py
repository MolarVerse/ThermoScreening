"""High-level thermochemistry functions: coordinate/frequency I/O and run helpers."""

import logging
import os
from contextlib import contextmanager
from pathlib import Path

import numpy as np
from PQAnalysis.io import XYZFrameReader, read_gen_file
from PQAnalysis.io.traj_file.exceptions import FrameReaderError

from ThermoScreening.exceptions import TSNotImplementedError, TSValueError
from ThermoScreening.utils.custom_logging import setup_logger
from ThermoScreening import __package_name__

from .inputFileReader import InputFileReader
from .system import System, dof
from .thermo import Thermo
from .atoms import Atom
from ..calculator import Geoopt, Hessian, Modes
from ..calculator.dftbplus import _spin_kwargs, _solvation_kwargs, _dispersion_kwargs, SPIN_CONSTANTS_3OB
from ..calculator.xtb import optimise_and_frequencies, xtb_calculator
from ..calculator.xtb_cli import run_xtb


logger = logging.getLogger(__package_name__).getChild("api")
logger = setup_logger(logger)


def _pq_atom_names(system):
    return np.array([atom.name for atom in system.atoms], dtype=object)


def _pq_positions(system):
    return np.asarray(system.pos, dtype=float)


def _xyz_float64_positions(text: str, n_atoms: int):
    """
    Parse the atom-block coordinates of an XYZ string as float64.

    ``XYZFrameReader`` stores frame positions as single precision (float32);
    re-reading the coordinate columns here preserves the full double precision
    of the source file, matching :func:`read_gen`. The frame has already
    validated that the file is well-formed, so the atom lines follow the count
    and comment lines.
    """
    atom_lines = text.splitlines()[2 : 2 + n_atoms]
    return np.array(
        [line.split()[1:4] for line in atom_lines],
        dtype=float,
    )


def _pq_xyz_cell(cell):
    if cell.is_vacuum:
        return None, False
    return (
        np.array([cell.x, cell.y, cell.z, cell.alpha, cell.beta, cell.gamma]),
        True,
    )


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

    text = Path(coord_file).read_text(encoding="utf-8")

    try:
        frame = XYZFrameReader().read(text, traj_format="xyz")
    except (FrameReaderError, ValueError, IndexError) as exc:
        raise TSValueError("Invalid XYZ coordinate file.") from exc

    cell, pbc = _pq_xyz_cell(frame.cell)
    positions = _xyz_float64_positions(text, frame.n_atoms)
    return [frame.n_atoms, _pq_atom_names(frame), positions, cell, pbc]


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
    data_atoms = _pq_atom_names(system)
    data_xyz = _pq_positions(system)

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
    vibrational_frequencies = []
    with open(vibrational_file, "r", encoding="utf-8") as f:
        for line_number, line in enumerate(f, start=1):
            fields = line.split()
            if not fields:
                continue
            if len(fields) < 2:
                raise TSValueError(
                    f"Invalid vibrational frequency line {line_number}."
                )
            try:
                vibrational_frequencies.append(float(fields[1]))
            except ValueError as exc:
                raise TSValueError(
                    f"Invalid vibrational frequency line {line_number}."
                ) from exc

    if not vibrational_frequencies:
        raise TSValueError("No vibrational frequencies found.")

    return np.asarray(vibrational_frequencies)


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
        if coord_file.endswith(".gen"):
            data_N, data_atoms, data_xyz, cell_vectors, pbc = read_gen(coord_file)
            return data_N, data_atoms, data_xyz, cell_vectors, pbc
        raise TSNotImplementedError("The input file is not supported.")

    raise TSNotImplementedError("The engine is not supported.")


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
        raise TSNotImplementedError("The engine is not supported.")

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
    raise TSNotImplementedError("The engine is not supported.")


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
    raise TSNotImplementedError("The engine is not supported.")


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
    raise TSNotImplementedError("The engine is not supported.")


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
    raise TSNotImplementedError("The engine is not supported.")


def _atoms_from_ase(atoms):
    """
    Build the Atom list, cell and pbc from an ASE Atoms object.

    Mirrors the (symbols, positions, cell, pbc) tuple that read_coord produces,
    so an optimized geometry can be passed straight through without a disk
    round-trip.
    """
    positions = np.asarray(atoms.get_positions(), dtype=float)
    atom_list = [
        Atom(symbol=symbol, position=positions[i])
        for i, symbol in enumerate(atoms.get_chemical_symbols())
    ]
    pbc = bool(np.any(atoms.get_pbc()))
    cell = np.asarray(atoms.get_cell(), dtype=float) if pbc else None
    return atom_list, cell, pbc


def run_thermo(
    vibrational_frequencies,
    coord_file="geo_opt.xyz",
    temperature=298.15,
    pressure=101325,
    energy=0.0,
    engine="dftb+",
    charge=0.0,
    atoms=None,
    spin=None,
    quasi_rrho=False,
):
    """
    Run the thermo calculation. Returns thermo calculation object.

    Parameters
    ----------
    vibrational_frequencies : np.ndarray
        Vibrational frequencies in cm^-1.
    coord_file : str
        The name of the coordinate file. Default is 'geo_opt.xyz'. Ignored when
        ``atoms`` is provided.
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
    atoms : ase.Atoms, optional
        An optimized geometry to use directly instead of reading ``coord_file``.
    quasi_rrho : bool
        If True, use Grimme's quasi-RRHO treatment for the vibrational entropy
        instead of the pure harmonic oscillator. Default False.

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
    # Build the atom list either from an optimized ASE Atoms object or by
    # reading the coordinate file.
    if atoms is not None:
        atom_list, cell, pbc = _atoms_from_ase(atoms)
    else:
        _, atom_symbol, xyz, cell, pbc = read_coord(coord_file, engine)
        atom_list = [
            Atom(symbol=symbol, position=xyz[i, :])
            for i, symbol in enumerate(atom_symbol)
        ]

    expected_dof = dof(atom_list)
    if len(vibrational_frequencies) < expected_dof:
        raise TSValueError(
            "The number of vibrational frequencies does not match with the degree of freedom."
        )

    system_info = System(
        atoms=atom_list,
        electronic_energy=energy,
        cell=cell,
        vibrational_frequencies=vibrational_frequencies,
        pbc=pbc,
        charge=charge,
        spin=spin,
    )

    thermo_setup = Thermo(
        system=system_info,
        temperature=temperature,
        pressure=pressure,
        engine=engine,
        quasi_rrho=quasi_rrho,
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

@contextmanager
def _run_in_directory(directory):
    """
    Run the enclosed block in ``directory`` (created if needed), restoring the
    previous working directory afterwards. A no-op when ``directory`` is None.

    The DFTB+ steps read and write fixed filenames (geo_opt.gen, hessian.out,
    vibrations.tag, ...) in the current directory, so giving each job its own
    directory keeps batch runs from clobbering each other.
    """
    if directory is None:
        yield
        return

    previous = Path.cwd()
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    os.chdir(directory)
    try:
        yield
    finally:
        os.chdir(previous)


def dftbplus_thermo(
        atoms,
        temperature=298.15,
        pressure=101325,
        charge=0.0,
        directory=None,
        spin=None,
        spin_constants=None,
        solvent=None,
        solvation_param_file=None,
        dispersion=None,
        quasi_rrho=False,
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
    directory : str, optional
        Working directory to run the DFTB+ steps in (created if needed). The
        geometry/Hessian/modes files are written here, so separate jobs can run
        without clobbering each other. Defaults to the current directory.

    spin : float, optional
        Spin quantum number S. Defaults to the minimum-spin electron-count guess
        (even -> 0, odd -> 0.5). When S > 0 the DFTB+ steps run colinear
        spin-polarised (so radicals are treated open-shell automatically); S = 0
        keeps the restricted closed-shell calculation.
    spin_constants : dict, optional
        Element -> spin-constant mapping matching the Slater-Koster set in use.
        Defaults to the 3ob constants; pass ``SPIN_CONSTANTS_MIO`` (or the value
        from :func:`resolve_parameter_set`) when running the mio set.
    solvent : str, optional
        Solvent name (e.g. ``"water"``) for GBSA/ALPB implicit solvation. When
        set, the optimisation, energy, and Hessian all run in solution. Requires
        the solvent's parameter file (``thermo setup-dftb --solvent <name>``).
        Defaults to a gas-phase calculation.
    solvation_param_file : str, optional
        Explicit path to a GBSA parameter file, overriding ``solvent`` (use a
        method-consistent set instead of the default GFN-fit one).
    dispersion : str, optional
        Dispersion correction to add to the Hamiltonian. ``"d3-bj"`` uses Grimme
        D3 with Becke-Johnson damping (3ob-recommended parameters). Defaults to
        no dispersion. Requires a DFTB+ build with D3 support.
    quasi_rrho : bool
        If True, use Grimme's quasi-RRHO treatment for the vibrational entropy,
        which tames the entropy of low-frequency modes. Default False.

    Other Parameters
    ----------------
    **kwargs : dict
        Additional keyword arguments to pass to the Geoopt and Hessian classes.

    Returns
    -------
    Thermo
        The thermo calculation object.

    Notes
    -----
    Reproducibility: identical inputs on the same machine give bit-identical
    results (DFTB+ is thread-invariant). Across machines the optimiser stopping
    point can vary by ~3e-5 A, giving ~2e-9 Hartree / ~3e-4 cal/mol/K scatter --
    physically negligible (~1e-6 of the method error).

    Solvation accuracy: with 3ob/mio the GBSA parameters are GFN-xTB-fit, so
    solvation free energies are only qualitative and for small neutral solutes
    can be positive/non-monotonic in the solvent dielectric.
    """

    # Resolve the spin (electron-count guess when not given) so the calculation
    # and the analysis use the same multiplicity; S > 0 enables spin polarisation.
    if spin is None:
        electrons = round(float(sum(atoms.get_atomic_numbers())) - charge)
        spin = 0.0 if electrons % 2 == 0 else 0.5

    if spin_constants is None:
        spin_constants = SPIN_CONSTANTS_3OB

    spin_kwargs = _spin_kwargs(atoms, spin, spin_constants)

    # Implicit solvation (empty for the gas-phase default). Applied to both the
    # optimisation and the Hessian so the geometry and frequencies are consistent.
    solvation_kwargs = _solvation_kwargs(solvent, solvation_param_file)
    dispersion_kwargs = _dispersion_kwargs(dispersion)
    engine_kwargs = {**spin_kwargs, **solvation_kwargs, **dispersion_kwargs, **kwargs}

    with _run_in_directory(directory):
        # run geometry optimization
        geoopt = Geoopt(atoms=atoms, charge=charge, **engine_kwargs)
        potential_energy = geoopt.potential_energy()
        optimized_atoms = geoopt.read()

        # run hessian calculation
        Hessian(atoms=optimized_atoms, charge=charge, **engine_kwargs)

        # run normal mode calculation
        modes = Modes()
        frequencies = modes.wave_numbers

        # run thermo calculation on the optimized geometry directly (DFTB+ writes
        # geo_opt.gen, not geo_opt.xyz, so avoid the disk round-trip entirely)
        thermo = run_thermo(
            frequencies,
            atoms=optimized_atoms,
            temperature=temperature,
            pressure=pressure,
            energy=potential_energy,
            engine='dftb+',
            charge=charge,
            spin=spin,
            quasi_rrho=quasi_rrho,
        )

    return thermo


def xtb_thermo(
        atoms,
        temperature=298.15,
        pressure=101325,
        charge=0.0,
        spin=None,
        method="GFN2-xTB",
        directory=None,
        quasi_rrho=False,
        fmax=0.01,
    ):
    """
    Run the thermo pipeline with GFN-xTB (tblite): optimise, compute the
    Hessian/frequencies, and evaluate the thermochemistry.

    Needs no Slater-Koster files or spin constants -- GFN-xTB parameters are
    built in for H-Rn and open-shell systems are handled natively via the number
    of unpaired electrons.

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
    spin : float, optional
        Spin quantum number S. Defaults to the minimum-spin electron-count guess
        (even -> 0, odd -> 0.5). ``round(2*S)`` unpaired electrons are passed to
        xTB, so radicals run open-shell automatically.
    method : str
        GFN-xTB parametrisation, ``"GFN2-xTB"`` (default) or ``"GFN1-xTB"``.
    directory : str, optional
        Working directory to run in (created if needed). Defaults to the current
        directory.
    quasi_rrho : bool
        If True, use Grimme's quasi-RRHO vibrational entropy. Default False.
    fmax : float
        Force convergence threshold for the optimisation (eV/A). Default 0.01.

    Returns
    -------
    Thermo
        The thermo calculation object.

    Notes
    -----
    Reproducibility: with ``OMP_NUM_THREADS=1`` identical inputs give
    bit-identical energies on the same machine; the entropy is thread-invariant
    regardless. The optimiser stopping point is geometry/tolerance sensitive, so
    tighten ``fmax`` (e.g. 1e-3) for more environment-stable geometries.
    """

    # Resolve the spin the same way as the DFTB+ path so the calculation and the
    # analysis use the same multiplicity.
    if spin is None:
        electrons = round(float(sum(atoms.get_atomic_numbers())) - charge)
        spin = 0.0 if electrons % 2 == 0 else 0.5

    prepared = atoms.copy()
    prepared.info["charge"] = int(round(charge))
    prepared.info["spin"] = int(round(2.0 * float(spin)))

    with _run_in_directory(directory):
        optimized_atoms, potential_energy, frequencies = optimise_and_frequencies(
            prepared, xtb_calculator(method), fmax=fmax
        )

        thermo = run_thermo(
            frequencies,
            atoms=optimized_atoms,
            temperature=temperature,
            pressure=pressure,
            energy=potential_energy,
            engine="xtb",
            charge=charge,
            spin=spin,
            quasi_rrho=quasi_rrho,
        )

    return thermo


def xtb_cli_thermo(
        atoms,
        temperature=298.15,
        pressure=101325,
        charge=0.0,
        spin=None,
        method="GFN2-xTB",
        solvent=None,
        directory=None,
        quasi_rrho=False,
    ):
    """
    Run the thermo pipeline with the native ``xtb`` program (GFN-xTB).

    Unlike :func:`xtb_thermo` (in-process tblite, gas-phase), this uses the
    ``xtb`` binary, which natively supports open-shell (``--uhf``), charge
    (``--chrg``) and implicit solvation (``--alpb``) -- so charged radicals in
    solution are handled with method-consistent parameters. A single
    ``xtb --ohess`` optimises and computes the Hessian.

    Parameters
    ----------
    atoms : ase.Atoms
        The initial geometry.
    temperature : float
        Temperature in K. Default 298.15.
    pressure : float
        Pressure in Pa. Default 101325.
    charge : float
        System charge. Default 0.0.
    spin : float, optional
        Spin quantum number S. Defaults to the minimum-spin electron-count guess
        (even -> 0, odd -> 0.5); ``round(2*S)`` unpaired electrons are passed to
        xtb, so radicals run open-shell automatically.
    method : str
        GFN parametrisation: ``"GFN2-xTB"`` (default), ``"GFN1-xTB"`` or
        ``"GFN0-xTB"``.
    solvent : str, optional
        ALPB implicit-solvation solvent (e.g. ``"water"``). Native to xtb, so the
        parameters are method-consistent. Defaults to gas phase.
    directory : str, optional
        Working directory (created if needed). Defaults to the current directory.
    quasi_rrho : bool
        If True, use Grimme's quasi-RRHO vibrational entropy. Default False.

    Returns
    -------
    Thermo
        The thermo calculation object.
    """

    if spin is None:
        electrons = round(float(sum(atoms.get_atomic_numbers())) - charge)
        spin = 0.0 if electrons % 2 == 0 else 0.5
    unpaired = int(round(2.0 * float(spin)))

    with _run_in_directory(directory):
        optimized_atoms, potential_energy, frequencies = run_xtb(
            atoms, charge=charge, unpaired=unpaired, method=method, solvent=solvent,
        )

        thermo = run_thermo(
            frequencies,
            atoms=optimized_atoms,
            temperature=temperature,
            pressure=pressure,
            energy=potential_energy,
            engine="xtb",
            charge=charge,
            spin=spin,
            quasi_rrho=quasi_rrho,
        )

    return thermo
