"""Reader for ORCA Hessian (``.hess``) files.

Lets ThermoScreening consume DFT-quality geometries, vibrational frequencies and
energies from an ORCA frequency calculation, so the RRHO thermochemistry is
computed on accurate data instead of a semiempirical Hamiltonian.
"""

import numpy as np
from ase import Atoms
from ase.units import Bohr

from ..exceptions import TSValueError


def _read_block(lines, name):
    """
    Return the lines of the ``$name`` block (between its tag and the next ``$``
    tag or end of file), or ``None`` if the block is absent.
    """
    tag = f"${name}"
    for index, line in enumerate(lines):
        if line.strip() == tag:
            block = []
            for following in lines[index + 1:]:
                if following.lstrip().startswith("$"):
                    break
                block.append(following)
            return block
    return None


def read_orca_hess(path):
    """
    Read geometry, vibrational frequencies and energy from an ORCA ``.hess`` file.

    Parameters
    ----------
    path : str
        Path to an ORCA ``.hess`` file (from an ORCA frequency calculation).

    Returns
    -------
    atoms : ase.Atoms
        The geometry, converted from Bohr to Angstrom.
    frequencies : np.ndarray
        The vibrational frequencies in cm^-1, including the near-zero
        translational/rotational modes as ORCA writes them.
    energy : float or None
        The electronic energy in Hartree from the ``$act_energy`` block, or
        ``None`` if the file has no such block.

    Raises
    ------
    TSValueError
        If the required ``$atoms`` or ``$vibrational_frequencies`` block is
        missing or malformed.
    """
    with open(path, "r", encoding="utf-8") as handle:
        lines = handle.readlines()

    atoms_block = _read_block(lines, "atoms")
    if atoms_block is None:
        raise TSValueError(f"No $atoms block in ORCA hess file '{path}'.")
    freq_block = _read_block(lines, "vibrational_frequencies")
    if freq_block is None:
        raise TSValueError(
            f"No $vibrational_frequencies block in ORCA hess file '{path}'."
        )

    # $atoms: <natoms>, then per atom "<symbol> <mass> <x> <y> <z>" (Bohr)
    atom_lines = [line for line in atoms_block if line.strip()]
    try:
        n_atoms = int(atom_lines[0].split()[0])
        symbols, positions = [], []
        for line in atom_lines[1:1 + n_atoms]:
            tokens = line.split()
            symbols.append(tokens[0])
            positions.append([float(tokens[2]), float(tokens[3]), float(tokens[4])])
    except (IndexError, ValueError) as exc:
        raise TSValueError(f"Malformed $atoms block in '{path}': {exc}")
    if len(symbols) != n_atoms:
        raise TSValueError(
            f"$atoms block in '{path}' declares {n_atoms} atoms but lists {len(symbols)}."
        )
    atoms = Atoms(symbols=symbols, positions=np.array(positions) * Bohr)

    # $vibrational_frequencies: <count>, then "<index> <frequency>" (cm^-1)
    freq_lines = [line for line in freq_block if line.strip()]
    try:
        n_freq = int(freq_lines[0].split()[0])
        frequencies = np.array(
            [float(line.split()[1]) for line in freq_lines[1:1 + n_freq]]
        )
    except (IndexError, ValueError) as exc:
        raise TSValueError(
            f"Malformed $vibrational_frequencies block in '{path}': {exc}"
        )
    if len(frequencies) != n_freq:
        raise TSValueError(
            f"$vibrational_frequencies block in '{path}' declares {n_freq} entries "
            f"but lists {len(frequencies)}."
        )

    energy = None
    energy_block = _read_block(lines, "act_energy")
    if energy_block is not None:
        for line in energy_block:
            if line.strip():
                energy = float(line.split()[0])
                break

    return atoms, frequencies, energy
