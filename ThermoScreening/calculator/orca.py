"""Reader for ORCA Hessian (``.hess``) files.

Lets ThermoScreening consume DFT-quality geometries, vibrational frequencies and
energies from an ORCA frequency calculation, so the RRHO thermochemistry is
computed on accurate data instead of a semiempirical Hamiltonian.
"""

import re
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.units import Bohr

from ..exceptions import TSValueError

# $Single_Point_Data / &FinalEnergy line in an ORCA .property.txt file, e.g.:
#   &FinalEnergy [&Type "Double"]      -7.4963023139027911e+01  "Final single point energy"
_FINAL_ENERGY_RE = re.compile(r'&FinalEnergy\s*\[&Type\s*"Double"\]\s*([-+0-9.eE]+)')


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


def _read_property_energy(hess_path):
    """
    The final single-point energy (Hartree) from the companion ``.property.txt``
    file ORCA writes alongside a ``.hess`` (``$Single_Point_Data`` /
    ``&FinalEnergy``), or ``None`` if that file is absent or has no such entry.

    For a multi-step job (e.g. a relaxed scan) several ``&FinalEnergy`` entries
    can be present; the last one is the final/most converged energy.
    """
    property_path = Path(hess_path).with_suffix(".property.txt")
    if not property_path.is_file():
        return None
    text = property_path.read_text(encoding="utf-8")
    matches = _FINAL_ENERGY_RE.findall(text)
    if not matches:
        return None
    return float(matches[-1])


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
        The electronic energy in Hartree, preferably from the companion
        ``<basename>.property.txt`` file ORCA writes alongside the ``.hess``
        (its ``$Single_Point_Data`` / ``&FinalEnergy`` entry); falls back to a
        nonzero ``$act_energy`` in the ``.hess`` itself (populated for relaxed
        surface scans); ``None`` if neither is available.

    Notes
    -----
    ``$act_energy`` (along with ``$act_atom``/``$act_coord``) is an ORCA
    relaxed-surface-scan field for the *active* point on the scan; for an
    ordinary optimization + frequency job it is ``0.0`` and does not carry the
    electronic energy at all, so it is not read unless it is genuinely nonzero.

    The per-atom mass column is not used; ASE's standard isotope masses are used
    for the rotational/translational terms (as with the other engines).

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

    energy = _read_property_energy(path)
    if energy is None:
        energy_block = _read_block(lines, "act_energy")
        if energy_block is not None:
            for line in energy_block:
                if line.strip():
                    # $act_energy is a relaxed-scan field: 0.0 for an ordinary
                    # job means "not set", not a real zero-Hartree energy.
                    act_energy = float(line.split()[0])
                    if act_energy != 0.0:
                        energy = act_energy
                    break

    return atoms, frequencies, energy
