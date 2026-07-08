"""Import QM outputs (Gaussian, Turbomole, ORCA, Psi4, NWChem, ...) via cclib.

`cclib <https://cclib.github.io/>`_ auto-detects the program that produced an
output file and exposes a uniform data model, so a single reader consumes a
frequency calculation from any of the programs it supports. cclib is an optional
dependency; install it with ``pip install thermoscreening[qm]``.
"""

import numpy as np
from ase import Atoms

from ..exceptions import TSValueError

# 1 Hartree in eV (cclib 1.x reports energies in eV; the qm extra pins < 2.0).
_EV_PER_HARTREE = 27.211386245988


def _import_cclib():
    """Import cclib, or raise a helpful error if the optional dep is missing."""
    try:
        import cclib  # noqa: F401
        import cclib.io
    except ImportError as exc:  # pragma: no cover - exercised via monkeypatch
        raise TSValueError(
            "cclib is required to import QM outputs; install it with "
            "'pip install thermoscreening[qm]'."
        ) from exc
    return cclib


def _best_energy_ev(data):
    """
    The highest-level electronic energy available (eV), or ``None``.

    Prefers coupled-cluster, then Moller-Plesset, then SCF -- matching the level
    the Hessian was most likely computed at.
    """
    ccenergies = getattr(data, "ccenergies", None)
    if ccenergies is not None and len(ccenergies):
        return float(ccenergies[-1])
    mpenergies = getattr(data, "mpenergies", None)
    if mpenergies is not None and len(mpenergies):
        return float(mpenergies[-1][-1])  # last step, highest MP order
    scfenergies = getattr(data, "scfenergies", None)
    if scfenergies is not None and len(scfenergies):
        return float(scfenergies[-1])
    return None


def _normalize_source(path):
    """
    Coerce ``path`` into what ``cclib.io.ccread`` expects: a single path string,
    or a list of path strings for a multi-file program.
    """
    if isinstance(path, (list, tuple)):
        return [str(p) for p in path]
    return str(path)


def _describe_source(path):
    """A short, readable description of ``path`` for error messages."""
    if isinstance(path, (list, tuple)):
        return ", ".join(str(p) for p in path)
    return str(path)


def read_cclib(path):
    """
    Read geometry, vibrational frequencies and energy from a QM output file.

    Uses cclib to parse any supported program's output (Gaussian, ORCA, Psi4,
    NWChem, ...). Most programs write one logfile per job; pass its path.

    **Turbomole** splits a job's output across many small files instead of one
    logfile (``control``, ``coord``, ``aoforce.out``, ...); pass a list of every
    relevant file's path (cclib's own multi-file mode) rather than a single path.

    Parameters
    ----------
    path : str or list of str
        Path to a QM frequency-calculation output file, or (for a multi-file
        program such as Turbomole) a list of paths to every relevant file from
        the same job.

    Returns
    -------
    atoms : ase.Atoms
        The final geometry (cclib reports coordinates in Angstrom).
    frequencies : np.ndarray
        The vibrational frequencies in cm^-1 (cclib gives the real vibrational
        modes; imaginary modes appear as negatives).
    energy : float or None
        The electronic energy in Hartree (converted from cclib's eV), or ``None``
        if the output has no parseable energy.

    Notes
    -----
    cclib returns only the real vibrational modes (3N-6, or 3N-5 for a linear
    molecule). The thermochemistry keeps the top ``dof`` modes for the geometry's
    own linearity classification, so this matches for clearly (non)linear
    molecules; a near-linear geometry where the classification disagrees with the
    QM program's mode count could drop or miss a mode.

    Raises
    ------
    TSValueError
        If cclib is not installed, the file(s) cannot be parsed, or there are no
        vibrational frequencies.
    """
    cclib = _import_cclib()

    data = cclib.io.ccread(_normalize_source(path))
    if data is None:
        raise TSValueError(f"cclib could not parse '{_describe_source(path)}' as a QM output.")
    if getattr(data, "vibfreqs", None) is None or not len(data.vibfreqs):
        raise TSValueError(
            f"'{_describe_source(path)}' has no vibrational frequencies "
            "(run a frequency calculation)."
        )

    atoms = Atoms(numbers=np.asarray(data.atomnos), positions=np.asarray(data.atomcoords[-1]))
    frequencies = np.asarray(data.vibfreqs, dtype=float)

    energy_ev = _best_energy_ev(data)
    energy = energy_ev / _EV_PER_HARTREE if energy_ev is not None else None

    return atoms, frequencies, energy
