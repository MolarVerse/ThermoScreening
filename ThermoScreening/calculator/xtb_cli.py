"""Native xTB (GFN) engine via the ``xtb`` command-line program.

Unlike the in-process tblite route (:mod:`ThermoScreening.calculator.xtb`), the
native ``xtb`` binary exposes open-shell (``--uhf``), charge (``--chrg``) and
implicit solvation (``--alpb``), so it covers charged radicals in solution. A
single ``xtb --ohess`` call optimises the geometry and computes the Hessian; the
energy, optimised geometry and frequencies are then read back and fed into
``run_thermo`` exactly like the other engines.
"""

import os
import shutil
import subprocess

import numpy as np
import ase.io

# GFN method name -> ``--gfn`` argument
_GFN_METHODS = {"GFN2-xTB": "2", "GFN1-xTB": "1", "GFN0-xTB": "0"}


def resolve_xtb(command=None):
    """
    Locate the ``xtb`` executable.

    Order: explicit ``command``, then ``$XTB_COMMAND``, then ``xtb`` on PATH.

    Raises
    ------
    FileNotFoundError
        If no ``xtb`` executable can be found.
    """
    candidate = command or os.environ.get("XTB_COMMAND") or shutil.which("xtb")
    if not candidate:
        raise FileNotFoundError(
            "The 'xtb' executable was not found. Install it (e.g. "
            "`conda install -c conda-forge xtb`) or set the XTB_COMMAND "
            "environment variable to its path."
        )
    return candidate


def _parse_optimised_energy(path="xtbopt.xyz"):
    """
    Read the optimised total energy (Hartree) from an ``xtbopt.xyz`` comment.

    The second line looks like ``energy: -5.070544 gnorm: ... xtb: 6.7.1``.
    """
    with open(path, "r", encoding="utf-8") as handle:
        handle.readline()  # atom count
        comment = handle.readline()

    tokens = comment.split()
    for index, token in enumerate(tokens):
        if token == "energy:":
            return float(tokens[index + 1])
    raise ValueError(f"No energy found in {path!r} comment line: {comment!r}")


def _parse_vibspectrum(path="vibspectrum"):
    """
    Parse frequencies (cm^-1) from a Turbomole ``vibspectrum`` file.

    Each data line is ``mode [symmetry] wavenumber IR-intensity selection``; the
    wavenumber is the third-from-last token whether or not a symmetry label is
    present. Imaginary modes are negative; translational/rotational modes are
    ~0. The result is sorted ascending so ``System`` keeps the highest ``dof``
    real vibrations (as with the other engines).
    """
    frequencies = []
    in_block = False
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if stripped.startswith("$vibrational"):
                in_block = True
                continue
            if stripped.startswith("$end"):
                break
            if not in_block or stripped.startswith("#") or not stripped:
                continue
            tokens = stripped.split()
            if len(tokens) >= 4 and tokens[0].isdigit():
                frequencies.append(float(tokens[-3]))
    return np.sort(np.array(frequencies))


def run_xtb(
    atoms,
    charge=0.0,
    unpaired=0,
    method="GFN2-xTB",
    solvent=None,
    command=None,
):
    """
    Run ``xtb --ohess`` and return the optimised geometry, energy, frequencies.

    Runs in the current working directory (xtb writes several files there), so
    call it inside a per-job directory.

    Parameters
    ----------
    atoms : ase.Atoms
        The initial geometry.
    charge : float
        Total charge (``--chrg``).
    unpaired : int
        Number of unpaired electrons, i.e. round(2*S) (``--uhf``).
    method : str
        ``"GFN2-xTB"`` (default), ``"GFN1-xTB"`` or ``"GFN0-xTB"``.
    solvent : str, optional
        ALPB implicit-solvation solvent (``--alpb``), e.g. ``"water"``.
    command : str, optional
        Path to the ``xtb`` executable (defaults to PATH / ``$XTB_COMMAND``).

    Returns
    -------
    tuple(ase.Atoms, float, np.ndarray)
        Optimised atoms, energy in Hartree, sorted real frequencies in cm^-1.

    Raises
    ------
    ValueError
        If ``method`` is unknown.
    RuntimeError
        If the xtb run fails.
    """
    try:
        gfn = _GFN_METHODS[method]
    except KeyError:
        known = ", ".join(sorted(_GFN_METHODS))
        raise ValueError(f"Unknown xTB method {method!r}; choose one of: {known}.")

    executable = resolve_xtb(command)
    ase.io.write("xtb_input.xyz", atoms)

    argv = [
        executable, "xtb_input.xyz", "--ohess", "--gfn", gfn,
        "--chrg", str(int(round(charge))), "--uhf", str(int(unpaired)),
    ]
    if solvent is not None:
        argv += ["--alpb", solvent]

    result = subprocess.run(argv, capture_output=True, text=True, check=False)
    if result.returncode != 0 or not os.path.isfile("xtbopt.xyz"):
        tail = (result.stdout or "")[-800:] + (result.stderr or "")[-800:]
        raise RuntimeError(f"xtb failed (exit {result.returncode}):\n{tail}")

    optimized = ase.io.read("xtbopt.xyz")
    energy_hartree = _parse_optimised_energy("xtbopt.xyz")
    frequencies = _parse_vibspectrum("vibspectrum")
    return optimized, energy_hartree, frequencies


def _parse_fukui(output):
    """
    Parse the per-atom ``Fukui functions:`` table from ``xtb --vfukui`` stdout.

    The table looks like::

        Fukui functions:
             #        f(+)     f(-)     f(0)
             1C       0.024    0.024    0.024
             8O       0.131    0.129    0.130

    Returns a list of ``(symbol, f_plus, f_minus, f_zero)`` tuples, in the
    input geometry's atom order.

    Raises
    ------
    ValueError
        If no ``Fukui functions:`` table is found in ``output``.
    """
    lines = output.splitlines()
    for index, line in enumerate(lines):
        if line.strip() == "Fukui functions:":
            break
    else:
        raise ValueError("No 'Fukui functions:' table found in xtb output.")

    rows = []
    for line in lines[index + 2:]:
        tokens = line.split()
        if len(tokens) != 4:
            break
        label, f_plus, f_minus, f_zero = tokens
        symbol = label.lstrip("0123456789")
        rows.append((symbol, float(f_plus), float(f_minus), float(f_zero)))
    return rows


def run_xtb_fukui(
    atoms,
    charge=0.0,
    unpaired=0,
    method="GFN2-xTB",
    solvent=None,
    command=None,
):
    """
    Run ``xtb --vfukui`` and return per-atom Fukui reactivity indices.

    A single-point (non-geometry-optimising) calculation: pass an already
    optimised geometry (e.g. from :func:`run_xtb`). Runs in the current
    working directory, so call it inside a per-job directory.

    Parameters
    ----------
    atoms : ase.Atoms
        The geometry to analyse.
    charge : float
        Total charge (``--chrg``).
    unpaired : int
        Number of unpaired electrons, i.e. round(2*S) (``--uhf``).
    method : str
        ``"GFN2-xTB"`` (default), ``"GFN1-xTB"`` or ``"GFN0-xTB"``.
    solvent : str, optional
        ALPB implicit-solvation solvent (``--alpb``), e.g. ``"water"``.
    command : str, optional
        Path to the ``xtb`` executable (defaults to PATH / ``$XTB_COMMAND``).

    Returns
    -------
    list of tuple(str, float, float, float)
        ``(symbol, f_plus, f_minus, f_zero)`` per atom, in input geometry order.
        ``f_plus`` (susceptibility to nucleophilic attack / electron gain),
        ``f_minus`` (electrophilic attack / electron loss), ``f_zero`` (radical
        attack, the average of the two).

    Raises
    ------
    ValueError
        If ``method`` is unknown, or xtb's output has no Fukui table.
    RuntimeError
        If the xtb run fails.
    """
    try:
        gfn = _GFN_METHODS[method]
    except KeyError:
        known = ", ".join(sorted(_GFN_METHODS))
        raise ValueError(f"Unknown xTB method {method!r}; choose one of: {known}.")

    executable = resolve_xtb(command)
    ase.io.write("xtb_input.xyz", atoms)

    argv = [
        executable, "xtb_input.xyz", "--vfukui", "--gfn", gfn,
        "--chrg", str(int(round(charge))), "--uhf", str(int(unpaired)),
    ]
    if solvent is not None:
        argv += ["--alpb", solvent]

    result = subprocess.run(argv, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        tail = (result.stdout or "")[-800:] + (result.stderr or "")[-800:]
        raise RuntimeError(f"xtb failed (exit {result.returncode}):\n{tail}")

    return _parse_fukui(result.stdout)
