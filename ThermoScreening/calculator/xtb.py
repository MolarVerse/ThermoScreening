"""xTB (GFN) calculator via tblite: geometry optimisation and vibrations.

Runs entirely in-process through the tblite ASE calculator (no external binary
and no Slater-Koster files -- GFN-xTB parameters are built in for H-Rn). The
geometry is optimised with an ASE optimiser and the frequencies come from ASE's
finite-difference vibrational analysis, so the result plugs into ``run_thermo``
exactly like the DFTB+ path.
"""

import numpy as np

from ..utils.physicalConstants import PhysicalConstants


def _eV_to_hartree(energy_ev):
    """
    Convert an energy in eV (ASE's unit) to Hartree (the thermo unit).
    """

    return energy_ev * PhysicalConstants["eV"] / PhysicalConstants["H"]


def _real_frequencies_cm(frequencies):
    """
    Convert ASE's complex vibrational frequencies to a sorted real array.

    ASE returns frequencies in cm^-1 as a complex array where imaginary modes
    carry the value in the imaginary part. They are mapped to negative real
    frequencies (the convention the ``System`` uses) and sorted ascending, so
    the translational/rotational and any imaginary modes sit at the bottom and
    ``System`` keeps the highest ``dof`` real vibrations.

    Parameters
    ----------
    frequencies : np.ndarray
        Complex vibrational frequencies in cm^-1 (from ``Vibrations``).

    Returns
    -------
    np.ndarray
        Sorted real frequencies in cm^-1 (imaginary modes negative).
    """
    frequencies = np.asarray(frequencies)
    real = np.where(
        np.abs(frequencies.imag) > 1e-6,
        -np.abs(frequencies.imag),
        frequencies.real,
    )
    return np.sort(real)


def optimise_and_frequencies(atoms, calc, fmax=0.01):
    """
    Optimise ``atoms`` with ``calc`` and return geometry, energy, frequencies.

    Parameters
    ----------
    atoms : ase.Atoms
        Initial geometry with total charge and unpaired electrons represented by
        ASE initial charges and magnetic moments for the calculator.
    calc : ase.calculators.calculator.Calculator
        The ASE calculator to attach (e.g. a tblite ``TBLite`` instance).
    fmax : float
        Force convergence threshold for the optimisation (eV/A).

    Returns
    -------
    tuple(ase.Atoms, float, np.ndarray)
        The optimised atoms, the energy in Hartree, and the sorted real
        vibrational frequencies in cm^-1.
    """
    from ase.optimize import BFGS
    from ase.vibrations import Vibrations

    atoms = atoms.copy()
    atoms.calc = calc

    BFGS(atoms, logfile=None).run(fmax=fmax)
    energy_hartree = _eV_to_hartree(atoms.get_potential_energy())

    vibrations = Vibrations(atoms, name="xtb_vib")
    vibrations.clean()
    try:
        vibrations.run()
        frequencies = _real_frequencies_cm(vibrations.get_frequencies())
    finally:
        vibrations.clean()

    return atoms, energy_hartree, frequencies


def xtb_calculator(method="GFN2-xTB"):
    """
    Build a tblite GFN-xTB ASE calculator.

    tblite reads total charge and unpaired electrons from the attached atoms'
    initial charges and magnetic moments, respectively.

    Parameters
    ----------
    method : str
        ``"GFN2-xTB"`` or ``"GFN1-xTB"``.
    """
    from tblite.ase import TBLite  # pragma: no cover - optional heavy dependency

    return TBLite(method=method, verbosity=0)  # pragma: no cover
