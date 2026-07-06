"""Conformer-ensemble thermochemistry from computed ``Thermo`` objects.

Pure post-processing helpers that combine the absolute Gibbs free energies
(``Thermo.total_EeGtot()`` = electronic + thermal, in Hartree) of several
conformers of the *same* molecule into ensemble properties: Boltzmann
populations, the Boltzmann-averaged (ensemble) free energy, and the
lowest-free-energy conformer. Compute each conformer with the same engine and
conditions (see :func:`ThermoScreening.thermo.conformers.generate`), then
combine them here.
"""

import math

from ..utils.physicalConstants import PhysicalConstants
from ._units import convert_from_hartree

# Boltzmann constant in Hartree per kelvin (k_B [J/K] converted to Hartree).
_BOLTZMANN_HARTREE_PER_K = PhysicalConstants["kB"] / PhysicalConstants["H"]


def _energies_and_kt(thermos, temperature):
    """Return the list of conformer free energies (Hartree) and k_B*T (Hartree)."""
    if not thermos:
        raise ValueError("thermos must contain at least one conformer.")
    if temperature <= 0:
        raise ValueError("temperature must be positive.")
    energies = [thermo.total_EeGtot() for thermo in thermos]
    return energies, _BOLTZMANN_HARTREE_PER_K * temperature


def boltzmann_weights(thermos, temperature=298.15):
    """
    Boltzmann populations of a conformer ensemble at ``temperature`` (K).

    Parameters
    ----------
    thermos : iterable of Thermo
        Conformers of the same molecule, computed with matching settings.
    temperature : float
        Temperature in kelvin. Default 298.15.

    Returns
    -------
    list of float
        Normalised populations (summing to 1) in the order of ``thermos``.

    Raises
    ------
    ValueError
        If ``thermos`` is empty or ``temperature`` is not positive.
    """
    thermos = list(thermos)
    energies, kt = _energies_and_kt(thermos, temperature)
    # Shift by the minimum for numerical stability (the ratios are unchanged).
    reference = min(energies)
    boltzmann = [math.exp(-(energy - reference) / kt) for energy in energies]
    partition = sum(boltzmann)
    return [weight / partition for weight in boltzmann]


def ensemble_free_energy(thermos, temperature=298.15, unit="H"):
    """
    Boltzmann-averaged (ensemble) free energy of a conformer set.

    ``G_ensemble = -k_B T ln( sum_i exp(-G_i / k_B T) )``, which lies at or below
    the lowest conformer free energy by the mixing (conformational) entropy.

    Parameters
    ----------
    thermos : iterable of Thermo
        Conformers of the same molecule, computed with matching settings.
    temperature : float
        Temperature in kelvin. Default 298.15.
    unit : str
        Output unit: ``"H"`` (Hartree per particle, default), ``"eV"`` (per
        particle), ``"kcal"`` (kcal/mol), or ``"kJ"`` (kJ/mol).

    Returns
    -------
    float
        The ensemble free energy in the requested unit.

    Raises
    ------
    ValueError
        If ``thermos`` is empty, ``temperature`` is not positive, or ``unit`` is
        not supported.
    """
    thermos = list(thermos)
    energies, kt = _energies_and_kt(thermos, temperature)
    reference = min(energies)
    partition = sum(math.exp(-(energy - reference) / kt) for energy in energies)
    free_energy = reference - kt * math.log(partition)
    return convert_from_hartree(free_energy, unit)


def lowest_gibbs(thermos):
    """
    Return the conformer with the lowest absolute Gibbs free energy.

    Parameters
    ----------
    thermos : iterable of Thermo
        Conformers of the same molecule, computed with matching settings.

    Returns
    -------
    Thermo
        The lowest-free-energy conformer.

    Raises
    ------
    ValueError
        If ``thermos`` is empty.
    """
    thermos = list(thermos)
    if not thermos:
        raise ValueError("thermos must contain at least one conformer.")
    return min(thermos, key=lambda thermo: thermo.total_EeGtot())
