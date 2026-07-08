"""Transition-state-theory rate constants from computed ``Thermo`` objects.

Pure post-processing helpers that combine the absolute Gibbs free energies
(``Thermo.total_EeGtot()``) of a reactant and a transition state
(``Thermo(..., transition_state=True)``, see :mod:`ThermoScreening.thermo.thermo`)
into an Eyring rate constant, plus a simple tunneling correction.
"""

import math

from ..utils.physicalConstants import PhysicalConstants
from .reactions import _total_gibbs


def eyring_rate_constant(reactants, ts, temperature=298.15, kappa=1.0):
    """
    Eyring transition-state-theory rate constant.

    ``k = kappa * (kB T / h) * exp(-dG-double-dagger / (R T))``, with the
    activation free energy ``dG-double-dagger = G(ts) - sum(G(reactants))``.

    Parameters
    ----------
    reactants : iterable
        The reactant species, in the same convention as
        :func:`ThermoScreening.thermo.reactions.reaction_free_energy`: each
        entry is a ``Thermo`` (stoichiometric coefficient 1) or a
        ``(coefficient, Thermo)`` tuple. Pass a single-entry list for a
        unimolecular reaction.
    ts : Thermo
        The transition state, computed with ``transition_state=True`` (see
        :class:`ThermoScreening.thermo.thermo.Thermo`).
    temperature : float
        Temperature in K, matching the temperature ``reactants`` and ``ts``
        were computed at. Default 298.15.
    kappa : float
        Transmission coefficient (a tunneling/recrossing correction). Default
        1.0 (no correction); see :func:`wigner_tunneling_correction` for a
        simple estimate from the transition state's imaginary frequency.

    Returns
    -------
    float
        The TST rate constant. For a single reactant this is a first-order
        rate constant in s^-1. For multiple reactants (a bimolecular or higher
        reaction) this is the pseudo-first-order rate assuming each reactant is
        at its computed standard state; converting to a true bimolecular+ rate
        constant (e.g. L mol^-1 s^-1) requires an additional standard-state
        correction that is not applied here.
    """
    delta_g_ddagger = ts.total_EeGtot() - sum(
        _total_gibbs(reactant) for reactant in reactants
    )  # Hartree per particle
    delta_g_j_per_mol = delta_g_ddagger * PhysicalConstants["H"] * PhysicalConstants["N_A"]
    return (
        kappa
        * (PhysicalConstants["kB"] * temperature / PhysicalConstants["h"])
        * math.exp(-delta_g_j_per_mol / (PhysicalConstants["R"] * temperature))
    )


def wigner_tunneling_correction(imaginary_wavenumber, temperature=298.15):
    """
    Wigner's first-order tunneling correction to a TST rate constant.

    ``kappa = 1 + (1/24) * (h c |nu_imag| / (kB T))^2`` (Wigner, Z. Phys. Chem.
    B 19, 203 (1932)), a small quantum correction for tunneling along the
    reaction coordinate. Valid only for small corrections; for deep tunneling
    (large kappa) use a more complete treatment (e.g. Eckart).

    Parameters
    ----------
    imaginary_wavenumber : float
        The transition state's imaginary-mode wavenumber in cm^-1 (as returned
        by ``Thermo.imaginary_mode_wavenumber()``; the sign is ignored).
    temperature : float
        Temperature in K. Default 298.15.

    Returns
    -------
    float
        The dimensionless transmission coefficient kappa (>= 1), for use as
        ``eyring_rate_constant(..., kappa=...)``.
    """
    u = (
        PhysicalConstants["h"]
        * PhysicalConstants["c"]
        * abs(imaginary_wavenumber)
        * 10**2  # cm^-1 -> m^-1
        / (PhysicalConstants["kB"] * temperature)
    )
    return 1.0 + u**2 / 24.0
