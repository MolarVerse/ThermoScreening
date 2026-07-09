"""pKa (acid dissociation) from computed ``Thermo`` objects.

A pure post-processing helper combining the absolute Gibbs free energies
(``Thermo.total_EeGtot()``) of an acid and its conjugate base into a pKa, via
the "direct method" thermodynamic cycle ``HA(soln) -> A-(soln) + H+(soln)``.
Compute both species with the same engine and conditions -- in solution (e.g.
``solvent="water"``) -- then combine them here.
"""

import math

from ..utils.physicalConstants import PhysicalConstants
from ._units import HARTREE_TO_KCAL_PER_MOL

# Absolute free energy of the proton in aqueous solution at 298.15 K, in
# kcal/mol, referenced to the same 1 atm gas-phase-style standard state
# ThermoScreening's own RRHO thermochemistry uses (so it adds directly to
# Thermo.total_EeGtot() with no further standard-state correction):
#
#   G(H+, gas, 1 atm)      = -6.28 kcal/mol   (ideal monatomic gas, Sackur-
#                                               Tetrode; Bartmess, J. Phys.
#                                               Chem. 1994, 98, 6420)
# + dG_solv(H+, aq, 1 atm) = -264.0 kcal/mol  (Tissandier et al., J. Phys.
#                                               Chem. A 1998, 102, 7787, and
#                                               Kelly, Cramer & Truhlar,
#                                               J. Phys. Chem. B 2006, 110,
#                                               16066, recommend -265.9
#                                               kcal/mol at a 1 mol/L
#                                               gas-phase reference state;
#                                               converting to the 1 atm
#                                               reference used throughout here
#                                               subtracts RT ln(24.46) = 1.89
#                                               kcal/mol at 298.15 K -- see
#                                               their Section 2 and Eq. 1-3)
# = G(H+, aq, 1 atm)       = -270.28 kcal/mol
#
# Convention-dependent like SHE_ABSOLUTE_POTENTIAL in reactions.py; override
# reference_free_energy for a value calibrated to your method (see
# calibrate_proton_reference), which the literature recommends for
# quantitative accuracy.
PROTON_AQUEOUS_FREE_ENERGY_KCAL = -270.28

_R_KCAL_PER_MOL_K = PhysicalConstants["R"] / (PhysicalConstants["cal"] * 1000.0)


def _delta_g_kcal(acid, base):
    """G(base) - G(acid), in kcal/mol (Hartree per particle -> kcal/mol)."""
    return (base.total_EeGtot() - acid.total_EeGtot()) * HARTREE_TO_KCAL_PER_MOL


def pKa(acid, base, temperature=298.15, reference_free_energy=PROTON_AQUEOUS_FREE_ENERGY_KCAL):
    """
    pKa for ``HA(soln) -> A-(soln) + H+(soln)`` (the "direct method").

    ``pKa = dG / (R T ln 10)``, with ``dG = G(base) - G(acid) +
    reference_free_energy`` the deprotonation free energy including the
    proton's aqueous reference free energy (the proton has no electronic
    structure, so no engine can compute it directly).

    Parameters
    ----------
    acid, base : Thermo
        The acid (HA) and its conjugate base (A-), computed consistently (same
        engine, same conditions -- ideally in solution, e.g. ``solvent="water"``
        -- and open-shell/charged as appropriate).
    temperature : float
        Temperature in K. Should match the temperature ``acid``/``base`` were
        computed at. Default 298.15.
    reference_free_energy : float
        The proton's absolute aqueous free energy in kcal/mol. Defaults to
        :data:`PROTON_AQUEOUS_FREE_ENERGY_KCAL`; pass a value calibrated to
        your method (see :func:`calibrate_proton_reference`) for quantitative
        accuracy.

    Returns
    -------
    float
        The (dimensionless) pKa.

    Notes
    -----
    The raw "direct method" with a literature proton reference is known to
    have several-pKa-unit systematic error even at DFT+continuum-solvent
    levels (e.g. Ho & Coote, Theor. Chem. Acc. 2010, 125, 3); GFN-xTB/DFTB
    absolute pKa is expected to be considerably less accurate still. Use for
    **relative** comparisons among structurally similar acids, or calibrate
    ``reference_free_energy`` against one experimental pKa of a similar
    reference acid.
    """
    delta_g_kcal = _delta_g_kcal(acid, base) + reference_free_energy
    return delta_g_kcal / (_R_KCAL_PER_MOL_K * temperature * math.log(10))


def calibrate_proton_reference(acid, base, experimental_pKa, temperature=298.15):
    """
    The ``reference_free_energy`` that reproduces ``experimental_pKa`` for a
    reference acid/base pair.

    Solves :func:`pKa` for ``reference_free_energy`` instead of ``pKa``. Use
    the result as ``pKa(..., reference_free_energy=...)`` for other acids
    computed the same way (same engine/conditions), ideally structurally
    similar to the calibration pair -- the literature-recommended approach for
    quantitative pKa accuracy, since the raw literature proton reference alone
    is not quantitatively reliable (see :func:`pKa`'s Notes).

    Parameters
    ----------
    acid, base : Thermo
        A reference acid/base pair with a known experimental pKa.
    experimental_pKa : float
        The reference acid's experimental pKa.
    temperature : float
        Temperature in K, matching ``acid``/``base``. Default 298.15.

    Returns
    -------
    float
        The calibrated ``reference_free_energy`` (kcal/mol).
    """
    target_delta_g_kcal = experimental_pKa * _R_KCAL_PER_MOL_K * temperature * math.log(10)
    return target_delta_g_kcal - _delta_g_kcal(acid, base)
