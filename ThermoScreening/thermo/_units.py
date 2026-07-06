"""Shared Hartree-based unit conversions for the post-processing helpers."""

from ..utils.physicalConstants import PhysicalConstants

HARTREE_TO_EV = PhysicalConstants["H"] / PhysicalConstants["eV"]
HARTREE_TO_KCAL_PER_MOL = (
    PhysicalConstants["H"] * PhysicalConstants["N_A"] / (PhysicalConstants["cal"] * 1000.0)
)
HARTREE_TO_KJ_PER_MOL = PhysicalConstants["H"] * PhysicalConstants["N_A"] / 1000.0

# Output units for energies given internally in Hartree (per particle).
# "eV" and "H" are per particle; "kcal" and "kJ" are per mole.
UNIT_FACTORS = {
    "H": 1.0,
    "eV": HARTREE_TO_EV,
    "kcal": HARTREE_TO_KCAL_PER_MOL,
    "kJ": HARTREE_TO_KJ_PER_MOL,
}


def convert_from_hartree(value, unit):
    """
    Convert an energy in Hartree to ``unit`` (one of ``H``/``eV``/``kcal``/``kJ``).

    Raises
    ------
    ValueError
        If ``unit`` is not supported.
    """
    try:
        factor = UNIT_FACTORS[unit]
    except KeyError:
        known = ", ".join(UNIT_FACTORS)
        raise ValueError(f"Unknown unit {unit!r}; choose one of: {known}.")
    return value * factor
