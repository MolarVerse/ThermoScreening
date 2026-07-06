"""Reaction and redox thermochemistry from computed ``Thermo`` objects.

These are pure post-processing helpers that combine the absolute Gibbs free
energies (``Thermo.total_EeGtot()`` = electronic + thermal) of individual species
into reaction free energies and reduction potentials. Compute each species with
any engine (``dftbplus_thermo`` / ``xtb_thermo`` / ``xtb_cli_thermo``) under the
same conditions -- ideally the same solvent, and open-shell/charged as
appropriate -- then combine them here.
"""

from ._units import HARTREE_TO_EV as _HARTREE_TO_EV, convert_from_hartree

# Absolute potential of the standard hydrogen electrode (V). Convention-dependent
# (values from ~4.28 to ~4.44 V are used); this is the IUPAC-recommended value.
# Override ``reference_potential`` for a calibrated reference.
SHE_ABSOLUTE_POTENTIAL = 4.44


def _total_gibbs(species):
    """
    Absolute Gibbs free energy (Hartree) of a species entry.

    ``species`` is a ``Thermo`` (coefficient 1) or a ``(coefficient, Thermo)``
    tuple.
    """
    if isinstance(species, tuple):
        coefficient, thermo = species
    else:
        coefficient, thermo = 1.0, species
    return coefficient * thermo.total_EeGtot()


def reaction_free_energy(reactants, products, unit="kcal"):
    """
    Reaction free energy dG = sum(products) - sum(reactants).

    Parameters
    ----------
    reactants, products : iterable
        Each entry is a ``Thermo`` (stoichiometric coefficient 1) or a
        ``(coefficient, Thermo)`` tuple. All species should be computed with the
        same engine and conditions for the difference to be meaningful.
    unit : str
        Output unit: ``"kcal"`` (kcal/mol, default), ``"kJ"`` (kJ/mol), ``"eV"``
        (per particle), or ``"H"`` (Hartree per particle).

    Returns
    -------
    float
        The reaction free energy in the requested unit.

    Raises
    ------
    ValueError
        If ``unit`` is not supported.
    """
    delta = (
        sum(_total_gibbs(product) for product in products)
        - sum(_total_gibbs(reactant) for reactant in reactants)
    )
    return convert_from_hartree(delta, unit)


def reduction_potential(
    oxidized, reduced, n_electrons=1, reference_potential=SHE_ABSOLUTE_POTENTIAL
):
    """
    Reduction potential (V) for ``Ox + n e- -> Red``.

    Uses ``E = -dG_red / (n F) - E_reference`` with ``dG_red = G(Red) - G(Ox)``
    and the free energy of the electron taken as zero. Because one electron-volt
    per electron is one volt, the absolute potential simplifies to
    ``-dG_red[eV] / n``, which is then referenced to ``reference_potential``.

    Parameters
    ----------
    oxidized, reduced : Thermo
        The oxidised and reduced species, computed consistently (same engine,
        same solvent, and open-shell/charged as appropriate).
    n_electrons : int
        Number of electrons transferred. Default 1.
    reference_potential : float
        Absolute potential (V) of the reference electrode to report against.
        Defaults to the SHE (:data:`SHE_ABSOLUTE_POTENTIAL`); pass ``0.0`` for the
        absolute reduction potential.

    Returns
    -------
    float
        The reduction potential in volts (versus ``reference_potential``).

    Raises
    ------
    ValueError
        If ``n_electrons`` is zero.
    """
    if n_electrons == 0:
        raise ValueError("n_electrons must be non-zero.")

    delta_g_hartree = reduced.total_EeGtot() - oxidized.total_EeGtot()
    absolute = -delta_g_hartree * _HARTREE_TO_EV / n_electrons
    return absolute - reference_potential
