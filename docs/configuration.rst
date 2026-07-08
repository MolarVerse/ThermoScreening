Configuration
=============

Engines and methods
-------------------

.. list-table::
   :header-rows: 1
   :widths: 15 30 55

   * - ``engine``
     - Backend
     - Notes
   * - ``dftb+`` (default)
     - DFTB+
     - ``--parameter-set 3ob`` (DFTB3) or ``mio`` (DFTB2); needs ``DFTB_PREFIX``
   * - ``xtb``
     - GFN-xTB via tblite
     - in-process, gas phase only
   * - ``xtb-cli``
     - native ``xtb`` binary
     - open-shell + charge + implicit solvation

For the xtb engines, ``--method`` selects ``GFN2-xTB`` (default) or ``GFN1-xTB``.

Options
-------

``--solvent <name>``
    GBSA/ALPB implicit solvation (``water``, ``dmso``, ``acetonitrile``,
    ``thf``, ...). Applies to ``dftb+`` and ``xtb-cli``. With DFTB+ the GBSA
    parameters are GFN-xTB-fit, so solvation free energies are qualitative.

``--quasi-rrho``
    Use Grimme's quasi-RRHO vibrational entropy, which tames the entropy of
    low-frequency modes (recommended for flexible molecules).

``--charge`` / manifest ``charge`` column
    Total charge. Applied to directory input, or per-row in a CSV manifest.

Spin / multiplicity
    Resolved automatically from the electron count (even -> singlet, odd ->
    doublet); an explicit spin ``S`` can be given per manifest row (``spin``
    column) or via the Python API for cases like triplets. ``round(2*S)``
    unpaired electrons drive both the (spin-polarised) calculation and the
    electronic entropy ``R ln(2S+1)``.

``--temperature`` / ``--pressure``
    Thermodynamic conditions (default 298.15 K, 101325 Pa = 1 atm). Note that
    tabulated standard entropies use 1 bar (100000 Pa).

``--resume``
    Reuse the successful records from a prior ``<out>.json`` and only (re)run the
    missing or previously-failed molecules. Results are written incrementally.

Screening input
---------------

``thermo screen`` accepts either:

- a **directory** of ``.xyz``/``.gen`` structures (all run at ``--charge``), or
- a **CSV manifest** with a ``path`` column and optional ``name``, ``charge``
  and ``spin`` columns.

Results are written to ``<out>.csv`` and ``<out>.json`` with the electronic
energy (``Eelec_hartree``), the thermal corrections, the absolute Gibbs free
energy (``G_total_hartree``), entropy and heat capacity.
