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
energy (``G_total_hartree``), entropy and heat capacity. ``<out>-run.json``
records the structure hashes and scientific settings used for safe resume.
Automatic Gibbs ordering is shown only when every structure has the same
molecular formula and charge. Absolute energies do not rank mixed molecular
compositions.

Cluster execution
-----------------

Batch screening can be divided into deterministic, zero-based shards on any
scheduler. Every shard writes isolated result files under
``<out>-shards/shard-NNNNN.*`` and calculations under
``<directory>/shards/shard-NNNNN``:

.. code-block:: bash

   thermo screen molecules.csv -o results \
       --shard-index "$TASK_INDEX" --shard-count 32

After every shard has finished, validate and combine them in original input
order:

.. code-block:: bash

   thermo collect results-shards -o results

Collection fails on missing or duplicate shards, mixed settings, missing jobs,
or fingerprint mismatches. Scientific calculation failures remain in the
combined result and produce a non-zero command exit status.

For Slurm, one command can generate an array script and optionally submit both
the array and a dependent collection job:

.. code-block:: bash

   thermo slurm --tasks 32 --cpus-per-task 8 \
       --time 04:00:00 --mem 16G --partition compute --submit -- \
       screen molecules.csv -o results --engine dftb+ --jobs 2

The same interface distributes the complete redox workflow by candidate:

.. code-block:: bash

   thermo slurm --tasks 32 --cpus-per-task 8 --submit -- \
       redox molecules.csv -o potentials --engine xtb-cli --jobs 2

``--cpus-per-task`` describes each allocation; ``--jobs`` controls independent
processes inside that allocation. The generated script divides available CPUs
between those processes through ``OMP_NUM_THREADS`` and common BLAS thread
variables. Use ``--preamble setup.sh`` for site-specific module loads or
environment variables. Omitting ``--submit`` only writes the script.

All array tasks and the collector require the same shared working directory,
Python environment, executables, and parameter files. The shard interface also
works with PBS, LSF, and other schedulers by supplying their array index and the
chosen shard count directly. Both ``thermo screen`` and ``thermo redox`` expose
``--shard-index`` and ``--shard-count`` for this scheduler-neutral path.

Redox workflow input
--------------------

``thermo redox`` accepts a SMILES string, one ``.xyz``/``.gen`` structure, a
directory of structures, or a CSV manifest. Each manifest row needs exactly one
of ``path`` or ``smiles``:

.. code-block:: text

   name,path,smiles,charge,oxidized_spin,reduced_once_spin,reduced_twice_spin
   AQ,aq.xyz,,0,0,0.5,0
   hydroxy-AQ,,OC1=CC2=C(C=C1)C(=O)C1=CC=CC(=O)C1=C2,0,0,0.5,0

``name`` and ``charge`` are optional. Spin columns are optional and otherwise
inferred from electron count for each state. Relative paths are resolved from
the manifest directory. Molecular charges must be integers and spins must be
non-negative integers or half-integers. Ionic SMILES supply their formal charge;
structure files without a charge default to zero.

For SMILES input, ``--max-conformers`` controls how many ETKDG structures are
embedded. The lowest-MMFF-energy structure becomes the common starting geometry;
each charge state is then optimized independently by the selected backend.
This is a single-starting-geometry screening approximation, not a charge-state
conformer search or a conformational ensemble free energy. Flexible molecules
require separate conformer sampling for every charge state.

The workflow writes:

- ``<out>.csv`` and ``<out>.json`` with E1, E2, E2e, the potential gap, and the
  three absolute Gibbs free energies;
- ``<out>-states.csv`` and ``<out>-states.json`` with the complete per-state
  thermochemistry; and
- ``<out>-run.json`` and ``<out>-states-run.json`` with input hashes,
  calculation settings, reference calibration, and reproducibility fingerprints;
- ``<directory>/states.csv`` as the generated, reproducible charge-state
  manifest.

``--jobs`` is the number of independent charge-state calculations executed in
parallel. A screen of *N* molecules contains ``3 * N`` jobs, plus three only when
a separately supplied reference is not already part of the candidate set.
``--resume`` reuses only successful states whose structure content, charge,
spin, engine, method, solvent, thermodynamic conditions, and other scientific
settings match the stored fingerprint. Missing, failed, legacy, or mismatched
states are rerun. On multi-core backends, set ``OMP_NUM_THREADS=1`` when using
several process jobs to avoid oversubscribing CPU cores. Launch process-based
parallel runs from the CLI or a script protected by ``if __name__ == "__main__"``;
use ``--jobs 1`` in notebooks and other interactive sessions.
