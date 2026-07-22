"""Conformer generation for molecule screening.

RDKit is the conformer backend: ETKDG distance-geometry embedding, optional
MMFF force-field optimisation, and RMSD + energy-window pruning. Conformers are
returned as ASE ``Atoms`` so they feed straight into the screening pipeline
(e.g. generate -> ``write_conformers`` -> ``screen`` that directory).
"""

import numpy as np

from ..exceptions import TSValueError


def _import_rdkit():
    """Import RDKit lazily, with a clear message if it is not installed."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError as exc:  # pragma: no cover - exercised only without rdkit
        raise ImportError(
            "Conformer generation requires RDKit. Install it with "
            "`conda install -c conda-forge rdkit` or `pip install rdkit`."
        ) from exc
    return Chem, AllChem


def _conformer_to_atoms(mol, conformer_id):
    from ase import Atoms

    conformer = mol.GetConformer(conformer_id)
    symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
    positions = np.asarray(conformer.GetPositions(), dtype=float)
    atoms = Atoms(symbols=symbols, positions=positions)
    atoms.info["formal_charge"] = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    return atoms


def generate(
    smiles,
    max_conformers=10,
    optimize=True,
    prune_rms_thresh=0.5,
    energy_window=None,
    random_seed=42,
):
    """
    Generate conformers for a molecule from its SMILES.

    Parameters
    ----------
    smiles : str
        The molecule as a SMILES string.
    max_conformers : int
        Maximum number of conformers to embed (ETKDG). Default 10.
    optimize : bool
        If True, optimise each conformer with the MMFF force field and sort the
        results by MMFF energy (lowest first). Default True.
    prune_rms_thresh : float
        RMSD threshold (A) for pruning duplicate embedded conformers. Default 0.5.
    energy_window : float, optional
        If given (and ``optimize``), keep only conformers within this many
        kcal/mol of the lowest-energy conformer.
    random_seed : int
        Random seed for reproducible embedding. Default 42.

    Returns
    -------
    list of ase.Atoms
        The generated conformers (energy-sorted when optimised).

    Raises
    ------
    ImportError
        If RDKit is not installed.
    ValueError
        If the SMILES cannot be parsed or no conformer could be embedded.
    """
    Chem, AllChem = _import_rdkit()

    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"Could not parse SMILES: {smiles!r}")
    molecule = Chem.AddHs(molecule)

    params = AllChem.ETKDGv3()
    params.randomSeed = random_seed
    params.pruneRmsThresh = prune_rms_thresh
    conformer_ids = list(AllChem.EmbedMultipleConfs(molecule, numConfs=max_conformers, params=params))
    if not conformer_ids:
        raise ValueError(f"Could not embed any conformer for SMILES: {smiles!r}")

    if optimize:
        results = AllChem.MMFFOptimizeMoleculeConfs(molecule)
        energies = [energy for _converged, energy in results]
        order = sorted(range(len(conformer_ids)), key=lambda i: energies[i])
        conformer_ids = [conformer_ids[i] for i in order]
        energies = [energies[i] for i in order]
        if energy_window is not None:
            lowest = energies[0]
            conformer_ids = [
                cid for cid, energy in zip(conformer_ids, energies)
                if energy - lowest <= energy_window
            ]

    return [_conformer_to_atoms(molecule, cid) for cid in conformer_ids]


def generate_thermo_ensemble(
    smiles,
    thermo_fn,
    charge=0.0,
    max_conformers=10,
    max_attempts=5,
    prune_rms_thresh=0.5,
    energy_window=None,
    random_seed=42,
    **thermo_kwargs,
):
    """
    Generate a conformer ensemble and compute ``Thermo`` for each, retrying
    past saddle points.

    A conformer that fails to converge to a true minimum (``TSValueError``,
    e.g. an imaginary frequency) or whose engine process fails outright
    (``RuntimeError``, e.g. an xtb CLI non-zero exit) is skipped; if fewer
    than ``max_conformers`` succeed, :func:`generate` is re-run with a new
    random seed, up to ``max_attempts`` seeds total. This is the loop a
    charged-species ensemble (radical anion, dianion, ...) otherwise
    requires by hand before it can be passed to :class:`EnsembleThermo`.

    Parameters
    ----------
    smiles : str
        The molecule as a SMILES string.
    thermo_fn : callable
        A ``Thermo``-computing engine, e.g. ``xtb_cli_thermo``, called as
        ``thermo_fn(atoms, charge=charge, **thermo_kwargs)``.
    charge : float
        System charge, forwarded to ``thermo_fn``. Default 0.0.
    max_conformers : int
        Target number of successfully-computed conformers. Default 10.
    max_attempts : int
        Number of distinct random seeds (each embedding up to
        ``max_conformers`` conformers) to try before giving up. Default 5.
    prune_rms_thresh, energy_window : see :func:`generate`.
    random_seed : int
        Random seed for the first attempt; later attempts use
        ``random_seed + 1``, ``random_seed + 2``, etc. Default 42.
    **thermo_kwargs
        Forwarded to ``thermo_fn`` (e.g. ``solvent="water"``).

    Returns
    -------
    list of Thermo
        The successfully computed conformers (fewer than ``max_conformers``
        if ``max_attempts`` is exhausted first).

    Raises
    ------
    ValueError
        If no conformer converges to a true minimum within ``max_attempts``.
    """
    thermos = []
    for attempt in range(max_attempts):
        atoms_list = generate(
            smiles,
            max_conformers=max_conformers,
            prune_rms_thresh=prune_rms_thresh,
            energy_window=energy_window,
            random_seed=random_seed + attempt,
        )
        for atoms in atoms_list:
            if len(thermos) >= max_conformers:
                break
            try:
                thermos.append(thermo_fn(atoms, charge=charge, **thermo_kwargs))
            except (TSValueError, RuntimeError):
                continue
        if len(thermos) >= max_conformers:
            break

    if not thermos:
        raise ValueError(
            f"No conformer of {smiles!r} converged to a true minimum after "
            f"{max_attempts} attempts."
        )
    return thermos


def write_conformers(conformers, directory, prefix="conformer"):
    """
    Write conformers to ``directory`` as ``<prefix>_<i>.xyz`` files.

    Returns the list of written paths, so the directory can be passed straight to
    :func:`ThermoScreening.thermo.screening.screen`.
    """
    import ase.io
    from pathlib import Path

    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    paths = []
    for index, atoms in enumerate(conformers):
        path = directory / f"{prefix}_{index}.xyz"
        ase.io.write(str(path), atoms)
        paths.append(path)
    return paths
