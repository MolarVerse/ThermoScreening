"""Conformer generation for molecule screening.

RDKit is the conformer backend: ETKDG distance-geometry embedding, optional
MMFF/UFF force-field optimisation, and RMSD + energy-window pruning. Conformers
are returned as ASE ``Atoms`` so they feed straight into the screening pipeline
(e.g. generate -> ``write_conformers`` -> ``screen`` that directory).
"""

import numpy as np


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
    return Atoms(symbols=symbols, positions=positions)


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
