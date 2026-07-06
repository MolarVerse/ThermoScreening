import numpy as np
import pytest
from ase import Atoms

from ThermoScreening.thermo import conformers


def test_generate_returns_ase_conformers():
    # n-butane is flexible enough to embed at least one conformer
    result = conformers.generate("CCCC", max_conformers=5)

    assert 1 <= len(result) <= 5
    assert all(isinstance(atoms, Atoms) for atoms in result)
    assert all(atoms.get_chemical_formula() == "C4H10" for atoms in result)


def test_generate_respects_max_conformers():
    result = conformers.generate("CCCCCCCC", max_conformers=3)
    assert len(result) <= 3


def test_generate_energy_window_filters():
    wide = conformers.generate("CCCCCC", max_conformers=15, energy_window=100.0)
    narrow = conformers.generate("CCCCCC", max_conformers=15, energy_window=0.1)
    assert len(narrow) <= len(wide)


def test_generate_without_optimize_runs():
    result = conformers.generate("CCO", max_conformers=3, optimize=False)
    assert len(result) >= 1


def test_generate_rejects_bad_smiles():
    with pytest.raises(ValueError, match="Could not parse SMILES"):
        conformers.generate("this-is-not-smiles")


def test_write_conformers_writes_readable_xyz(tmp_path):
    import ase.io

    result = conformers.generate("CCO", max_conformers=2)
    paths = conformers.write_conformers(result, tmp_path / "confs", prefix="c")

    assert len(paths) == len(result)
    assert all(p.name.startswith("c_") and p.suffix == ".xyz" for p in paths)
    back = ase.io.read(str(paths[0]))
    assert back.get_chemical_formula() == "C2H6O"


def test_public_api_is_exported():
    from ThermoScreening.thermo import generate_conformers, write_conformers

    assert generate_conformers is conformers.generate
    assert write_conformers is conformers.write_conformers
