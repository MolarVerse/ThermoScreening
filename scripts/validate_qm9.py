#!/usr/bin/env python3
"""Validate ThermoScreening against a deterministic sample of QM9."""

from __future__ import annotations

import argparse
import hashlib
import os
import shutil
import tarfile
import urllib.request
from pathlib import Path

import numpy as np
from ase import Atoms

from ThermoScreening.thermo.api import run_thermo


FIGSHARE_ARTICLE_ID = "1057646"
ARCHIVE_NAME = "dsgdb9nsd.xyz.tar.bz2"
DOWNLOAD_URL = "https://ndownloader.figshare.com/files/3195389"
ARCHIVE_MD5 = "ad1ebd51ee7f5b3a6e32e974e5d54012"
MOLECULE_COUNT = 133_885
SAMPLE_SIZE = 1_000
RANDOM_SEED = 20_260_724

# QM9 does not store Gaussian's rotational symmetry number. Its published
# thermochemistry identifies these three exact small geometries as sigma=2;
# the remaining records in the pinned sample use sigma=1.
SYMMETRY_NUMBER_OVERRIDES = {3: 2, 4: 2, 23: 2}

LIMITS = {
    "U0": 0.6,
    "U": 1.5,
    "H": 1.5,
    "G": 3.5,
    "Cv": 0.0007,
}


def _md5(path: Path) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _download_archive(cache_dir: Path) -> Path:
    cache_dir.mkdir(parents=True, exist_ok=True)
    destination = cache_dir / ARCHIVE_NAME
    if destination.exists() and _md5(destination) == ARCHIVE_MD5:
        return destination

    temporary = destination.with_suffix(destination.suffix + ".part")
    request = urllib.request.Request(
        DOWNLOAD_URL,
        headers={"User-Agent": "ThermoScreening-reference-validation"},
    )
    try:
        with (
            urllib.request.urlopen(request, timeout=120) as response,
            temporary.open("wb") as output,
        ):
            shutil.copyfileobj(response, output)
        if _md5(temporary) != ARCHIVE_MD5:
            raise RuntimeError(f"Checksum mismatch for {ARCHIVE_NAME}.")
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)
    return destination


def _sample_ids() -> set[int]:
    rng = np.random.default_rng(RANDOM_SEED)
    selected = set(range(1, 101))
    selected.update(np.linspace(101, MOLECULE_COUNT, 450, dtype=int).tolist())
    remaining = np.asarray(
        sorted(set(range(101, MOLECULE_COUNT + 1)) - selected)
    )
    selected.update(rng.choice(remaining, size=450, replace=False).tolist())
    if len(selected) != SAMPLE_SIZE:
        raise RuntimeError(f"Expected {SAMPLE_SIZE} sample identifiers.")
    return selected


def _number(value: str) -> float:
    return float(value.replace("*^", "e"))


def _parse_record(text: str):
    lines = text.splitlines()
    atom_count = int(lines[0])
    properties = lines[1].split()
    atom_rows = [line.split() for line in lines[2 : 2 + atom_count]]
    atoms = Atoms(
        [row[0] for row in atom_rows],
        positions=[
            [_number(value) for value in row[1:4]]
            for row in atom_rows
        ],
    )
    frequencies = np.asarray(
        [_number(value) for value in lines[2 + atom_count].split()]
    )
    expected = {
        "zpve": _number(properties[11]),
        "U0": _number(properties[12]),
        "U": _number(properties[13]),
        "H": _number(properties[14]),
        "G": _number(properties[15]),
        "Cv": _number(properties[16]),
    }
    return int(properties[1]), atoms, frequencies, expected


def _calculate_errors(archive: Path) -> dict[str, np.ndarray]:
    selected = _sample_ids()
    errors = {name: [] for name in LIMITS}
    processed = set()

    with tarfile.open(archive, mode="r:bz2") as tar:
        for member in tar:
            stem = (
                member.name.removeprefix("dsgdb9nsd_").removesuffix(".xyz")
            )
            if not stem.isdigit() or int(stem) not in selected:
                continue

            handle = tar.extractfile(member)
            if handle is None:
                raise RuntimeError(f"Cannot read {member.name}.")
            index, atoms, frequencies, expected = _parse_record(
                handle.read().decode("utf-8")
            )
            if index != int(stem):
                raise RuntimeError(f"Index mismatch in {member.name}.")

            thermo = run_thermo(
                frequencies,
                atoms=atoms,
                energy=expected["U0"] - expected["zpve"],
                temperature=298.15,
                pressure=101325,
                spin=0.0,
                symmetry_number=SYMMETRY_NUMBER_OVERRIDES.get(index, 1),
            )
            actual = {
                "U0": thermo.total_EeZPE(),
                "U": thermo.total_EeEtot(),
                "H": thermo.total_EeHtot(),
                "G": thermo.total_EeGtot(),
                "Cv": thermo.total_heat_capacity("cal/(mol*K)"),
            }
            for name in LIMITS:
                factor = 1.0 if name == "Cv" else 1e6
                errors[name].append((actual[name] - expected[name]) * factor)
            processed.add(index)

    missing = selected - processed
    if missing:
        raise RuntimeError(f"Missing {len(missing)} sampled QM9 records.")
    return {name: np.asarray(values) for name, values in errors.items()}


def _default_cache_dir() -> Path:
    root = Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache"))
    return root / "thermoscreening" / "qm9"


def main() -> int:
    """Download QM9, run the comparison, and report status."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=_default_cache_dir(),
        help="Archive cache directory.",
    )
    args = parser.parse_args()

    archive = _download_archive(args.cache_dir.expanduser().resolve())
    errors = _calculate_errors(archive)

    print(f"Figshare article: {FIGSHARE_ARTICLE_ID}")
    print(f"QM9 molecules sampled: {SAMPLE_SIZE}")
    passed = True
    for name, values in errors.items():
        absolute = np.abs(values)
        unit = "cal/(mol*K)" if name == "Cv" else "microhartree"
        maximum = float(np.max(absolute))
        mean = float(np.mean(absolute))
        print(f"{name}: MAE {mean:.6f}, maximum {maximum:.6f} {unit}")
        passed = passed and maximum <= LIMITS[name]

    print("Result: PASS" if passed else "Result: FAIL")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
