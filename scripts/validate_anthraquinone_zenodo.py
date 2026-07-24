#!/usr/bin/env python3
"""Reproduce the published DFTB/3ob/COSMO anthraquinone potentials."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import shutil
import tarfile
import urllib.request
from pathlib import Path

import numpy as np

from ThermoScreening.thermo.api import read_vibrational, run_thermo
from ThermoScreening.thermo.reactions import (
    calibrate_reduction_reference,
    reduction_potential,
)


RECORD_ID = "20796838"
BASE_URL = f"https://zenodo.org/records/{RECORD_ID}/files"
ARCHIVES = {
    "04_dftb_3ob_cosmo_raw_outputs.tar.xz": (
        "c0ef53aa96547e10cf8478192abcd0dc0d09a0464ff41877016aa9fe461436dd"
    ),
    "08_processed_tables.tar.xz": (
        "4a2a083fd1c702e882844bfca210848b7dcf8afbec09ffb87e4e31db4121ab47"
    ),
}
MAX_ERROR_MV = 0.5
MAX_MEAN_ERROR_MV = 0.05

LABEL_TO_ORDER = {
    "AQ": 1,
    "1,2-OH": 2,
    "1,4-OH": 3,
    "1,5-OH": 4,
    "1,8-OH": 5,
    "1,2-NH2": 6,
    "1,4-NH2": 7,
    "2,6-NH2": 8,
    "1-OH": 9,
    "2-OH": 10,
    "1-NH2": 11,
    "2-NH2": 12,
    "1-NH-4-OH": 13,
    "1-NH2-4-OH": 13,
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _download_archive(cache_dir: Path, name: str, checksum: str) -> Path:
    destination = cache_dir / name
    if destination.exists() and _sha256(destination) == checksum:
        return destination

    temporary = destination.with_suffix(destination.suffix + ".part")
    request = urllib.request.Request(
        f"{BASE_URL}/{name}?download=1",
        headers={"User-Agent": "ThermoScreening-reference-validation"},
    )
    try:
        with (
            urllib.request.urlopen(request, timeout=60) as response,
            temporary.open("wb") as output,
        ):
            shutil.copyfileobj(response, output)
        if _sha256(temporary) != checksum:
            raise RuntimeError(f"Checksum mismatch for {name}.")
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)
    return destination


def _prepare_data(cache_dir: Path) -> Path:
    cache_dir.mkdir(parents=True, exist_ok=True)
    for name, checksum in ARCHIVES.items():
        archive = _download_archive(cache_dir, name, checksum)
        with tarfile.open(archive, mode="r:xz") as tar:
            tar.extractall(cache_dir, filter="data")
    return cache_dir


def _read_published_potentials(path: Path) -> dict[int, float]:
    rows = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter=";"):
            if not row or row[0].startswith("#"):
                continue
            try:
                order = LABEL_TO_ORDER[row[0]]
            except KeyError as exc:
                raise RuntimeError(
                    f"Unknown molecule label {row[0]!r} in {path}."
                ) from exc
            if order in rows:
                raise RuntimeError(f"Duplicate molecule {order} in {path}.")
            rows[order] = float(row[1]) / 1000.0
    return rows


def _calculate_states(data_dir: Path):
    calculations = (
        data_dir / "dftb_3ob_cosmo_dmf" / "calculations"
    )
    states = {}
    for order in range(1, 14):
        matches = list(calculations.glob(f"{order:02d}_*"))
        if len(matches) != 1:
            raise RuntimeError(
                f"Expected one calculation directory for molecule {order}, "
                f"found {len(matches)}."
            )
        for charge in (0, -1, -2):
            thermo_dir = matches[0] / f"run_{charge}" / "thermo"
            frequencies = read_vibrational(
                str(thermo_dir / "frequency.txt"), "dftb+"
            )
            energy = float(
                (thermo_dir / "electronic_energy.txt").read_text(encoding="utf-8")
            )
            states[order, charge] = run_thermo(
                frequencies,
                coord_file=str(thermo_dir / "geo_opt.xyz"),
                energy=energy,
                charge=charge,
                spin=1.0 if charge == -1 else 0.0,
                symmetry_number=1,
            )
    return states


def _compare(data_dir: Path) -> tuple[np.ndarray, tuple[str, int, float]]:
    table_dir = (
        data_dir
        / "processed_results"
        / "tables"
        / "dft_and_dftb_paper_plot_tables"
        / "dftb_legacy_tables_from_plot_workflow"
    )
    published_e1 = _read_published_potentials(
        table_dir / "calc_exp_1_DFTB_COSMO_3ob.csv"
    )
    published_e2e = _read_published_potentials(
        table_dir / "calc_exp_2_DFTB_COSMO_3ob.csv"
    )
    if set(published_e1) != set(range(1, 14)):
        raise RuntimeError("The first-reduction table does not contain 13 molecules.")
    if set(published_e2e) != {*range(1, 12), 13}:
        raise RuntimeError("The two-electron table does not contain 12 molecules.")
    states = _calculate_states(data_dir)

    e1_reference = calibrate_reduction_reference(
        states[1, 0],
        states[1, -1],
        experimental_potential=published_e1[1],
    )
    e2e_reference = calibrate_reduction_reference(
        states[1, 0],
        states[1, -2],
        experimental_potential=published_e2e[1],
        n_electrons=2,
    )

    errors = []
    labels = []
    for order, published in published_e1.items():
        calculated = reduction_potential(
            states[order, 0],
            states[order, -1],
            reference_potential=e1_reference,
        )
        errors.append(calculated - published)
        labels.append(("E1", order))
    for order, published in published_e2e.items():
        calculated = reduction_potential(
            states[order, 0],
            states[order, -2],
            n_electrons=2,
            reference_potential=e2e_reference,
        )
        errors.append(calculated - published)
        labels.append(("E2e", order))

    errors_mv = np.asarray(errors) * 1000.0
    worst_index = int(np.argmax(np.abs(errors_mv)))
    worst = (*labels[worst_index], float(errors_mv[worst_index]))
    return errors_mv, worst


def _default_cache_dir() -> Path:
    root = Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache"))
    return root / "thermoscreening" / f"anthraquinone-{RECORD_ID}"


def main() -> int:
    """Download the reference data, run the comparison, and report status."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=_default_cache_dir(),
        help="Archive and extraction cache directory.",
    )
    args = parser.parse_args()

    data_dir = _prepare_data(args.cache_dir.expanduser().resolve())
    errors_mv, worst = _compare(data_dir)
    mean_error = float(np.mean(np.abs(errors_mv)))
    max_error = float(np.max(np.abs(errors_mv)))

    print(f"Zenodo record: {RECORD_ID}")
    print("Charge-state calculations: 39")
    print(f"Published potentials: {len(errors_mv)}")
    print(f"Mean absolute error: {mean_error:.6f} mV")
    print(f"Maximum absolute error: {max_error:.6f} mV")
    print(
        f"Worst case: {worst[0]}, molecule {worst[1]}, "
        f"signed error {worst[2]:+.6f} mV"
    )

    passed = max_error <= MAX_ERROR_MV and mean_error <= MAX_MEAN_ERROR_MV
    print("Result: PASS" if passed else "Result: FAIL")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
