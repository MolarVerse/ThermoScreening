"""
Thermochemistry screening over a set of molecules.

``screen`` runs the DFTB+ geometry-optimisation/Hessian/modes pipeline for each
molecule in its own working directory and collects the key thermochemical
quantities into a results table (CSV and JSON). A failing molecule is recorded
with an error status and does not stop the run.
"""

import csv
import json
import logging
from dataclasses import dataclass
from pathlib import Path

import ase.io

from ThermoScreening import __package_name__
from ThermoScreening.calculator.dftbplus import resolve_parameter_set
from ThermoScreening.exceptions import TSValueError
from ThermoScreening.utils.custom_logging import setup_logger

from .api import dftbplus_thermo

logger = logging.getLogger(__package_name__).getChild("screening")
logger = setup_logger(logger)

_STRUCTURE_SUFFIXES = (".xyz", ".gen")
_RESULT_FIELDS = [
    "name",
    "path",
    "charge",
    "status",
    "Eelec_hartree",
    "E_hartree",
    "H_hartree",
    "G_hartree",
    "G_total_hartree",
    "S_cal_per_mol_K",
    "Cv_cal_per_mol_K",
    "error",
]


@dataclass
class ScreeningJob:
    """A single molecule to screen."""

    name: str
    path: Path
    charge: float
    spin: float | None = None


def _jobs_from_directory(directory: Path, charge: float, spin=None):
    jobs = [
        ScreeningJob(name=path.stem, path=path, charge=charge, spin=spin)
        for path in sorted(directory.iterdir())
        if path.suffix.lower() in _STRUCTURE_SUFFIXES
    ]
    if not jobs:
        raise TSValueError(f"No .xyz or .gen structures found in '{directory}'.")
    return jobs


def _jobs_from_manifest(manifest: Path, charge: float, spin=None):
    base = manifest.parent
    jobs = []
    with open(manifest, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if "path" not in row or not row["path"]:
                raise TSValueError(f"Manifest '{manifest}' needs a 'path' column.")
            path = Path(row["path"])
            if not path.is_absolute():
                path = base / path
            row_charge = row.get("charge")
            row_spin = row.get("spin")
            jobs.append(
                ScreeningJob(
                    name=row.get("name") or path.stem,
                    path=path,
                    charge=float(row_charge) if row_charge not in (None, "") else charge,
                    spin=float(row_spin) if row_spin not in (None, "") else spin,
                )
            )
    if not jobs:
        raise TSValueError(f"Manifest '{manifest}' lists no molecules.")
    return jobs


def _load_jobs(source, charge: float, spin=None):
    source = Path(source)
    if source.is_dir():
        return _jobs_from_directory(source, charge, spin)
    if source.suffix.lower() == ".csv":
        return _jobs_from_manifest(source, charge, spin)
    raise TSValueError(
        f"Screening input must be a directory or a .csv manifest, got '{source}'."
    )


def _thermo_summary(thermo):
    # Eelec_hartree is the DFTB+ electronic (SCC) energy, which carries the
    # implicit-solvation term; E/H/G_hartree are the thermal corrections and
    # G_total_hartree = Eelec + G correction is the absolute Gibbs free energy.
    return {
        "Eelec_hartree": thermo.electronic_energy(),
        "E_hartree": thermo.total_energy("H"),
        "H_hartree": thermo.total_enthalpy("H"),
        "G_hartree": thermo.total_gibbs_free_energy("H"),
        "G_total_hartree": thermo.total_EeGtot(),
        "S_cal_per_mol_K": thermo.total_entropy("cal/(mol*K)"),
        "Cv_cal_per_mol_K": thermo.total_heat_capacity("cal/(mol*K)"),
    }


def _write_results(results, out):
    stem = Path(str(out)).with_suffix("")
    if stem.parent != Path(""):
        stem.parent.mkdir(parents=True, exist_ok=True)

    csv_path = stem.with_suffix(".csv")
    with open(csv_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=_RESULT_FIELDS)
        writer.writeheader()
        for record in results:
            writer.writerow({field: record.get(field, "") for field in _RESULT_FIELDS})

    json_path = stem.with_suffix(".json")
    json_path.write_text(json.dumps(results, indent=2), encoding="utf-8")

    return csv_path, json_path


def screen(
    source,
    out="results",
    charge=0.0,
    temperature=298.15,
    pressure=101325,
    directory="screening",
    parameters=None,
    spin=None,
    parameter_set="3ob",
    solvent=None,
):
    """
    Run a thermochemistry screen over a set of molecules.

    Parameters
    ----------
    source : str
        A directory of ``.xyz``/``.gen`` structures (all run at ``charge``), or
        a ``.csv`` manifest with ``path`` and optional ``name``/``charge``/``spin``
        columns.
    out : str
        Output stem; results are written to ``<out>.csv`` and ``<out>.json``.
    charge : float
        Charge used for directory input and for manifest rows without a charge.
    temperature : float
        Temperature in K.
    pressure : float
        Pressure in Pa.
    directory : str
        Root working directory; each molecule runs in ``<directory>/<name>``.
    parameters : dict, optional
        DFTB+ Hamiltonian parameters. Defaults to those of ``parameter_set``.
    spin : float, optional
        Spin quantum number S applied to molecules without a manifest ``spin``;
        defaults to an electron-count guess per molecule.
    parameter_set : str
        Slater-Koster parameter set to use (``"3ob"`` or ``"mio"``). Selects both
        the default Hamiltonian parameters and the matching spin constants.
    solvent : str, optional
        Solvent name for GBSA/ALPB implicit solvation applied to every molecule
        (e.g. ``"water"``). Defaults to gas phase. The solvent parameter file
        must be installed (``thermo setup-dftb --solvent <name>``).

    Returns
    -------
    list of dict
        One result record per molecule, including failed ones (status="error").
    """
    default_parameters, spin_constants = resolve_parameter_set(parameter_set)
    parameters = default_parameters if parameters is None else parameters
    jobs = _load_jobs(source, charge, spin)
    root = Path(directory)

    results = []
    for job in jobs:
        record = {
            "name": job.name,
            "path": str(job.path),
            "charge": job.charge,
            "status": "ok",
            "error": "",
        }
        logger.info(f"Screening {job.name} (charge {job.charge})")
        try:
            atoms = ase.io.read(str(job.path))
            thermo = dftbplus_thermo(
                atoms,
                temperature=temperature,
                pressure=pressure,
                charge=job.charge,
                directory=str(root / job.name),
                spin=job.spin,
                spin_constants=spin_constants,
                solvent=solvent,
                **parameters,
            )
            record.update(_thermo_summary(thermo))
        except Exception as exc:  # pylint: disable=broad-except
            # isolate failures so one bad molecule does not abort the screen
            record["status"] = "error"
            record["error"] = str(exc)
            # warning, not error: the active custom logger raises on error and
            # we must keep screening the remaining molecules
            logger.warning(f"Screening failed for {job.name}: {exc}")
        results.append(record)

    csv_path, json_path = _write_results(results, out)
    logger.info(f"Wrote {csv_path} and {json_path}")
    return results
