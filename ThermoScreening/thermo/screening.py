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
import os
from dataclasses import dataclass
from pathlib import Path

import ase.io

from ThermoScreening import __package_name__
from ThermoScreening.calculator.dftbplus import resolve_parameter_set
from ThermoScreening.exceptions import TSValueError
from ThermoScreening.utils.custom_logging import setup_logger

from .api import dftbplus_thermo, xtb_thermo, xtb_cli_thermo

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


def rank_by_gibbs(results):
    """
    Return the successful screen records sorted by Gibbs free energy.

    Parameters
    ----------
    results : list of dict
        Records as returned by :func:`screen`.

    Returns
    -------
    list of dict
        The records whose ``status`` is ``"ok"`` (and that carry a
        ``G_total_hartree``), sorted by absolute Gibbs free energy ascending --
        i.e. the most stable molecule first. Failed records are omitted.
    """
    ranked = [
        record
        for record in results
        if record.get("status") == "ok" and record.get("G_total_hartree") is not None
    ]
    return sorted(ranked, key=lambda record: record["G_total_hartree"])


def _load_completed(out):
    """
    Load the successfully-completed records from a prior ``<out>.json``.

    Returns a ``{name: record}`` mapping of the records whose ``status`` is
    ``"ok"``, so a resumed screen can skip them. A missing or unreadable file
    yields an empty mapping (nothing to resume).
    """
    json_path = Path(str(out)).with_suffix(".json")
    if not json_path.exists():
        return {}
    try:
        prior = json.loads(json_path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return {}
    if not isinstance(prior, list):
        return {}
    return {
        record["name"]: record
        for record in prior
        if isinstance(record, dict) and record.get("status") == "ok" and "name" in record
    }


def _atomic_write(path, write_fn):
    """
    Write via a temp file + ``os.replace`` so an interrupted write never leaves a
    truncated (unresumable) file behind.
    """
    tmp = path.with_name(path.name + ".tmp")
    with open(tmp, "w", newline="", encoding="utf-8") as handle:
        write_fn(handle)
    os.replace(tmp, path)


def _write_results(results, out):
    stem = Path(str(out)).with_suffix("")
    if stem.parent != Path(""):
        stem.parent.mkdir(parents=True, exist_ok=True)

    def write_csv(handle):
        writer = csv.DictWriter(handle, fieldnames=_RESULT_FIELDS)
        writer.writeheader()
        for record in results:
            writer.writerow({field: record.get(field, "") for field in _RESULT_FIELDS})

    csv_path = stem.with_suffix(".csv")
    json_path = stem.with_suffix(".json")
    _atomic_write(csv_path, write_csv)
    _atomic_write(json_path, lambda handle: handle.write(json.dumps(results, indent=2)))

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
    dispersion=None,
    quasi_rrho=False,
    engine="dftb+",
    method="GFN2-xTB",
    resume=False,
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
        must be installed (``thermo setup-dftb --solvent <name>``). With the
        DFTB+ engine the GBSA parameters are GFN-xTB-fit, so solvation free
        energies are only qualitative (and can be non-monotonic in the dielectric
        for small neutral solutes).
    dispersion : str, optional
        Dispersion correction for the DFTB+ engine. ``"d3-bj"`` adds Grimme
        D3(BJ) with the 3ob-recommended parameters. Defaults to none.
    quasi_rrho : bool
        If True, use Grimme's quasi-RRHO treatment for the vibrational entropy
        (recommended for flexible molecules with low-frequency modes). Default
        False (pure harmonic oscillator).
    engine : str
        Calculation engine: ``"dftb+"`` (default), ``"xtb"`` (GFN-xTB via tblite,
        gas-phase) or ``"xtb-cli"`` (native ``xtb`` binary, which additionally
        supports implicit solvation of charged radicals). ``parameter_set``
        applies only to DFTB+; ``solvent`` applies to DFTB+ and xtb-cli.
    method : str
        GFN-xTB parametrisation for the xtb engines (``"GFN2-xTB"`` default).
    resume : bool
        If True, reuse the successful (``status="ok"``) records from a prior
        ``<out>.json`` and skip those molecules; molecules that previously failed
        (or were never run) are (re-)run. Results are written incrementally after
        every molecule regardless, so an interrupted screen can be resumed.

    Returns
    -------
    list of dict
        One result record per molecule, including failed ones (status="error").
    """
    if engine not in ("dftb+", "xtb", "xtb-cli"):
        raise TSValueError(
            f"Unknown engine {engine!r}; choose 'dftb+', 'xtb' or 'xtb-cli'."
        )

    if engine == "dftb+":
        default_parameters, spin_constants = resolve_parameter_set(parameter_set)
        parameters = default_parameters if parameters is None else parameters

    jobs = _load_jobs(source, charge, spin)
    root = Path(directory)

    completed = _load_completed(out) if resume else {}

    results = []
    csv_path = json_path = None
    for job in jobs:
        if job.name in completed:
            logger.info(f"Skipping {job.name} (already completed)")
            results.append(completed[job.name])
            csv_path, json_path = _write_results(results, out)
            continue
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
            if engine == "xtb-cli":
                thermo = xtb_cli_thermo(
                    atoms,
                    temperature=temperature,
                    pressure=pressure,
                    charge=job.charge,
                    directory=str(root / job.name),
                    spin=job.spin,
                    method=method,
                    solvent=solvent,
                    quasi_rrho=quasi_rrho,
                )
            elif engine == "xtb":
                thermo = xtb_thermo(
                    atoms,
                    temperature=temperature,
                    pressure=pressure,
                    charge=job.charge,
                    directory=str(root / job.name),
                    spin=job.spin,
                    method=method,
                    quasi_rrho=quasi_rrho,
                )
            else:
                thermo = dftbplus_thermo(
                    atoms,
                    temperature=temperature,
                    pressure=pressure,
                    charge=job.charge,
                    directory=str(root / job.name),
                    spin=job.spin,
                    spin_constants=spin_constants,
                    solvent=solvent,
                    dispersion=dispersion,
                    quasi_rrho=quasi_rrho,
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
        # write after every molecule so an interrupted run stays resumable
        csv_path, json_path = _write_results(results, out)

    logger.info(f"Wrote {csv_path} and {json_path}")
    return results
