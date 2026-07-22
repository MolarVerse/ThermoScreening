"""
Thermochemistry screening over a set of molecules.

``screen`` runs the DFTB+ geometry-optimisation/Hessian/modes pipeline for each
molecule in its own working directory and collects the key thermochemical
quantities into a results table (CSV and JSON). A failing molecule is recorded
with an error status and does not stop the run.
"""

import csv
import functools
import hashlib
import json
import logging
import math
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

import ase.io

from ThermoScreening import __package_name__
from ThermoScreening.calculator.dftbplus import resolve_parameter_set
from ThermoScreening.exceptions import TSValueError
from ThermoScreening.utils.custom_logging import setup_logger
from ThermoScreening.version import __version__

from .api import dftbplus_thermo, xtb_thermo, xtb_cli_thermo
from ._units import HARTREE_TO_EV
from .conformers import generate as generate_conformers, write_conformers
from .reactions import SHE_ABSOLUTE_POTENTIAL

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
    "fingerprint",
    "error",
]
_REDOX_RESULT_FIELDS = [
    "name",
    "path",
    "charge",
    "status",
    "G_oxidized_hartree",
    "G_reduced_once_hartree",
    "G_reduced_twice_hartree",
    "E1_V",
    "E2_V",
    "E2e_V",
    "potential_gap_V",
    "potential_inversion",
    "potential_scale",
    "run_fingerprint",
    "error",
]
_REDOX_STATES = (
    ("oxidized", 0, 0),
    ("reduced_once", -1, 1),
    ("reduced_twice", -2, 2),
)


@dataclass
class ScreeningJob:
    """A single molecule to screen."""

    name: str
    path: Path
    charge: float
    spin: float | None = None


@dataclass
class _RedoxJob:
    """One molecule and the three electronic states in a redox workflow."""

    name: str
    path: Path
    charge: float
    spins: tuple[float | None, float | None, float | None]


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
        jobs = _jobs_from_directory(source, charge, spin)
    elif source.suffix.lower() == ".csv":
        jobs = _jobs_from_manifest(source, charge, spin)
    else:
        raise TSValueError(
            f"Screening input must be a directory or a .csv manifest, got '{source}'."
        )

    # Names key the result records and the per-molecule working directory, so a
    # collision would silently drop a result and (in parallel) race two jobs into
    # the same directory. Reject duplicates up front. (E.g. a directory holding
    # both mol.xyz and mol.gen, or a manifest listing a name twice.)
    names = [job.name for job in jobs]
    duplicates = sorted({name for name in names if names.count(name) > 1})
    if duplicates:
        raise TSValueError(
            f"Duplicate molecule name(s) in the screen: {', '.join(duplicates)}. "
            "Each molecule needs a unique name."
        )
    return jobs


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
        ``G_total_hartree``), sorted by total Gibbs free energy ascending --
        i.e. the most stable molecule first. Failed records are omitted.
    """
    ranked = [
        record
        for record in results
        if record.get("status") == "ok" and record.get("G_total_hartree") is not None
    ]
    return sorted(ranked, key=lambda record: record["G_total_hartree"])


def _stable_value(value):
    """Return a deterministic JSON-compatible representation."""
    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _stable_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_stable_value(item) for item in value]
    if hasattr(value, "tolist"):
        return _stable_value(value.tolist())
    return repr(value)


def _file_sha256(path):
    """Hash a structure without making a missing input abort the whole screen."""
    digest = hashlib.sha256()
    try:
        with open(path, "rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
    except OSError:
        return None
    return digest.hexdigest()


def _payload_fingerprint(payload):
    encoded = json.dumps(
        _stable_value(payload), sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _screen_settings(
    *,
    engine,
    temperature,
    pressure,
    method,
    solvent,
    dispersion,
    quasi_rrho,
    parameter_set,
    parameters,
    spin_constants,
):
    return {
        "thermoscreening_version": __version__,
        "engine": engine,
        "temperature_K": temperature,
        "pressure_Pa": pressure,
        "method": method if engine in ("xtb", "xtb-cli") else None,
        "solvent": solvent,
        "dispersion": dispersion if engine == "dftb+" else None,
        "quasi_rrho": quasi_rrho,
        "parameter_set": parameter_set if engine == "dftb+" else None,
        "parameters": parameters if engine == "dftb+" else None,
        "spin_constants": spin_constants if engine == "dftb+" else None,
    }


def _job_provenance(job, settings, input_index=None):
    structure_hash = _file_sha256(job.path)
    payload = {
        "name": job.name,
        "path": str(job.path),
        "structure_sha256": structure_hash,
        "charge": job.charge,
        "spin": job.spin,
        "settings": settings,
    }
    provenance = {**payload, "fingerprint": _payload_fingerprint(payload)}
    if input_index is not None:
        provenance["input_index"] = input_index
    return provenance


def _load_completed(out, expected_fingerprints=None):
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
    completed = {}
    for record in prior:
        if not (
            isinstance(record, dict)
            and record.get("status") == "ok"
            and "name" in record
        ):
            continue
        if expected_fingerprints is not None:
            expected = expected_fingerprints.get(record["name"])
            if expected is None or record.get("fingerprint") != expected:
                continue
        completed[record["name"]] = record
    return completed


def _atomic_write(path, write_fn):
    """
    Write via a temp file + ``os.replace`` so an interrupted write never leaves a
    truncated (unresumable) file behind.
    """
    tmp = path.with_name(path.name + ".tmp")
    with open(tmp, "w", newline="", encoding="utf-8") as handle:
        write_fn(handle)
    os.replace(tmp, path)


def _sidecar_path(out, suffix):
    stem = Path(str(out)).with_suffix("")
    return stem.parent / f"{stem.name}{suffix}"


def _write_json(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    _atomic_write(
        path,
        lambda handle: handle.write(
            json.dumps(_stable_value(payload), indent=2, sort_keys=True, allow_nan=False)
            + "\n"
        ),
    )


def _write_screen_run_metadata(out, source, settings, job_provenance, shard=None):
    payload = {
        "schema_version": 1,
        "workflow": "thermochemistry_screen",
        "source": str(source),
        "settings": settings,
        "jobs": job_provenance,
    }
    if shard is not None:
        payload["shard"] = shard
    path = _sidecar_path(out, "-run.json")
    _write_json(path, payload)
    return path


def _validate_shard(shard_index, shard_count):
    if shard_index is None and shard_count is None:
        return None
    if shard_index is None or shard_count is None:
        raise TSValueError("shard_index and shard_count must be supplied together.")
    if isinstance(shard_index, bool) or not isinstance(shard_index, int):
        raise TSValueError("shard_index must be an integer.")
    if isinstance(shard_count, bool) or not isinstance(shard_count, int):
        raise TSValueError("shard_count must be an integer.")
    if shard_count < 1:
        raise TSValueError("shard_count must be >= 1.")
    if not 0 <= shard_index < shard_count:
        raise TSValueError(
            f"shard_index must satisfy 0 <= index < {shard_count}, got {shard_index}."
        )
    return {"index": shard_index, "count": shard_count}


def screen_shard_directory(out):
    """Return the directory containing distributed screen result shards."""
    stem = Path(str(out)).with_suffix("")
    return stem.parent / f"{stem.name}-shards"


def _screen_shard_stem(out, shard_index):
    return screen_shard_directory(out) / f"shard-{shard_index:05d}"


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


def _run_job(
    job,
    *,
    fingerprint,
    engine,
    temperature,
    pressure,
    root,
    method,
    solvent,
    dispersion,
    quasi_rrho,
    parameters,
    spin_constants,
):
    """
    Run one screening job and return its result record.

    A module-level, picklable worker (usable from a process pool). Failures are
    isolated into the record (``status="error"``) rather than raised, so one bad
    molecule never aborts the screen; the caller does the logging.
    """
    record = {
        "name": job.name,
        "path": str(job.path),
        "charge": job.charge,
        "status": "ok",
        "fingerprint": fingerprint,
        "error": "",
    }
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
        record["status"] = "error"
        record["error"] = str(exc)
    return record


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
    jobs=1,
    shard_index=None,
    shard_count=None,
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
    jobs : int
        Number of molecules to run concurrently. Default 1 (serial). With
        ``jobs > 1`` the molecules run in a process pool (process-based because
        each job changes the working directory). Results and output ordering are
        unchanged; only the wall-clock time differs.
    shard_index, shard_count : int, optional
        Select one deterministic zero-based shard of the input for execution on
        a cluster. Both values are required together. Shards write to isolated
        output and calculation directories and can be combined with
        :func:`collect_screen_shards`.

    Returns
    -------
    list of dict
        One result record per molecule, including failed ones (status="error").
    """
    if engine not in ("dftb+", "xtb", "xtb-cli"):
        raise TSValueError(
            f"Unknown engine {engine!r}; choose 'dftb+', 'xtb' or 'xtb-cli'."
        )
    if jobs < 1:
        raise TSValueError(f"jobs must be >= 1, got {jobs}.")
    if not math.isfinite(float(temperature)) or float(temperature) <= 0:
        raise TSValueError("temperature must be finite and positive.")
    if not math.isfinite(float(pressure)) or float(pressure) <= 0:
        raise TSValueError("pressure must be finite and positive.")
    shard = _validate_shard(shard_index, shard_count)

    spin_constants = None
    if engine == "dftb+":
        default_parameters, spin_constants = resolve_parameter_set(parameter_set)
        parameters = default_parameters if parameters is None else parameters

    all_jobs = _load_jobs(source, charge, spin)
    root = Path(directory)
    settings = _screen_settings(
        engine=engine,
        temperature=temperature,
        pressure=pressure,
        method=method,
        solvent=solvent,
        dispersion=dispersion,
        quasi_rrho=quasi_rrho,
        parameter_set=parameter_set,
        parameters=parameters,
        spin_constants=spin_constants,
    )
    if shard is None:
        selected = list(range(len(all_jobs)))
        job_list = all_jobs
    else:
        selected = [
            index
            for index in range(len(all_jobs))
            if index % shard["count"] == shard["index"]
        ]
        job_list = [all_jobs[index] for index in selected]
        out = _screen_shard_stem(out, shard["index"])
        root = root / "shards" / f"shard-{shard['index']:05d}"
    provenance = [
        _job_provenance(all_jobs[index], settings, input_index=index)
        for index in selected
    ]
    fingerprints = {item["name"]: item["fingerprint"] for item in provenance}
    _write_screen_run_metadata(out, source, settings, provenance, shard=shard)
    completed = (
        _load_completed(out, expected_fingerprints=fingerprints) if resume else {}
    )

    run_one = functools.partial(
        _run_job,
        engine=engine,
        temperature=temperature,
        pressure=pressure,
        root=root,
        method=method,
        solvent=solvent,
        dispersion=dispersion,
        quasi_rrho=quasi_rrho,
        parameters=parameters,
        spin_constants=spin_constants,
    )

    # records keyed by name; the returned list follows the input job order
    records = {}
    csv_path = json_path = None

    def _ordered():
        return [records[job.name] for job in job_list if job.name in records]

    def _store(job, record):
        nonlocal csv_path, json_path
        if record["status"] != "ok":
            # warning, not error: the active custom logger raises on error and
            # we must keep screening the remaining molecules
            logger.warning(f"Screening failed for {job.name}: {record['error']}")
        records[job.name] = record
        # write after every molecule so an interrupted run stays resumable
        csv_path, json_path = _write_results(_ordered(), out)

    pending = []
    for job in job_list:
        if job.name in completed:
            logger.info(f"Skipping {job.name} (already completed)")
            records[job.name] = completed[job.name]
        else:
            pending.append(job)
    if records:
        csv_path, json_path = _write_results(_ordered(), out)
    elif not pending:
        csv_path, json_path = _write_results([], out)

    if jobs > 1 and len(pending) > 1:
        with ProcessPoolExecutor(max_workers=min(jobs, len(pending))) as executor:
            future_to_job = {}
            for job in pending:
                logger.info(f"Screening {job.name} (charge {job.charge})")
                future_to_job[
                    executor.submit(
                        run_one, job, fingerprint=fingerprints[job.name]
                    )
                ] = job
            for future in as_completed(future_to_job):
                job = future_to_job[future]
                try:
                    record = future.result()
                except Exception as exc:  # pylint: disable=broad-except
                    # the worker process itself died (BrokenProcessPool, OOM, a
                    # segfault in a C extension); record this job as failed and
                    # keep going instead of aborting the whole screen
                    record = {
                        "name": job.name,
                        "path": str(job.path),
                        "charge": job.charge,
                        "status": "error",
                        "fingerprint": fingerprints[job.name],
                        "error": str(exc),
                    }
                _store(job, record)
    else:
        for job in pending:
            logger.info(f"Screening {job.name} (charge {job.charge})")
            _store(job, run_one(job, fingerprint=fingerprints[job.name]))

    logger.info(f"Wrote {csv_path} and {json_path}")
    return _ordered()


def _read_json(path, expected_type):
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError) as exc:
        raise TSValueError(f"Could not read cluster result '{path}': {exc}") from exc
    if not isinstance(payload, expected_type):
        raise TSValueError(f"Cluster result '{path}' has an invalid JSON structure.")
    return payload


def collect_screen_shards(shard_directory, out="results"):
    """Validate and combine all shards from a distributed screen."""
    shard_directory = Path(shard_directory)
    metadata_paths = sorted(shard_directory.glob("shard-*-run.json"))
    if not metadata_paths:
        raise TSValueError(f"No screen shards found in '{shard_directory}'.")

    expected_count = None
    expected_settings = None
    expected_source = None
    seen_shards = set()
    records_by_name = {}
    provenance_by_name = {}

    for metadata_path in metadata_paths:
        metadata = _read_json(metadata_path, dict)
        if metadata.get("workflow") != "thermochemistry_screen":
            raise TSValueError(f"'{metadata_path}' is not screen run metadata.")
        shard = metadata.get("shard")
        if not isinstance(shard, dict) or not {"index", "count"} <= set(shard):
            raise TSValueError(f"'{metadata_path}' has no valid shard metadata.")
        shard_index = shard["index"]
        shard_count = shard["count"]
        _validate_shard(shard_index, shard_count)
        if shard_index in seen_shards:
            raise TSValueError(f"Duplicate cluster shard index {shard_index}.")
        seen_shards.add(shard_index)

        settings = metadata.get("settings")
        source = metadata.get("source")
        if expected_count is None:
            expected_count = shard_count
            expected_settings = settings
            expected_source = source
        elif shard_count != expected_count:
            raise TSValueError("Cluster shards disagree on shard_count.")
        elif settings != expected_settings or source != expected_source:
            raise TSValueError("Cluster shards were produced from different inputs or settings.")

        result_name = metadata_path.name.removesuffix("-run.json") + ".json"
        result_path = metadata_path.with_name(result_name)
        if not result_path.is_file():
            raise TSValueError(f"Result file is missing for cluster shard {shard_index}.")
        shard_records = _read_json(result_path, list)
        shard_jobs = metadata.get("jobs")
        if not isinstance(shard_jobs, list):
            raise TSValueError(f"'{metadata_path}' has invalid job metadata.")
        jobs_by_name = {
            job.get("name"): job for job in shard_jobs if isinstance(job, dict)
        }
        if (
            len(jobs_by_name) != len(shard_jobs)
            or any(not isinstance(name, str) or not name for name in jobs_by_name)
        ):
            raise TSValueError(f"'{metadata_path}' contains invalid or duplicate jobs.")
        shard_records_by_name = {
            record.get("name"): record
            for record in shard_records
            if isinstance(record, dict)
        }
        if (
            len(shard_records_by_name) != len(shard_records)
            or set(shard_records_by_name) != set(jobs_by_name)
        ):
            raise TSValueError(
                f"Cluster shard {shard_index} results do not match its job metadata."
            )

        for name, job in jobs_by_name.items():
            record = shard_records_by_name[name]
            provenance_payload = {
                key: job.get(key)
                for key in (
                    "name",
                    "path",
                    "structure_sha256",
                    "charge",
                    "spin",
                    "settings",
                )
            }
            if (
                job.get("settings") != settings
                or _payload_fingerprint(provenance_payload) != job.get("fingerprint")
            ):
                raise TSValueError(
                    f"Invalid provenance fingerprint for {name!r} in "
                    f"cluster shard {shard_index}."
                )
            if record.get("fingerprint") != job.get("fingerprint"):
                raise TSValueError(
                    f"Fingerprint mismatch for {name!r} in cluster shard {shard_index}."
                )
            if name in records_by_name:
                raise TSValueError(f"Duplicate molecule {name!r} across cluster shards.")
            if not isinstance(job.get("input_index"), int):
                raise TSValueError(f"Missing input index for {name!r}.")
            if job["input_index"] % shard_count != shard_index:
                raise TSValueError(
                    f"Input position for {name!r} does not belong to "
                    f"cluster shard {shard_index}."
                )
            records_by_name[name] = record
            provenance_by_name[name] = job

    missing = sorted(set(range(expected_count)) - seen_shards)
    if missing:
        formatted = ", ".join(str(index) for index in missing)
        raise TSValueError(f"Missing cluster shard(s): {formatted}.")

    provenance = sorted(provenance_by_name.values(), key=lambda job: job["input_index"])
    input_indices = [job["input_index"] for job in provenance]
    if len(input_indices) != len(set(input_indices)):
        raise TSValueError("Duplicate input positions across cluster shards.")
    if input_indices != list(range(len(input_indices))):
        raise TSValueError("Cluster shards do not cover every input position.")
    results = [records_by_name[job["name"]] for job in provenance]
    csv_path, json_path = _write_results(results, out)
    _write_json(
        _sidecar_path(out, "-run.json"),
        {
            "schema_version": 1,
            "workflow": "thermochemistry_screen",
            "source": expected_source,
            "settings": expected_settings,
            "jobs": provenance,
            "collection": {
                "shard_count": expected_count,
                "shard_directory": str(shard_directory),
            },
        },
    )
    logger.info(f"Collected {len(results)} molecules into {csv_path} and {json_path}")
    return results


def _optional_float(value, default, field, row_number=None):
    if value in (None, ""):
        return default
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        location = f" in row {row_number}" if row_number is not None else ""
        raise TSValueError(f"{field}{location} must be numeric, got {value!r}.") from exc
    if not math.isfinite(number):
        location = f" in row {row_number}" if row_number is not None else ""
        raise TSValueError(f"{field}{location} must be finite.")
    return number


def _integer_charge(value, default, field, row_number=None):
    number = _optional_float(value, default, field, row_number)
    if number is None:
        return None
    if not float(number).is_integer():
        location = f" in row {row_number}" if row_number is not None else ""
        raise TSValueError(f"{field}{location} must be an integer.")
    return int(number)


def _valid_spin(value, default, field, row_number=None):
    spin = _optional_float(value, default, field, row_number)
    if spin is None:
        return None
    doubled = 2.0 * spin
    if spin < 0 or not math.isclose(doubled, round(doubled), abs_tol=1e-9):
        location = f" in row {row_number}" if row_number is not None else ""
        raise TSValueError(
            f"{field}{location} must be a non-negative integer or half-integer."
        )
    return spin


def _validate_redox_name(name):
    if not name or name in (".", "..") or Path(name).name != name:
        raise TSValueError(
            f"Invalid molecule name {name!r}; use a non-empty name without path separators."
        )
    return name


def _redox_job_from_structure(path, name, charge, spins):
    path = Path(path)
    if not path.is_file():
        raise TSValueError(f"Structure file not found: '{path}'.")
    if path.suffix.lower() not in _STRUCTURE_SUFFIXES:
        raise TSValueError(f"Redox structures must be .xyz or .gen files, got '{path}'.")
    charge = _integer_charge(charge, 0, "charge")
    spins = tuple(
        _valid_spin(spin, None, f"{state}_spin")
        for (state, _offset, _index), spin in zip(_REDOX_STATES, spins)
    )
    try:
        atoms = ase.io.read(str(path))
    except Exception as exc:  # pylint: disable=broad-except
        raise TSValueError(f"Could not read redox structure '{path}': {exc}") from exc
    nuclear_charge = int(sum(atoms.get_atomic_numbers()))
    for (state, charge_offset, spin_index) in _REDOX_STATES:
        spin = spins[spin_index]
        if spin is None:
            continue
        electrons = nuclear_charge - (charge + charge_offset)
        unpaired = int(round(2.0 * spin))
        if unpaired > electrons or (electrons - unpaired) % 2:
            raise TSValueError(
                f"{state}_spin={spin:g} is incompatible with {electrons} electrons "
                f"for '{name}'."
            )
    return _RedoxJob(_validate_redox_name(name), path.resolve(), charge, spins)


def _redox_job_from_smiles(
    smiles,
    name,
    charge,
    spins,
    generated_directory,
    max_conformers,
):
    name = _validate_redox_name(name)
    try:
        conformers = generate_conformers(smiles, max_conformers=max_conformers)
    except (ImportError, RuntimeError, ValueError) as exc:
        raise TSValueError(f"Could not generate {name!r} from SMILES: {exc}") from exc
    if not conformers:
        raise TSValueError(f"Could not generate {name!r} from SMILES: no conformers returned.")
    formal_charge = int(getattr(conformers[0], "info", {}).get("formal_charge", 0))
    charge = _integer_charge(charge, formal_charge, "charge")
    if formal_charge and charge != formal_charge:
        raise TSValueError(
            f"charge={charge} conflicts with formal charge {formal_charge} in "
            f"SMILES for {name!r}."
        )
    path = write_conformers(
        [conformers[0]], generated_directory, prefix=name
    )[0]
    return _redox_job_from_structure(path, name, charge, spins)


def _redox_jobs_from_manifest(
    manifest,
    charge,
    spin,
    generated_directory,
    max_conformers,
):
    jobs = []
    with open(manifest, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row_number, row in enumerate(reader, start=2):
            path_value = (row.get("path") or "").strip()
            smiles = (row.get("smiles") or "").strip()
            if bool(path_value) == bool(smiles):
                raise TSValueError(
                    f"Row {row_number} of '{manifest}' needs exactly one of "
                    "'path' or 'smiles'."
                )

            fallback_name = Path(path_value).stem if path_value else f"molecule-{row_number - 1}"
            name = _validate_redox_name((row.get("name") or fallback_name).strip())
            row_charge = _integer_charge(
                row.get("charge"), charge, "charge", row_number
            )
            oxidized_spin = _valid_spin(
                row.get("oxidized_spin") or row.get("spin"),
                spin,
                "oxidized_spin",
                row_number,
            )
            spins = (
                oxidized_spin,
                _valid_spin(
                    row.get("reduced_once_spin"),
                    None,
                    "reduced_once_spin",
                    row_number,
                ),
                _valid_spin(
                    row.get("reduced_twice_spin"),
                    None,
                    "reduced_twice_spin",
                    row_number,
                ),
            )

            if path_value:
                path = Path(path_value)
                if not path.is_absolute():
                    path = Path(manifest).parent / path
                job = _redox_job_from_structure(path, name, row_charge, spins)
            else:
                job = _redox_job_from_smiles(
                    smiles,
                    name,
                    row_charge,
                    spins,
                    generated_directory,
                    max_conformers,
                )
            jobs.append(job)
    return jobs


def _load_redox_jobs(
    source,
    charge,
    spin,
    generated_directory,
    max_conformers,
    default_name="molecule",
):
    source_path = Path(source)
    try:
        exists = source_path.exists()
    except OSError:
        exists = False

    if exists and source_path.is_dir():
        jobs = [
            _redox_job_from_structure(path, path.stem, charge, (spin, None, None))
            for path in sorted(source_path.iterdir())
            if path.suffix.lower() in _STRUCTURE_SUFFIXES
        ]
    elif exists and source_path.suffix.lower() == ".csv":
        jobs = _redox_jobs_from_manifest(
            source_path,
            charge,
            spin,
            generated_directory,
            max_conformers,
        )
    elif exists:
        jobs = [
            _redox_job_from_structure(
                source_path, source_path.stem, charge, (spin, None, None)
            )
        ]
    else:
        if isinstance(source, Path) or source_path.suffix.lower() in (
            ".csv",
            *_STRUCTURE_SUFFIXES,
        ):
            raise TSValueError(f"Redox input not found: '{source}'.")
        jobs = [
            _redox_job_from_smiles(
                str(source),
                default_name,
                charge,
                (spin, None, None),
                generated_directory,
                max_conformers,
            )
        ]

    if not jobs:
        raise TSValueError(f"No molecules found in redox input '{source}'.")
    names = [job.name for job in jobs]
    duplicates = sorted({name for name in names if names.count(name) > 1})
    if duplicates:
        raise TSValueError(
            f"Duplicate molecule name(s) in the redox screen: {', '.join(duplicates)}."
        )
    return jobs


def _resolve_redox_reference(
    reference,
    candidates,
    reference_charge,
    generated_directory,
    max_conformers,
):
    for candidate in candidates:
        if reference == candidate.name:
            if (
                reference_charge is not None
                and _integer_charge(reference_charge, None, "reference_charge")
                != candidate.charge
            ):
                raise TSValueError(
                    f"reference_charge does not match candidate {candidate.name!r}."
                )
            return candidate

    references = _load_redox_jobs(
        reference,
        reference_charge,
        None,
        generated_directory,
        max_conformers,
        default_name="reference",
    )
    if len(references) != 1:
        raise TSValueError("The redox reference must identify exactly one molecule.")
    selected = references[0]
    for candidate in candidates:
        if (
            candidate.path == selected.path
            and candidate.charge == selected.charge
            and candidate.spins == selected.spins
        ):
            return candidate
    return selected


def _write_redox_state_manifest(jobs, path):
    mapping = {}

    def write(handle):
        writer = csv.DictWriter(handle, fieldnames=("name", "path", "charge", "spin"))
        writer.writeheader()
        for job in jobs:
            for state, charge_offset, spin_index in _REDOX_STATES:
                state_name = f"{job.name}--{state}"
                mapping[(job.name, state)] = state_name
                writer.writerow(
                    {
                        "name": state_name,
                        "path": str(job.path),
                        "charge": job.charge + charge_offset,
                        "spin": "" if job.spins[spin_index] is None else job.spins[spin_index],
                    }
                )

    path.parent.mkdir(parents=True, exist_ok=True)
    _atomic_write(path, write)
    return mapping


def _redox_state_records(job, mapping, records):
    state_records = {}
    for state, _charge_offset, _spin_index in _REDOX_STATES:
        name = mapping[(job.name, state)]
        state_records[state] = records.get(
            name,
            {"status": "error", "error": "state calculation produced no result"},
        )
    return state_records


def _state_failure(states):
    failures = []
    for state, record in states.items():
        if record.get("status") != "ok":
            failures.append(f"{state}: {record.get('error') or 'calculation failed'}")
        elif record.get("G_total_hartree") is None:
            failures.append(f"{state}: total Gibbs free energy is missing")
        else:
            try:
                finite = math.isfinite(float(record["G_total_hartree"]))
            except (TypeError, ValueError):
                finite = False
            if not finite:
                failures.append(f"{state}: total Gibbs free energy is not finite")
    return "; ".join(failures)


def _absolute_stepwise_potentials(states):
    oxidized = float(states["oxidized"]["G_total_hartree"])
    reduced_once = float(states["reduced_once"]["G_total_hartree"])
    reduced_twice = float(states["reduced_twice"]["G_total_hartree"])
    return (
        -(reduced_once - oxidized) * HARTREE_TO_EV,
        -(reduced_twice - reduced_once) * HARTREE_TO_EV,
    )


def _write_redox_results(results, out):
    stem = Path(str(out)).with_suffix("")
    if stem.parent != Path(""):
        stem.parent.mkdir(parents=True, exist_ok=True)

    def write_csv(handle):
        writer = csv.DictWriter(handle, fieldnames=_REDOX_RESULT_FIELDS)
        writer.writeheader()
        for result in results:
            writer.writerow(
                {field: result.get(field, "") for field in _REDOX_RESULT_FIELDS}
            )

    csv_path = stem.with_suffix(".csv")
    json_path = stem.with_suffix(".json")
    _atomic_write(csv_path, write_csv)
    _atomic_write(json_path, lambda handle: handle.write(json.dumps(results, indent=2)))
    return csv_path, json_path


def _write_redox_run_metadata(
    *,
    out,
    source,
    candidates,
    reference_job,
    reference_e1,
    reference_e2,
    potential_scale,
    settings,
    max_conformers,
):
    def describe(job):
        if job is None:
            return None
        return {
            "name": job.name,
            "path": str(job.path),
            "structure_sha256": _file_sha256(job.path),
            "oxidized_charge": job.charge,
            "spins": {
                state: job.spins[index]
                for state, _offset, index in _REDOX_STATES
            },
        }

    payload = {
        "schema_version": 1,
        "workflow": "stepwise_reduction_screen",
        "source": str(source),
        "single_starting_geometry_approximation": True,
        "smiles_embedding": {
            "method": "ETKDGv3 followed by MMFF; lowest sampled conformer retained",
            "max_conformers": max_conformers,
        },
        "states": [
            {"name": state, "charge_offset": offset}
            for state, offset, _index in _REDOX_STATES
        ],
        "settings": settings,
        "candidates": [describe(job) for job in candidates],
        "reference": {
            "calculation": describe(reference_job),
            "experimental_E1_V": reference_e1,
            "experimental_E2_V": reference_e2,
            "potential_scale": potential_scale,
        },
    }
    run_fingerprint = _payload_fingerprint(payload)
    payload["run_fingerprint"] = run_fingerprint
    _write_json(_sidecar_path(out, "-run.json"), payload)
    return run_fingerprint


def redox_screen(
    source,
    out="redox-results",
    charge=None,
    temperature=298.15,
    pressure=101325,
    directory="redox-screening",
    parameters=None,
    spin=None,
    parameter_set=None,
    solvent=None,
    dispersion=None,
    quasi_rrho=False,
    engine="dftb+",
    method=None,
    resume=False,
    jobs=1,
    reference=None,
    reference_e1=None,
    reference_e2=None,
    reference_charge=None,
    potential_scale=None,
    max_conformers=20,
):
    """
    Run a two-step molecular reduction screen from one starting geometry.

    The workflow prepares one lowest-MMFF starting conformer from SMILES when
    needed, computes the oxidized, one-electron-reduced and two-electron-reduced
    states, and reports the first, second and overall two-electron reduction
    potentials. Each electronic state is optimized independently, but this is
    not a charge-state conformer search or ensemble-free-energy calculation. A
    reference molecule and measured E1/E2 may calibrate both steps.
    """

    calibration = (reference, reference_e1, reference_e2)
    if any(value is not None for value in calibration) and not all(
        value is not None for value in calibration
    ):
        raise TSValueError(
            "reference, reference_e1 and reference_e2 must be supplied together."
        )
    if potential_scale is not None and reference is None:
        raise TSValueError("potential_scale requires reference calibration.")
    if (
        isinstance(max_conformers, bool)
        or not isinstance(max_conformers, int)
        or max_conformers < 1
    ):
        raise TSValueError("max_conformers must be an integer >= 1.")
    if engine not in ("dftb+", "xtb", "xtb-cli"):
        raise TSValueError(
            f"Unknown engine {engine!r}; choose 'dftb+', 'xtb' or 'xtb-cli'."
        )
    if engine == "xtb" and solvent is not None:
        raise TSValueError("The tblite xTB engine does not support solvent.")
    if engine != "dftb+" and dispersion is not None:
        raise TSValueError("dispersion is supported only by the DFTB+ engine.")
    if engine != "dftb+" and parameters is not None:
        raise TSValueError("parameters are supported only by the DFTB+ engine.")
    if engine == "dftb+" and method is not None:
        raise TSValueError("method applies only to xTB engines.")
    if engine != "dftb+" and parameter_set is not None:
        raise TSValueError("parameter_set applies only to the DFTB+ engine.")

    resolved_parameter_set = parameter_set or "3ob"
    resolved_method = method or "GFN2-xTB"
    charge = _integer_charge(charge, None, "charge")
    spin = _valid_spin(spin, None, "spin")
    reference_charge = _integer_charge(
        reference_charge, None, "reference_charge"
    )
    reference_e1 = _optional_float(reference_e1, None, "reference_e1")
    reference_e2 = _optional_float(reference_e2, None, "reference_e2")
    root = Path(directory)
    generated_directory = root / "inputs"
    candidates = _load_redox_jobs(
        source,
        charge,
        spin,
        generated_directory,
        max_conformers,
    )

    reference_job = None
    calculation_jobs = list(candidates)
    if reference is not None:
        reference_job = _resolve_redox_reference(
            reference,
            candidates,
            reference_charge,
            generated_directory,
            max_conformers,
        )
        if reference_job not in calculation_jobs:
            occupied_names = {job.name for job in calculation_jobs}
            if reference_job.name in occupied_names:
                index = 1
                reference_name = "reference"
                while reference_name in occupied_names:
                    index += 1
                    reference_name = f"reference-{index}"
                reference_job = _RedoxJob(
                    reference_name,
                    reference_job.path,
                    reference_job.charge,
                    reference_job.spins,
                )
            calculation_jobs.append(reference_job)

    scale = (
        "SHE"
        if reference_job is None
        else potential_scale or f"calibrated:{reference_job.name}"
    )
    redox_settings = {
        "thermoscreening_version": __version__,
        "engine": engine,
        "temperature_K": temperature,
        "pressure_Pa": pressure,
        "method": resolved_method if engine in ("xtb", "xtb-cli") else None,
        "solvent": solvent,
        "dispersion": dispersion if engine == "dftb+" else None,
        "quasi_rrho": quasi_rrho,
        "parameter_set": resolved_parameter_set if engine == "dftb+" else None,
        "parameters": parameters if engine == "dftb+" else None,
    }
    run_fingerprint = _write_redox_run_metadata(
        out=out,
        source=source,
        candidates=candidates,
        reference_job=reference_job,
        reference_e1=reference_e1,
        reference_e2=reference_e2,
        potential_scale=scale,
        settings=redox_settings,
        max_conformers=max_conformers,
    )

    state_manifest = root / "states.csv"
    mapping = _write_redox_state_manifest(calculation_jobs, state_manifest)
    state_out = f"{Path(str(out)).with_suffix('')}-states"
    state_results = screen(
        state_manifest,
        out=state_out,
        temperature=temperature,
        pressure=pressure,
        directory=root / "calculations",
        parameters=parameters,
        parameter_set=resolved_parameter_set,
        solvent=solvent,
        dispersion=dispersion,
        quasi_rrho=quasi_rrho,
        engine=engine,
        method=resolved_method,
        resume=resume,
        jobs=jobs,
    )
    records = {record["name"]: record for record in state_results}

    reference_error = ""
    e1_reference = e2_reference = SHE_ABSOLUTE_POTENTIAL
    if reference_job is not None:
        reference_states = _redox_state_records(reference_job, mapping, records)
        reference_error = _state_failure(reference_states)
        if not reference_error:
            reference_abs_e1, reference_abs_e2 = _absolute_stepwise_potentials(
                reference_states
            )
            e1_reference = reference_abs_e1 - reference_e1
            e2_reference = reference_abs_e2 - reference_e2

    results = []
    for candidate in candidates:
        states = _redox_state_records(candidate, mapping, records)
        error = _state_failure(states)
        if reference_error:
            error = f"reference: {reference_error}" + (f"; {error}" if error else "")
        result = {
            "name": candidate.name,
            "path": str(candidate.path),
            "charge": candidate.charge,
            "status": "error" if error else "ok",
            "potential_scale": scale,
            "run_fingerprint": run_fingerprint,
            "error": error,
        }
        for state, field in (
            ("oxidized", "G_oxidized_hartree"),
            ("reduced_once", "G_reduced_once_hartree"),
            ("reduced_twice", "G_reduced_twice_hartree"),
        ):
            value = states[state].get("G_total_hartree")
            try:
                value = float(value)
            except (TypeError, ValueError):
                continue
            if math.isfinite(value):
                result[field] = value

        if not error:
            absolute_e1, absolute_e2 = _absolute_stepwise_potentials(states)
            e1 = absolute_e1 - e1_reference
            e2 = absolute_e2 - e2_reference
            result.update(
                {
                    "E1_V": e1,
                    "E2_V": e2,
                    "E2e_V": (e1 + e2) / 2.0,
                    "potential_gap_V": e1 - e2,
                    "potential_inversion": e2 >= e1,
                }
            )
        results.append(result)

    csv_path, json_path = _write_redox_results(results, out)
    logger.info(f"Wrote {csv_path}, {json_path} and the per-state results")
    return results
