"""Slurm job-array support for distributed thermochemistry screens."""

import os
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

from ThermoScreening.exceptions import TSValueError


_SLURM_TOKEN = re.compile(r"^[A-Za-z0-9_.:,/-]+$")


def _positive_integer(value, name):
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise TSValueError(f"{name} must be an integer >= 1.")
    return value


def _directive_value(value, name):
    if value is None:
        return None
    if not _SLURM_TOKEN.fullmatch(str(value)):
        raise TSValueError(f"{name} contains unsupported characters.")
    return str(value)


def write_slurm_array_script(
    command_args,
    *,
    tasks,
    script,
    local_jobs=1,
    cpus_per_task=1,
    job_name="thermoscreening",
    walltime=None,
    memory=None,
    partition=None,
    account=None,
    preamble=None,
    working_directory=None,
    python_executable=None,
):
    """Write a Slurm array that runs one deterministic workflow shard per task."""
    tasks = _positive_integer(tasks, "tasks")
    local_jobs = _positive_integer(local_jobs, "local_jobs")
    cpus_per_task = _positive_integer(cpus_per_task, "cpus_per_task")
    if local_jobs > cpus_per_task:
        raise TSValueError("local --jobs cannot exceed cpus_per_task.")
    if not command_args or command_args[0] not in {"screen", "redox"}:
        raise TSValueError("Slurm arrays support the screen and redox commands.")
    if "--shard-index" in command_args or "--shard-count" in command_args:
        raise TSValueError("Do not pass shard options to the Slurm generator.")

    job_name = _directive_value(job_name, "job_name")
    optional_directives = (
        ("time", _directive_value(walltime, "time")),
        ("mem", _directive_value(memory, "memory")),
        ("partition", _directive_value(partition, "partition")),
        ("account", _directive_value(account, "account")),
    )
    script = Path(script).resolve()
    script.parent.mkdir(parents=True, exist_ok=True)
    working_directory = Path(working_directory or os.getcwd()).resolve()
    executable = str(python_executable or sys.executable)
    threads_per_job = max(1, cpus_per_task // local_jobs)

    lines = [
        "#!/usr/bin/env bash",
        f"#SBATCH --job-name={job_name}",
        f"#SBATCH --array=0-{tasks - 1}",
        f"#SBATCH --cpus-per-task={cpus_per_task}",
        f"#SBATCH --output={job_name}-%A_%a.out",
    ]
    lines.extend(
        f"#SBATCH --{option}={value}"
        for option, value in optional_directives
        if value is not None
    )
    lines.extend(
        [
            "",
            "set -euo pipefail",
            f"cd -- {shlex.quote(str(working_directory))}",
            f"export OMP_NUM_THREADS={threads_per_job}",
            f"export OPENBLAS_NUM_THREADS={threads_per_job}",
            f"export MKL_NUM_THREADS={threads_per_job}",
            f"export NUMEXPR_NUM_THREADS={threads_per_job}",
        ]
    )
    if preamble is not None:
        preamble_text = Path(preamble).read_text(encoding="utf-8").rstrip()
        if preamble_text:
            lines.extend(["", preamble_text])

    command = shlex.join([executable, "-m", "ThermoScreening", *command_args])
    lines.extend(
        [
            "",
            command
            + ' --shard-index "$SLURM_ARRAY_TASK_ID"'
            + ' --shard-count "$SLURM_ARRAY_TASK_COUNT"',
            "",
        ]
    )
    script.write_text("\n".join(lines), encoding="utf-8")
    script.chmod(script.stat().st_mode | 0o111)
    return script


def submit_slurm_array(
    script,
    *,
    shard_directory,
    out,
    job_name="thermoscreening",
    partition=None,
    account=None,
    working_directory=None,
    python_executable=None,
):
    """Submit a Slurm array and a dependent result-collection job."""
    sbatch = shutil.which("sbatch")
    if sbatch is None:
        raise TSValueError("sbatch was not found on PATH.")
    job_name = _directive_value(job_name, "job_name")
    partition = _directive_value(partition, "partition")
    account = _directive_value(account, "account")
    working_directory = Path(working_directory or os.getcwd()).resolve()
    executable = str(python_executable or sys.executable)

    try:
        array = subprocess.run(
            [sbatch, "--parsable", str(Path(script).resolve())],
            check=True,
            capture_output=True,
            text=True,
        )
        array_job_id = array.stdout.strip().split(";", maxsplit=1)[0]
        if not array_job_id:
            raise TSValueError("sbatch returned no array job ID.")
        collect_command = shlex.join(
            [
                executable,
                "-m",
                "ThermoScreening",
                "collect",
                str(shard_directory),
                "-o",
                str(out),
            ]
        )
        collector_args = [
            sbatch,
            "--parsable",
            f"--dependency=afterany:{array_job_id}",
            f"--job-name={job_name}-collect",
            f"--chdir={working_directory}",
            f"--output={job_name}-collect-%j.out",
        ]
        if partition is not None:
            collector_args.append(f"--partition={partition}")
        if account is not None:
            collector_args.append(f"--account={account}")
        collector_args.append(f"--wrap={collect_command}")
        collector = subprocess.run(
            collector_args,
            check=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError as exc:
        detail = (exc.stderr or exc.stdout or str(exc)).strip()
        raise TSValueError(f"Slurm submission failed: {detail}") from exc

    collector_job_id = collector.stdout.strip().split(";", maxsplit=1)[0]
    if not collector_job_id:
        raise TSValueError("sbatch returned no collector job ID.")
    return array_job_id, collector_job_id
