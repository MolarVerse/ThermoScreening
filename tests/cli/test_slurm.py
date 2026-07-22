from argparse import Namespace
from types import SimpleNamespace

import pytest

from ThermoScreening.cli import slurm
from ThermoScreening.exceptions import TSValueError


def test_write_slurm_array_script_sets_resources_and_shards(tmp_path):
    preamble = tmp_path / "preamble.sh"
    preamble.write_text("module load dftbplus\nsource env/bin/activate\n", encoding="utf-8")

    script = slurm.write_slurm_array_script(
        ["screen", "molecules.csv", "-o", "results", "--jobs", "2"],
        tasks=8,
        script=tmp_path / "screen.slurm",
        local_jobs=2,
        cpus_per_task=8,
        job_name="aq-screen",
        walltime="04:00:00",
        memory="16G",
        partition="compute",
        account="chemistry",
        preamble=preamble,
        working_directory=tmp_path,
        python_executable="/shared/python",
    )

    text = script.read_text(encoding="utf-8")
    assert "#SBATCH --array=0-7" in text
    assert "#SBATCH --cpus-per-task=8" in text
    assert "#SBATCH --time=04:00:00" in text
    assert "#SBATCH --mem=16G" in text
    assert "export OMP_NUM_THREADS=4" in text
    assert "module load dftbplus" in text
    assert "/shared/python -m ThermoScreening screen molecules.csv" in text
    assert '--shard-index "$SLURM_ARRAY_TASK_ID"' in text
    assert '--shard-count "$SLURM_ARRAY_TASK_COUNT"' in text
    assert script.stat().st_mode & 0o111


def test_write_slurm_array_script_supports_redox(tmp_path):
    script = slurm.write_slurm_array_script(
        ["redox", "molecules.csv", "-o", "potentials"],
        tasks=4,
        script=tmp_path / "redox.slurm",
        working_directory=tmp_path,
        python_executable="/shared/python",
    )

    text = script.read_text(encoding="utf-8")
    assert "/shared/python -m ThermoScreening redox molecules.csv" in text
    assert '--shard-index "$SLURM_ARRAY_TASK_ID"' in text


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"tasks": 0}, "tasks must be"),
        ({"tasks": 2, "local_jobs": 3, "cpus_per_task": 2}, "cannot exceed"),
        ({"tasks": 2, "job_name": "bad\nname"}, "unsupported characters"),
    ],
)
def test_write_slurm_array_script_rejects_invalid_resources(tmp_path, kwargs, message):
    options = {"tasks": 2, "script": tmp_path / "screen.slurm", **kwargs}
    with pytest.raises(TSValueError, match=message):
        slurm.write_slurm_array_script(["screen", "molecules.csv"], **options)


def test_write_slurm_array_script_rejects_manual_shard_options(tmp_path):
    with pytest.raises(TSValueError, match="Do not pass shard options"):
        slurm.write_slurm_array_script(
            ["screen", "molecules.csv", "--shard-index", "0"],
            tasks=2,
            script=tmp_path / "screen.slurm",
        )


def test_submit_slurm_array_schedules_dependent_collector(monkeypatch, tmp_path):
    calls = []
    responses = iter(
        [SimpleNamespace(stdout="12345;cluster\n"), SimpleNamespace(stdout="12346\n")]
    )

    monkeypatch.setattr(slurm.shutil, "which", lambda command: "/usr/bin/sbatch")

    def fake_run(args, **kwargs):
        calls.append((args, kwargs))
        return next(responses)

    monkeypatch.setattr(slurm.subprocess, "run", fake_run)
    script = tmp_path / "screen.slurm"
    script.write_text("#!/bin/bash\n", encoding="utf-8")

    result = slurm.submit_slurm_array(
        script,
        shard_directory=tmp_path / "results-shards",
        out=tmp_path / "results",
        job_name="screen",
        partition="compute",
        account="chemistry",
        working_directory=tmp_path,
        python_executable="/shared/python",
    )

    assert result == ("12345", "12346")
    assert calls[0][0] == ["/usr/bin/sbatch", "--parsable", str(script.resolve())]
    collector = calls[1][0]
    assert "--dependency=afterany:12345" in collector
    assert "--partition=compute" in collector
    assert "--account=chemistry" in collector
    assert any("ThermoScreening collect" in argument for argument in collector)


def test_submit_slurm_array_requires_sbatch(monkeypatch, tmp_path):
    monkeypatch.setattr(slurm.shutil, "which", lambda command: None)
    with pytest.raises(TSValueError, match="sbatch was not found"):
        slurm.submit_slurm_array(
            tmp_path / "screen.slurm",
            shard_directory=tmp_path / "shards",
            out=tmp_path / "results",
        )


def test_cli_generates_slurm_script(tmp_path, capsys):
    import ThermoScreening.cli.thermo as cli

    args = Namespace(
        tasks=4,
        script=str(tmp_path / "screen.slurm"),
        job_name="screen",
        cpus_per_task=2,
        walltime=None,
        memory=None,
        partition=None,
        account=None,
        preamble=None,
        submit=False,
        command_args=["--", "screen", "molecules.csv", "-o", "results"],
    )

    assert cli.run_slurm(args) == 0
    output = capsys.readouterr().out
    assert "Slurm script:" in output
    assert "thermo collect results-shards -o results" in output


def test_cli_collect_reports_failures(monkeypatch, capsys):
    import ThermoScreening.cli.thermo as cli

    monkeypatch.setattr(
        cli,
        "collect_shards",
        lambda *args, **kwargs: [
            {"name": "ok", "status": "ok"},
            {"name": "bad", "status": "error", "error": "failed"},
        ],
    )
    args = Namespace(shard_directory="results-shards", out="results")

    assert cli.run_collect(args) == 1
    assert "Collected 2 molecules (1 failed)" in capsys.readouterr().out


def test_slurm_parser_accepts_nested_screen_command():
    import ThermoScreening.cli.thermo as cli

    args = cli.parse_args(
        [
            "slurm",
            "--tasks",
            "8",
            "--cpus-per-task",
            "4",
            "--",
            "screen",
            "molecules.csv",
            "--jobs",
            "2",
        ]
    )

    assert args.command == "slurm"
    assert args.tasks == 8
    assert args.command_args[-3:] == ["molecules.csv", "--jobs", "2"]


def test_slurm_parser_accepts_nested_redox_command():
    import ThermoScreening.cli.thermo as cli

    args = cli.parse_args(
        ["slurm", "--tasks", "8", "--", "redox", "molecules.csv", "--jobs", "2"]
    )

    assert args.command == "slurm"
    assert args.command_args[-3:] == ["molecules.csv", "--jobs", "2"]
