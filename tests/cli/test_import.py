"""Command-line package import behavior."""

import subprocess
import sys


def test_importing_cli_is_silent():
    completed = subprocess.run(
        [sys.executable, "-c", "import ThermoScreening.cli"],
        check=True,
        capture_output=True,
        text=True,
    )

    assert completed.stdout == ""
    assert completed.stderr == ""
