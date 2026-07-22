"""
Helpers for installing and validating external DFTB+ dependencies.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
import importlib.util
import os
import shutil
import subprocess
import tarfile
import tempfile
from urllib.request import urlopen

from ThermoScreening.dftb_data import (
    DEFAULT_PARAMETER_SET,
    REQUIRED_PARAMETER_FILE,
    canonical_parameter_set,
    default_parameter_dir as calculator_parameter_dir,
    default_parameter_root,
)


# Download URLs for the supported Slater-Koster sets (dftbparams GitHub releases).
PARAMETER_SET_URLS = {
    "3ob-3-1": (
        "https://github.com/dftbparams/3ob/releases/latest/download/3ob-3-1.tar.xz"
    ),
    "mio-1-1": (
        "https://github.com/dftbparams/mio/releases/latest/download/mio-1-1.tar.xz"
    ),
}

def _canonical_set_name(parameter_set: str) -> str:
    """
    Map a short set name (e.g. ``"mio"``) to its canonical name (``"mio-1-1"``).
    """

    return canonical_parameter_set(parameter_set)


def slako_url(parameter_set: str = DEFAULT_PARAMETER_SET) -> str:
    """
    Return the download URL for a Slater-Koster parameter set.

    Raises
    ------
    ValueError
        If ``parameter_set`` is not a known set.
    """

    name = _canonical_set_name(parameter_set)
    try:
        return PARAMETER_SET_URLS[name]
    except KeyError:
        known = ", ".join(sorted(PARAMETER_SET_URLS))
        raise ValueError(
            f"Unknown parameter set {parameter_set!r}; choose one of: {known}."
        )


DEFAULT_SLAKO_URL = PARAMETER_SET_URLS[DEFAULT_PARAMETER_SET]


@dataclass(frozen=True)
class Diagnostic:
    """
    Result for one DFTB+ setup check.
    """

    name: str
    ok: bool
    detail: str
    optional: bool = False


def default_install_root() -> Path:
    """
    Return the default user-local Slater-Koster install directory.
    """

    return default_parameter_root()


def default_parameter_dir(
    install_root: str | Path | None = None,
    parameter_set: str = DEFAULT_PARAMETER_SET,
) -> Path:
    """
    Return the parameter directory for ``parameter_set`` under an install root.
    """

    return calculator_parameter_dir(parameter_set, install_root)


def _download_file(url: str, destination: Path, timeout: int = 60) -> None:
    """
    Download a URL to a local path.
    """

    with urlopen(url, timeout=timeout) as response:  # nosec B310 - user-visible setup helper
        with destination.open("wb") as output_file:
            shutil.copyfileobj(response, output_file)


def _safe_extract_tar(archive_path: Path, destination: Path) -> None:
    """
    Extract a tar archive while rejecting path traversal entries.
    """

    destination_resolved = destination.resolve()

    with tarfile.open(archive_path, "r:xz") as archive:
        for member in archive.getmembers():
            member_path = destination_resolved / member.name
            if not member_path.resolve().is_relative_to(destination_resolved):
                raise ValueError(f"Unsafe archive member path: {member.name}")

        archive.extractall(destination_resolved, filter="data")


def install_slakos(
    install_root: str | Path | None = None,
    url: str | None = None,
    force: bool = False,
    parameter_set: str = DEFAULT_PARAMETER_SET,
) -> Path:
    """
    Download and extract a Slater-Koster parameter set.

    Parameters
    ----------
    install_root : str or Path, optional
        Directory where parameter sets are installed. Defaults to the user-local
        share directory.
    url : str, optional
        Archive URL. Defaults to the release URL for ``parameter_set``.
    force : bool
        Re-download even when the set already exists.
    parameter_set : str
        Set to install (``"3ob"``/``"3ob-3-1"`` or ``"mio"``/``"mio-1-1"``).
    """

    name = _canonical_set_name(parameter_set)
    if url is None:
        url = slako_url(name)

    root = Path(install_root).expanduser() if install_root is not None else default_install_root()
    parameter_dir = default_parameter_dir(root, name)
    marker_file = parameter_dir / REQUIRED_PARAMETER_FILE

    if marker_file.exists() and not force:
        return parameter_dir.resolve()

    root.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="thermoscreening-dftb-") as tmp_dir:
        archive_path = Path(tmp_dir) / f"{name}.tar.xz"
        _download_file(url, archive_path)
        _safe_extract_tar(archive_path, root)

    if not marker_file.exists():
        raise FileNotFoundError(
            f"Downloaded parameter set is missing {REQUIRED_PARAMETER_FILE}: "
            f"{parameter_dir}"
        )

    return parameter_dir.resolve()


# GBSA/ALPB implicit-solvation parameter files (grimme-lab/gbsa-parameters). The
# published sets are fit for GFN-xTB; used with DFTB (3ob/mio) they are an
# approximation (see the note in ``calculator.dftbplus._solvation_kwargs``).
GBSA_PARAM_METHOD = "gfn2-0-1"
GBSA_BASE_URL = "https://raw.githubusercontent.com/grimme-lab/gbsa-parameters/main"

# User-facing solvent name -> parameter-file stem in the repository.
GBSA_SOLVENTS = {
    "acetone": "acetone",
    "acetonitrile": "acetonitrile",
    "benzene": "benzene",
    "ch2cl2": "ch2cl2",
    "dichloromethane": "ch2cl2",
    "chcl3": "chcl3",
    "chloroform": "chcl3",
    "cs2": "cs2",
    "dmf": "dmf",
    "dmso": "dmso",
    "ether": "ether",
    "diethylether": "ether",
    "water": "h2o",
    "h2o": "h2o",
    "methanol": "methanol",
    "hexane": "nhexane",
    "nhexane": "nhexane",
    "thf": "thf",
    "toluene": "toluene",
}


def _solvent_stem(solvent: str) -> str:
    """
    Map a user-facing solvent name to its parameter-file stem.

    Raises
    ------
    ValueError
        If ``solvent`` is not a known solvent.
    """

    try:
        return GBSA_SOLVENTS[solvent.lower()]
    except KeyError:
        known = ", ".join(sorted(set(GBSA_SOLVENTS)))
        raise ValueError(
            f"Unknown solvent {solvent!r}; choose one of: {known}."
        )


def default_gbsa_dir(install_root: str | Path | None = None) -> Path:
    """
    Return the directory holding downloaded GBSA solvation parameter files.
    """

    root = Path(install_root).expanduser() if install_root is not None else default_install_root()
    return root / "gbsa" / GBSA_PARAM_METHOD


def gbsa_param_path(
    solvent: str,
    install_root: str | Path | None = None,
) -> Path:
    """
    Return the local path where the GBSA parameter file for ``solvent`` lives
    (whether or not it has been downloaded yet).
    """

    return default_gbsa_dir(install_root) / f"param_gbsa_{_solvent_stem(solvent)}.txt"


def install_gbsa_param(
    solvent: str,
    install_root: str | Path | None = None,
    url: str | None = None,
    force: bool = False,
) -> Path:
    """
    Download the GBSA implicit-solvation parameter file for ``solvent``.

    Parameters
    ----------
    solvent : str
        Solvent name (e.g. ``"water"``); see :data:`GBSA_SOLVENTS`.
    install_root : str or Path, optional
        Directory root for downloaded parameters. Defaults to the user-local
        share directory.
    url : str, optional
        Explicit file URL override. Defaults to the grimme-lab release file.
    force : bool
        Re-download even when the file already exists.
    """

    stem = _solvent_stem(solvent)
    destination = gbsa_param_path(solvent, install_root)

    if destination.exists() and not force:
        return destination.resolve()

    if url is None:
        url = f"{GBSA_BASE_URL}/{GBSA_PARAM_METHOD}/param_gbsa_{stem}.txt"

    destination.parent.mkdir(parents=True, exist_ok=True)
    _download_file(url, destination)

    if not destination.exists():
        raise FileNotFoundError(
            f"Downloaded GBSA parameter file is missing: {destination}"
        )

    return destination.resolve()


def dftb_prefix_export(parameter_dir: str | Path) -> str:
    """
    Return the shell export line for a Slater-Koster parameter directory.
    """

    return f'export DFTB_PREFIX="{Path(parameter_dir).expanduser().resolve()}{os.sep}"'


def _xtb_diagnostics(env: Mapping[str, str]) -> list[Diagnostic]:
    """
    Optional checks for the xTB engines (``--engine xtb`` and ``xtb-cli``).
    """

    # native xtb binary (for --engine xtb-cli): honour XTB_COMMAND, then PATH
    xtb_command = env.get("XTB_COMMAND") or "xtb"
    xtb = _executable_diagnostic(
        "xtb",
        xtb_command,
        expected="xtb",
        optional=True,
        missing_detail=(
            "not found (set XTB_COMMAND or `conda install -c conda-forge xtb`; "
            "needed for --engine xtb-cli)"
        ),
    )

    # tblite python package (for the in-process --engine xtb)
    tblite_ok = importlib.util.find_spec("tblite") is not None
    tblite = Diagnostic(
        "tblite",
        tblite_ok,
        "importable" if tblite_ok
        else "not importable (`conda install -c conda-forge tblite-python`; "
        "needed for --engine xtb)",
        optional=True,
    )

    return [xtb, tblite]


def _executable_diagnostic(
    name,
    command=None,
    *,
    expected=None,
    optional=False,
    missing_detail="not found on PATH",
):
    """Check that an executable is present and can start."""
    executable = shutil.which(command or name)
    if executable is None:
        return Diagnostic(name, False, missing_detail, optional=optional)
    try:
        completed = subprocess.run(
            [executable, "--help"],
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return Diagnostic(name, False, f"could not start: {exc}", optional=optional)

    output = f"{completed.stdout}\n{completed.stderr}"
    loader_errors = ("Library not loaded", "error while loading shared libraries")
    starts = not any(marker in output for marker in loader_errors)
    if expected is not None:
        starts = starts and expected.lower() in output.lower()
    else:
        starts = starts and completed.returncode == 0
    detail = executable if starts else "found but could not start correctly"
    return Diagnostic(name, starts, detail, optional=optional)


def check_dftb_setup(
    env: dict[str, str] | None = None,
    parameter_set: str = "3ob",
) -> list[Diagnostic]:
    """
    Check whether the calculation backends are available.

    Reports the DFTB+ toolchain (dftb+, modes, DFTB_PREFIX + a Slater-Koster
    file) as required, and the xTB toolchain (xtb binary, tblite) as optional.
    """

    current_env = os.environ if env is None else env
    prefix = current_env.get("DFTB_PREFIX")
    diagnostics = [
        _executable_diagnostic(
            "dftb+",
            expected="DFTB+",
        ),
        _executable_diagnostic(
            "modes",
            expected="DFTB+",
        ),
    ]

    if not prefix:
        parameter_dir = default_parameter_dir(parameter_set=parameter_set)
        source = f"downloaded {canonical_parameter_set(parameter_set)}"
    else:
        parameter_dir = Path(prefix).expanduser()
        source = "DFTB_PREFIX"
    parameter_file = parameter_dir / REQUIRED_PARAMETER_FILE
    diagnostics.extend(
        [
            Diagnostic(
                "parameters",
                parameter_dir.is_dir(),
                f"{source}: {parameter_dir.resolve()}"
                if parameter_dir.is_dir()
                else f"{source} directory not found",
            ),
            Diagnostic(
                REQUIRED_PARAMETER_FILE,
                parameter_file.is_file(),
                str(parameter_file.resolve()) if parameter_file.is_file() else "not found",
            ),
        ]
    )

    diagnostics.extend(_xtb_diagnostics(current_env))
    return diagnostics


def format_diagnostics(diagnostics: list[Diagnostic]) -> str:
    """
    Format setup diagnostics for terminal output.
    """

    lines = []
    width = max(len(item.name) for item in diagnostics)

    for item in diagnostics:
        status = "found" if item.ok else "missing"
        suffix = "  (optional)" if item.optional and not item.ok else ""
        lines.append(f"{item.name:<{width}}  {status:<7}  {item.detail}{suffix}")

    return "\n".join(lines)
