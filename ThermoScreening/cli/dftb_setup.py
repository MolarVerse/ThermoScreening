"""
Helpers for installing and validating external DFTB+ dependencies.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import os
import shutil
import tarfile
import tempfile
from urllib.request import urlopen


DEFAULT_PARAMETER_SET = "3ob-3-1"
DEFAULT_SLAKO_URL = (
    "https://github.com/dftbparams/3ob/releases/latest/download/"
    f"{DEFAULT_PARAMETER_SET}.tar.xz"
)
REQUIRED_PARAMETER_FILE = "C-C.skf"


@dataclass(frozen=True)
class Diagnostic:
    """
    Result for one DFTB+ setup check.
    """

    name: str
    ok: bool
    detail: str


def default_install_root() -> Path:
    """
    Return the default user-local Slater-Koster install directory.
    """

    return Path.home() / ".local" / "share" / "thermoscreening" / "slakos"


def default_parameter_dir(install_root: str | Path | None = None) -> Path:
    """
    Return the default 3ob parameter directory under an install root.
    """

    root = Path(install_root).expanduser() if install_root is not None else default_install_root()
    return root / DEFAULT_PARAMETER_SET


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
    url: str = DEFAULT_SLAKO_URL,
    force: bool = False,
) -> Path:
    """
    Download and extract the default Slater-Koster parameter set.
    """

    root = Path(install_root).expanduser() if install_root is not None else default_install_root()
    parameter_dir = default_parameter_dir(root)
    marker_file = parameter_dir / REQUIRED_PARAMETER_FILE

    if marker_file.exists() and not force:
        return parameter_dir.resolve()

    root.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="thermoscreening-dftb-") as tmp_dir:
        archive_path = Path(tmp_dir) / f"{DEFAULT_PARAMETER_SET}.tar.xz"
        _download_file(url, archive_path)
        _safe_extract_tar(archive_path, root)

    if not marker_file.exists():
        raise FileNotFoundError(
            f"Downloaded parameter set is missing {REQUIRED_PARAMETER_FILE}: "
            f"{parameter_dir}"
        )

    return parameter_dir.resolve()


def dftb_prefix_export(parameter_dir: str | Path) -> str:
    """
    Return the shell export line for a Slater-Koster parameter directory.
    """

    return f'export DFTB_PREFIX="{Path(parameter_dir).expanduser().resolve()}{os.sep}"'


def check_dftb_setup(env: dict[str, str] | None = None) -> list[Diagnostic]:
    """
    Check whether DFTB+ executables and parameters are available.
    """

    current_env = os.environ if env is None else env
    prefix = current_env.get("DFTB_PREFIX")
    diagnostics = [
        Diagnostic(
            "dftb+",
            shutil.which("dftb+") is not None,
            shutil.which("dftb+") or "not found on PATH",
        ),
        Diagnostic(
            "modes",
            shutil.which("modes") is not None,
            shutil.which("modes") or "not found on PATH",
        ),
    ]

    if not prefix:
        diagnostics.extend(
            [
                Diagnostic("DFTB_PREFIX", False, "not set"),
                Diagnostic(REQUIRED_PARAMETER_FILE, False, "DFTB_PREFIX is not set"),
            ]
        )
        return diagnostics

    parameter_dir = Path(prefix).expanduser()
    parameter_file = parameter_dir / REQUIRED_PARAMETER_FILE
    diagnostics.extend(
        [
            Diagnostic(
                "DFTB_PREFIX",
                parameter_dir.is_dir(),
                str(parameter_dir.resolve()) if parameter_dir.is_dir() else "directory not found",
            ),
            Diagnostic(
                REQUIRED_PARAMETER_FILE,
                parameter_file.is_file(),
                str(parameter_file.resolve()) if parameter_file.is_file() else "not found",
            ),
        ]
    )

    return diagnostics


def format_diagnostics(diagnostics: list[Diagnostic]) -> str:
    """
    Format setup diagnostics for terminal output.
    """

    lines = []
    width = max(len(item.name) for item in diagnostics)

    for item in diagnostics:
        status = "found" if item.ok else "missing"
        lines.append(f"{item.name:<{width}}  {status:<7}  {item.detail}")

    return "\n".join(lines)
