"""Shared locations and names for external DFTB+ parameter data."""

from pathlib import Path


DEFAULT_PARAMETER_SET = "3ob-3-1"
REQUIRED_PARAMETER_FILE = "C-C.skf"
PARAMETER_SET_ALIASES = {"3ob": "3ob-3-1", "mio": "mio-1-1"}


def canonical_parameter_set(parameter_set):
    """Return the canonical directory name for a parameter-set alias."""
    return PARAMETER_SET_ALIASES.get(parameter_set, parameter_set)


def default_parameter_root():
    """Return the user-local directory containing parameter sets."""
    return Path.home() / ".local" / "share" / "thermoscreening" / "slakos"


def default_parameter_dir(parameter_set=DEFAULT_PARAMETER_SET, install_root=None):
    """Return the user-local directory for a parameter set."""
    root = Path(install_root).expanduser() if install_root else default_parameter_root()
    return root / canonical_parameter_set(parameter_set)
