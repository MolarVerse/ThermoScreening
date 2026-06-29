"""Package version helpers."""

from importlib.metadata import PackageNotFoundError, version


def get_version() -> str:
    try:
        from .__version__ import __version__
    except ModuleNotFoundError:
        try:
            return version("ThermoScreening")
        except PackageNotFoundError:
            return "0+unknown"
    return __version__


__version__ = get_version()
