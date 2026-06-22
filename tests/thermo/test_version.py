import builtins
import sys
import types

from importlib.metadata import PackageNotFoundError

import ThermoScreening.version as version_module


def test_get_version_uses_generated_module(monkeypatch):
    generated_version = types.ModuleType("ThermoScreening.__version__")
    generated_version.__version__ = "1.2.3"
    monkeypatch.setitem(sys.modules, "ThermoScreening.__version__", generated_version)

    assert version_module.get_version() == "1.2.3"


def test_get_version_falls_back_when_metadata_missing(monkeypatch):
    monkeypatch.delitem(sys.modules, "ThermoScreening.__version__", raising=False)
    real_import = builtins.__import__

    def import_without_generated_version(name, globals=None, locals=None, fromlist=(), level=0):
        if (level == 1 and name == "__version__") or name == "ThermoScreening.__version__":
            raise ModuleNotFoundError(name)
        return real_import(name, globals, locals, fromlist, level)

    def missing_version(package_name):
        raise PackageNotFoundError(package_name)

    monkeypatch.setattr(builtins, "__import__", import_without_generated_version)
    monkeypatch.setattr(version_module, "version", missing_version)

    assert version_module.get_version() == "0+unknown"
