import importlib

import ThermoScreening
import ThermoScreening.config as config


def test_log_file_env_works_without_logging_level(monkeypatch):
    monkeypatch.setenv("THERMOSCREENING_LOG_FILE", "on")
    monkeypatch.delenv("THERMOSCREENING_LOGGING_LEVEL", raising=False)
    monkeypatch.setattr(config, "use_log_file", False)
    monkeypatch.setattr(config, "log_file_name", None)

    importlib.reload(ThermoScreening)

    assert config.use_log_file is True
    assert config.log_file_name.startswith("ThermoScreening_")

    monkeypatch.delenv("THERMOSCREENING_LOG_FILE", raising=False)
    monkeypatch.setattr(config, "use_log_file", False)
    monkeypatch.setattr(config, "log_file_name", None)
    importlib.reload(ThermoScreening)
