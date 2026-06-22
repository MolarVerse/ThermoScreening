import logging
import sys

import pytest

from ThermoScreening.utils.custom_logging import (
    CustomColorFormatter,
    CustomFormatter,
    CustomLogger,
    setup_logger,
)


def test_setup_logger_adds_stream_handler_once():
    logger = logging.getLogger("thermoscreening-test-setup")
    logger.handlers.clear()

    configured = setup_logger(logger)
    configured_again = setup_logger(logging.getLogger(logger.name))

    assert configured is configured_again
    assert len(configured.handlers) == 1
    assert configured.propagate is False
    assert isinstance(configured.handlers[0].formatter, CustomColorFormatter)

    configured.handlers.clear()


def test_custom_formatter_includes_exception_name():
    record = logging.LogRecord(
        name="thermoscreening-test",
        level=logging.ERROR,
        pathname=__file__,
        lineno=1,
        msg="  first line\nsecond line",
        args=(),
        exc_info=None,
    )
    record.custom_exception = ValueError

    formatted = CustomFormatter().format(record)

    assert "ERROR:" in formatted
    assert "thermoscreening-test - ValueError" in formatted
    assert "first line" in formatted
    assert "second line" in formatted


def test_custom_formatter_without_exception_name():
    record = logging.LogRecord(
        name="thermoscreening-test",
        level=logging.INFO,
        pathname=__file__,
        lineno=1,
        msg="  message",
        args=(),
        exc_info=None,
    )

    formatted = CustomFormatter().format(record)

    assert "INFO:" in formatted
    assert "thermoscreening-test" in formatted
    assert "message" in formatted


def test_color_formatter_uses_custom_format_for_unknown_level():
    record = logging.LogRecord(
        name="thermoscreening-test",
        level=35,
        pathname=__file__,
        lineno=1,
        msg="message",
        args=(),
        exc_info=None,
    )

    assert CustomColorFormatter().format_level(record).startswith("\x1b[36;1m")


def test_custom_logger_error_raises_in_debug_mode():
    logger = CustomLogger("thermoscreening-test-debug", level=logging.DEBUG)
    logger.addHandler(logging.NullHandler())

    with pytest.raises(ValueError, match="bad value"):
        logger.error("bad value", exception=ValueError)


def test_custom_logger_critical_raises_and_restores_hook():
    logger = CustomLogger("thermoscreening-test-critical", level=logging.INFO)
    logger.addHandler(logging.NullHandler())
    original_hook = sys.excepthook

    try:
        with pytest.raises(RuntimeError, match="critical value"):
            logger.critical("critical value", exception=RuntimeError)
    finally:
        sys.excepthook = original_hook


def test_custom_logger_original_methods_delegate(monkeypatch):
    logger = CustomLogger("thermoscreening-test-original", level=logging.DEBUG)
    calls = []

    def fake_original_log(level, msg, args, **kwargs):
        calls.append((level, msg, args, kwargs))

    monkeypatch.setattr(logger, "_original_log", fake_original_log)

    logger.original_error("error message", 1)
    logger.original_critical("critical message", 2)

    assert calls == [
        (logging.ERROR, "error message", (1,), {}),
        (logging.CRITICAL, "critical message", (2,), {}),
    ]
