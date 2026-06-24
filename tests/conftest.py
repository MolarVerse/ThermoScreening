import shutil
from pathlib import Path

import pytest


@pytest.fixture(scope="function")
def tmpdir(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    return str(tmp_path)


@pytest.fixture(scope="function")
def test_with_data_dir(example_dir, tmp_path, monkeypatch):
    source = Path(__file__).parent / "data" / example_dir
    workdir = tmp_path / example_dir

    shutil.copytree(source, workdir)
    monkeypatch.chdir(workdir)

    return str(workdir)
