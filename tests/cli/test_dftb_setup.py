import io
import tarfile

import pytest

from ThermoScreening.cli import dftb_setup
from ThermoScreening.cli.dftb_setup import (
    DEFAULT_PARAMETER_SET,
    REQUIRED_PARAMETER_FILE,
    Diagnostic,
    _solvent_stem,
    check_dftb_setup,
    default_parameter_dir,
    dftb_prefix_export,
    format_diagnostics,
    gbsa_param_path,
    install_gbsa_param,
    install_slakos,
    slako_url,
)


def _parameter_archive(tmp_path, files=None):
    files = files or {REQUIRED_PARAMETER_FILE: "parameter data"}
    source_dir = tmp_path / "source" / DEFAULT_PARAMETER_SET
    source_dir.mkdir(parents=True)

    for name, content in files.items():
        target = source_dir / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(content, encoding="utf-8")

    archive_path = tmp_path / f"{DEFAULT_PARAMETER_SET}.tar.xz"
    with tarfile.open(archive_path, "w:xz") as archive:
        archive.add(source_dir, arcname=DEFAULT_PARAMETER_SET)

    return archive_path


def test_install_slakos_extracts_archive(tmp_path):
    archive_path = _parameter_archive(tmp_path)

    installed_dir = install_slakos(
        install_root=tmp_path / "install",
        url=archive_path.as_uri(),
    )

    assert installed_dir == (tmp_path / "install" / DEFAULT_PARAMETER_SET).resolve()
    assert (installed_dir / REQUIRED_PARAMETER_FILE).read_text(
        encoding="utf-8"
    ) == "parameter data"


def test_slako_url_resolves_known_sets():
    assert "3ob" in slako_url("3ob")
    assert slako_url("3ob") == slako_url("3ob-3-1")
    assert "mio" in slako_url("mio")
    assert slako_url("mio") == slako_url("mio-1-1")


def test_slako_url_rejects_unknown_set():
    with pytest.raises(ValueError, match="Unknown parameter set"):
        slako_url("nope")


def test_default_parameter_dir_uses_canonical_name(tmp_path):
    assert default_parameter_dir(tmp_path, "mio") == tmp_path / "mio-1-1"
    assert default_parameter_dir(tmp_path, "3ob") == tmp_path / "3ob-3-1"


def test_install_slakos_installs_mio_set(tmp_path):
    source_dir = tmp_path / "source" / "mio-1-1"
    source_dir.mkdir(parents=True)
    (source_dir / REQUIRED_PARAMETER_FILE).write_text("mio data", encoding="utf-8")

    archive_path = tmp_path / "mio-1-1.tar.xz"
    with tarfile.open(archive_path, "w:xz") as archive:
        archive.add(source_dir, arcname="mio-1-1")

    installed_dir = install_slakos(
        install_root=tmp_path / "install",
        url=archive_path.as_uri(),
        parameter_set="mio",
    )

    assert installed_dir == (tmp_path / "install" / "mio-1-1").resolve()
    assert (installed_dir / REQUIRED_PARAMETER_FILE).read_text(
        encoding="utf-8"
    ) == "mio data"


def test_install_slakos_reuses_existing_directory(monkeypatch, tmp_path):
    marker_file = tmp_path / "install" / DEFAULT_PARAMETER_SET / REQUIRED_PARAMETER_FILE
    marker_file.parent.mkdir(parents=True)
    marker_file.write_text("existing", encoding="utf-8")

    def fail_download(*args, **kwargs):
        raise AssertionError("download should not run")

    monkeypatch.setattr(dftb_setup, "_download_file", fail_download)

    installed_dir = install_slakos(install_root=tmp_path / "install")

    assert installed_dir == marker_file.parent.resolve()
    assert marker_file.read_text(encoding="utf-8") == "existing"


def test_install_slakos_requires_marker_file(tmp_path):
    archive_path = _parameter_archive(tmp_path, files={"H-H.skf": "parameter data"})

    with pytest.raises(FileNotFoundError, match=REQUIRED_PARAMETER_FILE):
        install_slakos(install_root=tmp_path / "install", url=archive_path.as_uri())


def test_install_slakos_rejects_unsafe_archive_member(tmp_path):
    archive_path = tmp_path / f"{DEFAULT_PARAMETER_SET}.tar.xz"
    data = b"unsafe"

    with tarfile.open(archive_path, "w:xz") as archive:
        member = tarfile.TarInfo("../unsafe.txt")
        member.size = len(data)
        archive.addfile(member, io.BytesIO(data))

    with pytest.raises(ValueError, match="Unsafe archive member path"):
        install_slakos(install_root=tmp_path / "install", url=archive_path.as_uri())

    assert not (tmp_path / "unsafe.txt").exists()


def test_solvent_stem_resolves_names_and_aliases():
    assert _solvent_stem("water") == "h2o"
    assert _solvent_stem("Water") == "h2o"
    assert _solvent_stem("dichloromethane") == "ch2cl2"
    assert _solvent_stem("chloroform") == "chcl3"


def test_solvent_stem_rejects_unknown():
    with pytest.raises(ValueError, match="Unknown solvent"):
        _solvent_stem("unobtainium")


def test_gbsa_param_path_uses_solvent_stem(tmp_path):
    path = gbsa_param_path("water", install_root=tmp_path)
    assert path.name == "param_gbsa_h2o.txt"
    assert path.parent.name == "gfn2-0-1"


def test_install_gbsa_param_downloads_file(tmp_path):
    source = tmp_path / "param_gbsa_h2o.txt"
    source.write_text("solvent parameters", encoding="utf-8")

    installed = install_gbsa_param(
        "water", install_root=tmp_path / "install", url=source.as_uri()
    )

    assert installed == gbsa_param_path("water", tmp_path / "install").resolve()
    assert installed.read_text(encoding="utf-8") == "solvent parameters"


@pytest.mark.parametrize(
    "solvent,stem",
    [("water", "h2o"), ("dmso", "dmso"), ("acetonitrile", "acetonitrile"), ("thf", "thf")],
)
def test_gbsa_plumbing_for_multiple_solvents(monkeypatch, tmp_path, solvent, stem):
    assert _solvent_stem(solvent) == stem

    path = gbsa_param_path(solvent, install_root=tmp_path)
    assert path.name == f"param_gbsa_{stem}.txt"
    assert path.parent.name == "gfn2-0-1"

    # a mocked download lands the file at the resolved path
    def fake_download(url, destination):
        destination.write_text("params", encoding="utf-8")

    monkeypatch.setattr(dftb_setup, "_download_file", fake_download)
    installed = install_gbsa_param(solvent, install_root=tmp_path)
    assert installed == path.resolve()


def test_install_gbsa_param_requires_downloaded_file(monkeypatch, tmp_path):
    # a "download" that writes nothing must not silently succeed
    monkeypatch.setattr(dftb_setup, "_download_file", lambda url, destination: None)

    with pytest.raises(FileNotFoundError):
        install_gbsa_param("water", install_root=tmp_path / "install")


def test_install_gbsa_param_reuses_existing(monkeypatch, tmp_path):
    destination = gbsa_param_path("water", tmp_path / "install")
    destination.parent.mkdir(parents=True)
    destination.write_text("existing", encoding="utf-8")

    def fail_download(*args, **kwargs):
        raise AssertionError("download should not run")

    monkeypatch.setattr(dftb_setup, "_download_file", fail_download)

    installed = install_gbsa_param("water", install_root=tmp_path / "install")

    assert installed == destination.resolve()
    assert destination.read_text(encoding="utf-8") == "existing"


def test_dftb_prefix_export_adds_trailing_separator(tmp_path):
    expected = f'export DFTB_PREFIX="{tmp_path.resolve()}/"'

    assert dftb_prefix_export(tmp_path) == expected


def test_check_dftb_setup_reports_ready_environment(monkeypatch, tmp_path):
    marker_file = tmp_path / REQUIRED_PARAMETER_FILE
    marker_file.write_text("parameter data", encoding="utf-8")

    def fake_which(command):
        return f"/usr/bin/{command}"

    monkeypatch.setattr(dftb_setup.shutil, "which", fake_which)

    diagnostics = check_dftb_setup({"DFTB_PREFIX": str(tmp_path)})

    required = [(item.name, item.ok) for item in diagnostics if not item.optional]
    assert required == [
        ("dftb+", True),
        ("modes", True),
        ("DFTB_PREFIX", True),
        (REQUIRED_PARAMETER_FILE, True),
    ]


def test_check_dftb_setup_reports_missing_environment(monkeypatch):
    monkeypatch.setattr(dftb_setup.shutil, "which", lambda command: None)

    diagnostics = check_dftb_setup({})

    required = [(item.name, item.ok) for item in diagnostics if not item.optional]
    assert required == [
        ("dftb+", False),
        ("modes", False),
        ("DFTB_PREFIX", False),
        (REQUIRED_PARAMETER_FILE, False),
    ]


def test_check_dftb_setup_reports_xtb_toolchain_as_optional(monkeypatch):
    # xtb resolved via XTB_COMMAND; tblite importable
    monkeypatch.setattr(dftb_setup.shutil, "which", lambda command: "/bin/" + command)
    monkeypatch.setattr(
        dftb_setup.importlib.util, "find_spec",
        lambda name: object() if name == "tblite" else None,
    )

    optional = {
        item.name: item
        for item in check_dftb_setup({"XTB_COMMAND": "/opt/xtb"})
        if item.optional
    }
    assert set(optional) == {"xtb", "tblite"}
    assert optional["xtb"].ok is True
    assert optional["tblite"].ok is True

    # missing xtb toolchain -> reported, still optional
    monkeypatch.setattr(dftb_setup.shutil, "which", lambda command: None)
    monkeypatch.setattr(dftb_setup.importlib.util, "find_spec", lambda name: None)
    missing = {item.name: item for item in check_dftb_setup({}) if item.optional}
    assert missing["xtb"].ok is False and missing["xtb"].optional
    assert missing["tblite"].ok is False


def test_format_diagnostics_aligns_statuses():
    output = format_diagnostics(
        [
            Diagnostic("dftb+", True, "/usr/bin/dftb+"),
            Diagnostic("DFTB_PREFIX", False, "not set"),
            Diagnostic("xtb", False, "not found", optional=True),
        ]
    )

    assert "dftb+        found" in output
    assert "DFTB_PREFIX  missing" in output
    # a missing optional backend is marked so it doesn't read as a hard failure
    assert "xtb          missing  not found  (optional)" in output
