import pytest

import os
import subprocess
from pathlib import Path
import numpy as np
import ase.io as ase_io
from ase import Atoms
import ThermoScreening.calculator.dftbplus as dftbplus_module
from ThermoScreening.calculator import Geoopt, Hessian, Modes
from ThermoScreening.calculator.dftbplus import _slako_dir, dftb_3ob_parameters
from ThermoScreening.thermo.api import dftbplus_thermo

# --------------------------------------------------------------------------- #


def executable_starts(command):
    try:
        result = subprocess.run(
            [command, "--help"],
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False
    output = result.stdout + result.stderr
    return "DFTB+" in output and "Library not loaded" not in output


def slako_dir_ready():
    try:
        _slako_dir()
    except FileNotFoundError:
        return False
    return True


dftbplus_ready = (
    executable_starts("dftb+")
    and executable_starts("modes")
    and slako_dir_ready()
)


def test_slako_dir_uses_explicit_directory(tmp_path):
    assert _slako_dir(str(tmp_path)) == str(tmp_path.resolve()) + os.sep


def test_slako_dir_uses_environment(monkeypatch, tmp_path):
    monkeypatch.setenv("DFTB_PREFIX", str(tmp_path))

    assert _slako_dir() == str(tmp_path.resolve()) + os.sep


def test_slako_dir_requires_external_parameters(monkeypatch):
    monkeypatch.delenv("DFTB_PREFIX", raising=False)

    with pytest.raises(FileNotFoundError, match="Slater-Koster files are not bundled"):
        _slako_dir()


def test_slako_dir_rejects_missing_directory(monkeypatch, tmp_path):
    missing_dir = tmp_path / "missing"
    monkeypatch.setenv("DFTB_PREFIX", str(missing_dir))

    with pytest.raises(FileNotFoundError, match="Slater-Koster directory does not exist"):
        _slako_dir()


def test_modes_missing_executable(monkeypatch):
    monkeypatch.setattr("ThermoScreening.calculator.dftbplus.shutil.which", lambda command: None)
    modes = Modes.__new__(Modes)

    with pytest.raises(FileNotFoundError, match="modes executable"):
        modes.calculate()


def _stub_dftb_base(monkeypatch):
    def fake_init(self, atoms, label, slako_dir, **kwargs):
        self.atoms = atoms
        self.label = label
        self.slako_dir = slako_dir
        self.parameters = kwargs

    def fake_calculate(self, atoms):
        self.calculated_atoms = atoms

    monkeypatch.setattr(dftbplus_module.Dftb, "__init__", fake_init)
    monkeypatch.setattr(dftbplus_module.Dftb, "calculate", fake_calculate)


def test_geoopt_configures_dftb_calculator(monkeypatch, tmp_path):
    _stub_dftb_base(monkeypatch)
    atoms = Atoms("H2", positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 0.7]])

    geoopt = Geoopt(
        atoms,
        label="custom_geo",
        charge=-1,
        slako_dir=str(tmp_path),
        max_force=2.0e-5,
        Hamiltonian_SCC="Yes",
    )

    assert geoopt.atoms == atoms
    assert geoopt.calculated_atoms == atoms
    assert geoopt.label == "custom_geo"
    assert geoopt.slako_dir == str(tmp_path.resolve()) + os.sep
    assert geoopt.parameters["Hamiltonian_Charge"] == -1
    assert geoopt.parameters["Driver_"] == "GeometryOptimisation"
    assert geoopt.parameters["Driver_Optimiser"] == "Rational {}"
    assert geoopt.parameters["Driver_MaxForceComponent"] == 2.0e-5
    assert geoopt.parameters["Driver_OutputPrefix"] == "geo_opt"
    assert geoopt.parameters["Hamiltonian_SCC"] == "Yes"


def test_hessian_configures_dftb_calculator(monkeypatch, tmp_path):
    _stub_dftb_base(monkeypatch)
    atoms = Atoms("H2", positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 0.7]])

    hessian = Hessian(
        atoms,
        label="custom_hessian",
        charge=1,
        delta=0.002,
        slako_dir=str(tmp_path),
        Hamiltonian_SCC="Yes",
    )

    assert hessian.atoms == atoms
    assert hessian.calculated_atoms == atoms
    assert hessian.label == "custom_hessian"
    assert hessian.slako_dir == str(tmp_path.resolve()) + os.sep
    assert hessian.parameters["Hamiltonian_Charge"] == 1
    assert hessian.parameters["Driver_"] == "SecondDerivatives"
    assert hessian.parameters["Driver_Delta"] == 0.002
    assert hessian.parameters["Hamiltonian_SCC"] == "Yes"


def test_geoopt_potential_energy_converts_ev_to_hartree():
    class FakeAtoms:
        calc = None

        def get_potential_energy(self):
            return 2.5

    geoopt = Geoopt.__new__(Geoopt)
    geoopt.atoms = FakeAtoms()

    assert geoopt.potential_energy() == pytest.approx(
        2.5
        * dftbplus_module.PhysicalConstants["eV"]
        / dftbplus_module.PhysicalConstants["H"]
    )
    assert geoopt.atoms.calc is geoopt


def test_geoopt_read_uses_gen_reader(monkeypatch):
    expected_atoms = object()
    calls = []

    def fake_read(path, format=None):
        calls.append((path, format))
        return expected_atoms

    monkeypatch.setattr(dftbplus_module, "read", fake_read)
    geoopt = Geoopt.__new__(Geoopt)

    assert geoopt.read() is expected_atoms
    assert geoopt.atoms is expected_atoms
    assert calls == [("geo_opt.gen", "gen")]


def test_hessian_read_parses_square_grid(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    np.savetxt("hessian.out", np.arange(36, dtype=float).reshape(6, 6))

    class FakeAtoms:
        def get_global_number_of_atoms(self):
            return 2

    hessian = Hessian.__new__(Hessian)
    hessian.atoms = FakeAtoms()

    result = hessian.read()

    assert result.shape == (6, 6)
    np.testing.assert_array_equal(result, hessian.hessian)
    np.testing.assert_array_equal(result.ravel(), np.arange(36, dtype=float))


def test_hessian_read_parses_wrapped_dftbplus_output(monkeypatch, tmp_path):
    # DFTB+ writes hessian.out as a flat stream wrapped at a fixed number of
    # values per matrix row, producing ragged line widths (e.g. 4 then 2) that
    # numpy.loadtxt cannot parse. The reader must read all values and reshape.
    monkeypatch.chdir(tmp_path)
    expected = np.arange(36, dtype=float).reshape(6, 6)
    with open("hessian.out", "w", encoding="utf-8") as handle:
        for row in expected:
            for start in range(0, row.size, 4):
                handle.write(
                    "  ".join(f"{value:.10f}" for value in row[start:start + 4])
                    + "\n"
                )

    # Sanity-check the fixture is genuinely the ragged layout DFTB+ emits.
    line_widths = {
        len(line.split())
        for line in open("hessian.out", encoding="utf-8")
        if line.strip()
    }
    assert line_widths == {4, 2}
    with pytest.raises(ValueError):
        np.loadtxt("hessian.out")

    class FakeAtoms:
        def get_global_number_of_atoms(self):
            return 2

    hessian = Hessian.__new__(Hessian)
    hessian.atoms = FakeAtoms()

    result = hessian.read()

    assert result.shape == (6, 6)
    np.testing.assert_allclose(result, expected)


def test_hessian_read_rejects_wrong_matrix_size(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    # 2 atoms -> a 6x6 (36-value) Hessian is expected; write only 30 values.
    with open("hessian.out", "w", encoding="utf-8") as handle:
        handle.write("\n".join("0.0 0.0 0.0 0.0 0.0" for _ in range(6)) + "\n")

    class FakeAtoms:
        def get_global_number_of_atoms(self):
            return 2

    hessian = Hessian.__new__(Hessian)
    hessian.atoms = FakeAtoms()

    with pytest.raises(ValueError, match="Hessian matrix size"):
        hessian.read()


def test_modes_initialization_runs_steps(monkeypatch):
    calls = []

    monkeypatch.setattr(Modes, "write", lambda self: calls.append("write"))
    monkeypatch.setattr(Modes, "calculate", lambda self: calls.append("calculate"))
    monkeypatch.setattr(Modes, "read", lambda self: np.array([1.0, 2.0]))

    modes = Modes(geometry="optimized.gen", hessian="second.out")

    assert calls == ["write", "calculate"]
    assert modes.geometry == "optimized.gen"
    assert modes.hessian == "second.out"
    np.testing.assert_array_equal(modes.wave_numbers, np.array([1.0, 2.0]))


def test_modes_write_creates_input_file(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    modes = Modes.__new__(Modes)
    modes.geometry = "optimized.gen"
    modes.hessian = "second.out"

    assert modes.write() is None

    assert (
        tmp_path / "modes_in.hsd"
    ).read_text(encoding="utf-8") == (
        "Geometry = GenFormat {\n"
        "    <<< optimized.gen\n"
        "}\n"
        "\n"
        "Hessian = {\n"
        "    <<< second.out\n"
        "}\n"
        "\n"
        "Atoms = 1:-1\n"
        "\n"
    )


def test_modes_calculate_runs_modes_executable(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(dftbplus_module.shutil, "which", lambda command: "/usr/bin/modes")

    def fake_run(command, stdout, check):
        stdout.write("modes output")
        assert command == ["modes"]
        assert check is True

    monkeypatch.setattr(dftbplus_module.subprocess, "run", fake_run)
    modes = Modes.__new__(Modes)

    assert modes.calculate() is None
    assert (tmp_path / "modes.out").read_text(encoding="utf-8") == "modes output"


def test_modes_read_converts_hartree_to_wavenumbers(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "vibrations.tag").write_text(
        "header\n"
        "0.0 0.1\n"
        "0.2\n",
        encoding="utf-8",
    )
    modes = Modes.__new__(Modes)

    wave_numbers = modes.read()

    np.testing.assert_allclose(wave_numbers, np.array([0.0, 0.1, 0.2]) * 219474.63)
    np.testing.assert_array_equal(wave_numbers, modes.wave_numbers)


def test_modes_read_stops_at_trailing_tag_section(monkeypatch, tmp_path):
    # DFTB+ may append further tag sections (e.g. saved_modes) after the
    # frequencies; read() must stop there instead of parsing the header as a
    # float.
    monkeypatch.chdir(tmp_path)
    (tmp_path / "vibrations.tag").write_text(
        "frequencies         :real:1:3\n"
        " 0.100000000000000E-002  0.200000000000000E-002  0.300000000000000E-002\n"
        "saved_modes         :integer:1:3\n"
        "                                1                                2\n"
        "                                3\n",
        encoding="utf-8",
    )
    modes = Modes.__new__(Modes)

    wave_numbers = modes.read()

    assert wave_numbers.shape == (3,)
    np.testing.assert_allclose(
        wave_numbers, np.array([1.0e-3, 2.0e-3, 3.0e-3]) * 219474.63
    )


def test_modes_pipeline_parses_authentic_dftbplus_layout(monkeypatch, tmp_path):
    # End-to-end Hessian -> modes -> frequency path against the authentic DFTB+
    # file layouts, with only the (CI-unavailable) `modes` binary mocked.
    data_dir = Path(__file__).resolve().parents[1] / "data" / "calculator" / "modes"
    for name in ("geo_opt.gen", "hessian.out", "vibrations.tag"):
        (tmp_path / name).write_text(
            (data_dir / name).read_text(encoding="utf-8"), encoding="utf-8"
        )
    monkeypatch.chdir(tmp_path)

    # The committed hessian.out is the real wrapped/ragged layout that a
    # square-grid reader (np.loadtxt) cannot parse.
    with pytest.raises(ValueError):
        np.loadtxt("hessian.out")

    # Mock only the binary: it must run after write() created modes_in.hsd, and
    # it leaves vibrations.tag (already staged) as its output.
    monkeypatch.setattr(
        dftbplus_module.shutil, "which", lambda command: "/usr/bin/modes"
    )

    def fake_run(command, stdout, check):
        assert command == ["modes"]
        assert check is True
        assert Path("modes_in.hsd").exists()
        stdout.write("modes finished")

    monkeypatch.setattr(dftbplus_module.subprocess, "run", fake_run)

    modes = Modes(geometry="geo_opt.gen", hessian="hessian.out")

    expected = []
    with open("vibrations.tag", encoding="utf-8") as handle:
        handle.readline()
        for line in handle:
            expected += line.split()
    expected = np.array(expected, dtype=float) * 219474.63

    assert modes.wave_numbers.shape == (9,)
    np.testing.assert_allclose(modes.wave_numbers, expected)

    # The Hessian stage parses the same authentic ragged layout into (3N, 3N).
    class FakeAtoms:
        def get_global_number_of_atoms(self):
            return 3

    hessian = Hessian.__new__(Hessian)
    hessian.atoms = FakeAtoms()
    assert hessian.read().shape == (9, 9)


def test_spin_kwargs_restricted_and_fractional():
    h2 = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.74]])
    assert dftbplus_module._spin_kwargs(h2, None) == {}
    assert dftbplus_module._spin_kwargs(h2, 0.0) == {}
    assert dftbplus_module._spin_kwargs(h2, 0.1) == {}  # rounds to 0 unpaired


def test_spin_kwargs_builds_colinear_block():
    kw = dftbplus_module._spin_kwargs(
        Atoms("OH", positions=[[0, 0, 0], [0.97, 0, 0]]), 0.5
    )
    assert kw["Hamiltonian_SpinPolarisation"] == "Colinear {"
    assert kw["Hamiltonian_SpinPolarisation_UnpairedElectrons"] == 1
    assert kw["Hamiltonian_SpinConstants_ShellResolvedSpin"] == "No"
    assert kw["Hamiltonian_SpinConstants_H"] == "{ -0.07174 }"
    assert kw["Hamiltonian_SpinConstants_O"] == "{ -0.02785 }"


def test_spin_kwargs_rejects_element_without_constant():
    with pytest.raises(ValueError, match="not available for element"):
        dftbplus_module._spin_kwargs(
            Atoms("Br2", positions=[[0, 0, 0], [0, 0, 2.3]]), 0.5
        )


@pytest.mark.skipif(
    not dftbplus_ready,
    reason="DFTB+ executables are not installed or cannot start.",
)
class TestDftbplus:

    # Test the DFTB+ calculator
    def test_dftb(self):
        assert executable_starts("dftb+")
        assert executable_starts("modes")

    # Test the Geoopt class
    @pytest.mark.parametrize("example_dir", ["calculator"])
    def test_geoopt(self, test_with_data_dir):
        """
        Test the Geoopt class.
        """

        # read the atoms object
        atoms = ase_io.read("water.xyz")

        # create an instance of the Geoopt class
        geoopt = Geoopt(atoms, charge=0, **dftb_3ob_parameters)

        atoms.calc = geoopt

        assert geoopt.atoms == atoms
        assert np.allclose(geoopt.potential_energy(), -4.06231229, atol=1e-5)
        assert geoopt.label == "geo_opt"
        assert geoopt.slako_dir == _slako_dir()

        # read the optimized geometry
        geoopt.read()

        # check the optimized geometry
        assert geoopt.atoms == ase_io.read("geo_opt.gen", format="gen")
        assert atoms != geoopt.atoms
        assert atoms != ase_io.read("geo_opt.gen", format="gen")

    @pytest.mark.parametrize("example_dir", ["calculator"])
    def test_hessian(self, test_with_data_dir):
        """
        Test the Hessian class.
        """

        # read the atoms object
        atoms = ase_io.read("water.xyz")

        # create an instance of the Hessian class
        optimizer = Geoopt(atoms, charge=0, **dftb_3ob_parameters)
        atoms = optimizer.read()
        print(atoms.get_positions())
        hessian = Hessian(optimizer.read(),
                          charge=0,
                          **dftb_3ob_parameters,
                          delta=0.0005)

        assert hessian.atoms == optimizer.read()
        assert hessian.label == "second_derivative"
        assert hessian.slako_dir == _slako_dir()

        assert os.path.exists("hessian.out")
        assert np.allclose(hessian.read(), hessian.hessian, atol=1e-5)

    @pytest.mark.parametrize("example_dir", ["calculator"])
    def test_modes(self, test_with_data_dir):
        """
        Test the Modes class.
        """

        # read the atoms object
        atoms = ase_io.read("water.xyz")

        # create an instance of the Modes class
        optimizer = Geoopt(atoms, charge=0, **dftb_3ob_parameters)
        hessian = Hessian(optimizer.read(),
                          charge=0,
                          **dftb_3ob_parameters,
                          delta=0.0005)

        assert os.path.exists("hessian.out")
        assert os.path.exists("geo_opt.gen")

        modes = Modes()

        assert modes.geometry == "geo_opt.gen"
        assert modes.hessian == "hessian.out"

        assert os.path.exists("modes_in.hsd")
        assert os.path.exists("modes.out")
        assert os.path.exists("vibrations.tag")

        wave_numbers = modes.read()

        assert np.allclose(
            wave_numbers[6:],
            [1462.29964226, 3604.10862879, 3876.6277397],
            atol=1e-5,
        )

    @pytest.mark.parametrize("example_dir", ["calculator"])
    def test_dftbplus_thermo(self, test_with_data_dir):
        atoms = ase_io.read("water.xyz")
        thermo = dftbplus_thermo(atoms, **dftb_3ob_parameters, delta=0.0005)
        assert np.allclose(thermo.total_EeGtot(),
                           -4.0595598803365744,
                           atol=1e-5)
