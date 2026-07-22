import pytest

import unittest
import numpy as np
import pytest
from types import SimpleNamespace
from pathlib import Path
from unittest.mock import patch, mock_open
from beartype.roar import BeartypeCallHintParamViolation
from ThermoScreening.calculator.dftbplus import dftb_3ob_parameters
from ThermoScreening.thermo.api import (
    dftbplus_thermo,
    read_coord,
    read_vibrational,
    read_gen,
    read_xyz, 
    run_thermo,
    unit_length,
    unit_mass,
    unit_energy,
    unit_frequency,
)
from ase.build import molecule
from ase import Atoms
from ThermoScreening.exceptions import TSNotImplementedError, TSValueError

class TestApi(unittest.TestCase):
    
    def test_read_coord_not_implemented(self):
        with pytest.raises(TSNotImplementedError) as e:
            read_coord("test.xyz", engine="not_implemented")
        assert str(e.value) == "The engine is not supported."
        
        with pytest.raises(TSNotImplementedError) as e:
            read_coord("test.com", engine="dftb+")
        assert str(e.value) == "The input file is not supported."
        
    def test_read_vibrational_not_implemented(self):
        with pytest.raises(TSNotImplementedError) as e:
            read_vibrational("test.vib", engine="not_implemented")
        assert str(e.value) == "The engine is not supported."
        
    def test_unit_length(self):
        with pytest.raises(TSNotImplementedError) as e:
            unit_length(engine="not_implemented")
        assert str(e.value) == "The engine is not supported."
        
        assert unit_length(engine="dftb+") == "Angstrom"
        
    def test_unit_mass(self):
        with pytest.raises(TSNotImplementedError) as e:
            unit_mass(engine="not_implemented")
        assert str(e.value) == "The engine is not supported."
        
        assert unit_mass(engine="dftb+") == "amu"
        
    def test_unit_energy(self):
        with pytest.raises(TSNotImplementedError) as e:
            unit_energy(engine="not_implemented")
        assert str(e.value) == "The engine is not supported."
        
        assert unit_energy(engine="dftb+") == "Hartree"
        
    def test_unit_frequency(self):
        with pytest.raises(TSNotImplementedError) as e:
            unit_frequency(engine="not_implemented")
        assert str(e.value) == "The engine is not supported."
        
        assert unit_frequency(engine="dftb+") == "cm^-1"

    def test_public_api_rejects_wrong_argument_type(self):
        with pytest.raises(BeartypeCallHintParamViolation):
            unit_length(engine=object())
            
   
    @patch(
        "pathlib.Path.read_text",
        return_value="24\n\nO     0.00000003     -0.00000060     -2.14906255      6.45170211\n    O     -0.00000001      0.00000032     -7.56375986      6.45170117\n    C     -3.70943130     -0.00000041     -5.55541507      4.06573254\n    C     -2.50358224      0.00000055     -6.25344004      4.05635404\n    C     -1.28152402      0.00000064     -5.56329596      4.05319169\n    C     -1.28152226     -0.00000048     -4.14952765      4.05319299\n    C     -2.50358591     -0.00000218     -3.45938372      4.05635329\n    C     -3.70942822     -0.00000196     -4.15740715      4.06573042\n    C      0.00000004      0.00000235     -6.33321909      3.57768382\n    C     -0.00000003      0.00000103     -3.37960436      3.57768340\n    C      1.28152286      0.00000564     -4.14952745      4.05319251\n    C      1.28152346      0.00000677     -5.56329581      4.05319217\n    C      2.50358346      0.00001200     -6.25343979      4.05635377\n    H      2.48494730      0.00001319     -7.34030311      0.89492446\n    C      3.70943026      0.00001569     -5.55541540      4.06573185\n    C      3.70942928      0.00001410     -4.15740743      4.06573112\n    C      2.50358473      0.00000924     -3.45938351      4.05635352\n    H     -4.65151464      0.00000015     -6.09787561      0.91510467\n    H     -2.48494736      0.00000157     -7.34030301      0.89492435\n    H     -2.48494723     -0.00000334     -2.37251966      0.89492456\n    H     -4.65151524     -0.00000279     -3.61494721      0.91510621\n    H      4.65151483      0.00001973     -6.09787570      0.91510518\n    H      4.65151507      0.00001672     -3.61494731      0.91510569\n    H      2.48494732      0.00000821     -2.37251974      0.89492448",
    )
    def test_read_coord(self, read_text_mock):
        
        coord = np.array([[ 0.00000003  ,   -0.00000060  ,   -2.14906255],
                            [-0.00000001 ,    0.00000032  ,   -7.56375986],
                            [-3.70943130 ,   -0.00000041  ,   -5.55541507],
                            [-2.50358224 ,    0.00000055  ,   -6.25344004],
                            [-1.28152402 ,    0.00000064  ,   -5.56329596],
                            [-1.28152226 ,   -0.00000048  ,   -4.14952765],
                            [-2.50358591 ,   -0.00000218  ,   -3.45938372],
                            [-3.70942822 ,   -0.00000196  ,   -4.15740715],
                            [ 0.00000004 ,    0.00000235  ,   -6.33321909],
                            [-0.00000003 ,    0.00000103  ,   -3.37960436],
                            [ 1.28152286 ,    0.00000564  ,   -4.14952745],
                            [ 1.28152346 ,    0.00000677  ,   -5.56329581],
                            [ 2.50358346 ,    0.00001200  ,   -6.25343979],
                            [ 2.48494730 ,    0.00001319  ,   -7.34030311],
                            [ 3.70943026 ,    0.00001569  ,   -5.55541540],
                            [ 3.70942928 ,    0.00001410  ,   -4.15740743],
                            [ 2.50358473 ,    0.00000924  ,   -3.45938351],
                            [-4.65151464 ,    0.00000015  ,   -6.09787561],
                            [-2.48494736 ,    0.00000157  ,   -7.34030301],
                            [-2.48494723 ,   -0.00000334  ,   -2.37251966],
                            [-4.65151524 ,   -0.00000279  ,   -3.61494721],
                            [ 4.65151483 ,    0.00001973  ,   -6.09787570],
                            [ 4.65151507 ,    0.00001672  ,   -3.61494731],
                            [ 2.48494732 ,    0.00000821  ,   -2.37251974]])
        atom_number = 24
        atoms = np.array(
            [
                "O",
                "O",
                "C",
                "C",
                "C",
                "C",
                "C",
                "C",
                "C",
                "C",
                "C",
                "C",
                "C",
                "H",
                "C",
                "C",
                "C",
                "H",
                "H",
                "H",
                "H",
                "H",
                "H",
                "H",
            ]
        )

        cell = None

        data_N, data_atoms, data_xyz, cell, pbc = read_xyz("test.xyz")

        assert data_N == atom_number
        np.testing.assert_array_equal(data_atoms, atoms)
        np.testing.assert_array_almost_equal(data_xyz[:, 0], coord[:, 0], decimal=3)
        assert cell == None
        assert pbc == False

    @patch(
        "pathlib.Path.read_text",
        return_value=(
            "2 10.0 11.0 12.0 90.0 90.0 90.0\n"
            "\n"
            "Cl 0.0 0.0 0.0\n"
            "Na 1.0 2.0 3.0\n"
        ),
    )
    def test_read_xyz_keeps_multichar_symbols(self, read_text_mock):
        data_N, data_atoms, data_xyz, cell, pbc = read_xyz("test.xyz")

        assert data_N == 2
        np.testing.assert_array_equal(data_atoms, np.array(["Cl", "Na"], dtype=object))
        np.testing.assert_allclose(data_xyz, np.array([[0.0, 0.0, 0.0], [1.0, 2.0, 3.0]]))
        np.testing.assert_allclose(cell, np.array([10.0, 11.0, 12.0, 90.0, 90.0, 90.0]))
        assert pbc is True

    @patch(
        "pathlib.Path.read_text",
        return_value=(
            "1 bad 11.0 12.0 90.0 90.0 90.0\n"
            "\n"
            "Cl 0.0 0.0 0.0\n"
        ),
    )
    def test_read_xyz_rejects_invalid_cell_header(self, read_text_mock):
        with pytest.raises(TSValueError, match="Invalid XYZ coordinate file"):
            read_xyz("test.xyz")

    @patch(
        "pathlib.Path.read_text",
        return_value="1\n\nH 0.12345678901234 1.0 2.0\n",
    )
    def test_read_xyz_positions_keep_double_precision(self, read_text_mock):
        # XYZFrameReader stores positions as float32; read_xyz must preserve the
        # full double precision of the source coordinates.
        _, _, data_xyz, _, _ = read_xyz("test.xyz")

        assert data_xyz.dtype == np.float64
        np.testing.assert_allclose(
            data_xyz[0],
            np.array([0.12345678901234, 1.0, 2.0]),
            rtol=0.0,
            atol=1e-15,
        )

    def test_read_xyz_missing_file_raises_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            read_xyz("/nonexistent/path/does_not_exist.xyz")

    @patch(
        "builtins.open",
        new_callable=mock_open,
        read_data="1.0  2.0\n2.0  3.0\n3.0  4.0\n4.0  5.0\n5.0  6.0\n6.0  7.0\n",
    )
    def test_read_vib_file(self, mock_open):
        vibrational_frequencies = np.array([2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
        data_vib = read_vibrational("test.vib", "dftb+")
        np.testing.assert_array_almost_equal(
            data_vib, vibrational_frequencies, decimal=8
        )

    @patch("builtins.open", new_callable=mock_open, read_data="1.0\n")
    def test_read_vib_file_rejects_missing_frequency_column(self, mock_open):
        with pytest.raises(TSValueError, match="Invalid vibrational frequency line 1"):
            read_vibrational("test.vib", "dftb+")

    @patch("builtins.open", new_callable=mock_open, read_data="\n")
    def test_read_vib_file_rejects_empty_file(self, mock_open):
        with pytest.raises(TSValueError, match="No vibrational frequencies found"):
            read_vibrational("test.vib", "dftb+")

    def test_read_gen(self):
        gen_file = Path(__file__).resolve().parents[1] / "data/thermo/geo_opt.gen"

        data_N, data_atoms, data_xyz, cell, pbc = read_gen(str(gen_file))

        assert data_N == 24
        assert pbc is True
        np.testing.assert_array_equal(data_atoms[:4], np.array(["O", "O", "C", "C"]))
        np.testing.assert_allclose(
            data_xyz[0],
            np.array([3.0e-08, -6.0e-07, -2.14906255]),
            atol=1.0e-12,
        )
        np.testing.assert_allclose(
            cell,
            np.array(
                [
                    [10000.0, 0.0, 0.0],
                    [0.0, 10000.0, 0.0],
                    [0.0, 0.0, 10000.0],
                ]
            ),
        )

        coord_data = read_coord(str(gen_file), engine="dftb+")
        assert coord_data[0] == data_N
        np.testing.assert_array_equal(coord_data[1], data_atoms)
        np.testing.assert_allclose(coord_data[2], data_xyz)
        np.testing.assert_allclose(coord_data[3], cell)
        assert coord_data[4] is True

    @patch("ThermoScreening.thermo.api.read_gen_file")
    def test_read_gen_uses_pqanalysis(self, read_gen_file_mock):
        read_gen_file_mock.return_value = SimpleNamespace(
            n_atoms=2,
            atoms=[SimpleNamespace(name="Cl"), SimpleNamespace(name="Na")],
            pos=np.array([[0.0, 0.0, 0.0], [1.0, 2.0, 3.0]]),
            cell=SimpleNamespace(is_vacuum=True),
        )

        data_N, data_atoms, data_xyz, cell, pbc = read_gen("test.gen")

        read_gen_file_mock.assert_called_once_with("test.gen")
        assert data_N == 2
        np.testing.assert_array_equal(data_atoms, np.array(["Cl", "Na"], dtype=object))
        np.testing.assert_allclose(data_xyz, np.array([[0.0, 0.0, 0.0], [1.0, 2.0, 3.0]]))
        assert cell is None
        assert pbc is False

    def test_run_thermo_transition_state(self):
        # bent (non-linear) triatomic, dof=3: 6 near-zero trans/rot modes then
        # one imaginary (reaction-coordinate) mode and two real vibrations
        water = Atoms(
            "H2O",
            positions=[[0.0, 0.76, -0.48], [0.0, -0.76, -0.48], [0.0, 0.0, 0.12]],
        )
        frequencies = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -300.0, 1600.0, 3800.0])

        thermo = run_thermo(
            frequencies, atoms=water, engine="dftb+", energy=-76.0,
            transition_state=True,
        )

        assert thermo.imaginary_mode_wavenumber() == pytest.approx(-300.0)

    def test_run_thermo_non_transition_state_rejects_imaginary_frequency(self):
        water = Atoms(
            "H2O",
            positions=[[0.0, 0.76, -0.48], [0.0, -0.76, -0.48], [0.0, 0.0, 0.12]],
        )
        frequencies = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -300.0, 1600.0, 3800.0])

        with pytest.raises(TSValueError, match="geometry is not a minimum"):
            run_thermo(frequencies, atoms=water, engine="dftb+", energy=-76.0)

    def test_run_thermo_rejects_wrong_frequency_count(self):
        coord_file = Path(__file__).resolve().parents[1] / "data/thermo/geo_opt.xyz"

        with pytest.raises(TSValueError, match="number of vibrational frequencies"):
            run_thermo(
                np.array([1.0]),
                coord_file=str(coord_file),
                engine="dftb+",
            )

    def test_run_thermo_ase_atoms_matches_file_path(self):
        import ase.io as ase_io

        data_dir = Path(__file__).resolve().parents[1] / "data" / "thermo"
        gen_file = data_dir / "geo_opt.gen"
        frequencies = read_vibrational(str(data_dir / "frequency.txt"), "dftb+")

        from_file = run_thermo(
            frequencies, coord_file=str(gen_file), engine="dftb+", energy=-1.0
        )
        from_atoms = run_thermo(
            frequencies,
            atoms=ase_io.read(str(gen_file), format="gen"),
            engine="dftb+",
            energy=-1.0,
        )

        assert from_atoms.total_energy("H") == pytest.approx(from_file.total_energy("H"))
        assert from_atoms.total_entropy("cal/(mol*K)") == pytest.approx(
            from_file.total_entropy("cal/(mol*K)")
        )
        assert from_atoms.total_gibbs_free_energy("H") == pytest.approx(
            from_file.total_gibbs_free_energy("H")
        )


def test_run_in_directory_isolates_and_restores(tmp_path):
    from ThermoScreening.thermo.api import _run_in_directory

    start = Path.cwd()
    job = tmp_path / "job1"

    with _run_in_directory(str(job)):
        assert Path.cwd().resolve() == job.resolve()
        Path("artifact.txt").write_text("x", encoding="utf-8")

    assert Path.cwd() == start
    assert (job / "artifact.txt").exists()


def test_run_in_directory_restores_on_error(tmp_path):
    from ThermoScreening.thermo.api import _run_in_directory

    start = Path.cwd()

    with pytest.raises(RuntimeError):
        with _run_in_directory(str(tmp_path / "job2")):
            raise RuntimeError("boom")

    assert Path.cwd() == start


def test_run_in_directory_none_is_noop():
    from ThermoScreening.thermo.api import _run_in_directory

    start = Path.cwd()
    with _run_in_directory(None):
        assert Path.cwd() == start
    assert Path.cwd() == start


def test_dftbplus_thermo_runs_pipeline_in_directory(monkeypatch, tmp_path):
    # Drive dftbplus_thermo with the DFTB+ steps mocked, so the directory wiring
    # is exercised without the binaries.
    import ThermoScreening.thermo.api as api

    seen = {}

    class FakeGeoopt:
        def __init__(self, atoms, charge, **kwargs):
            seen["cwd_during"] = Path.cwd().resolve()
            seen["geoopt_kwargs"] = kwargs

        def potential_energy(self):
            return -1.0

        def read(self):
            return "optimized-atoms"

    class FakeHessian:
        def __init__(self, atoms, charge, **kwargs):
            seen["hessian_atoms"] = atoms

    class FakeModes:
        def __init__(self):
            self.wave_numbers = np.array([1.0, 2.0, 3.0])

    def fake_run_thermo(frequencies, atoms=None, **kwargs):
        seen["run_thermo_atoms"] = atoms
        return "thermo-result"

    monkeypatch.setattr(api, "Geoopt", FakeGeoopt)
    monkeypatch.setattr(api, "Hessian", FakeHessian)
    monkeypatch.setattr(api, "Modes", FakeModes)
    monkeypatch.setattr(api, "run_thermo", fake_run_thermo)

    job = tmp_path / "job"
    start = Path.cwd()
    water = Atoms("OH2", positions=[[0, 0, 0.12], [0, 0.76, -0.48], [0, -0.76, -0.48]])
    result = api.dftbplus_thermo(water, directory=str(job))

    assert result == "thermo-result"
    assert seen["run_thermo_atoms"] == "optimized-atoms"
    assert seen["hessian_atoms"] == "optimized-atoms"
    assert seen["cwd_during"] == job.resolve()  # pipeline ran inside the job dir
    assert Path.cwd() == start  # working directory restored afterwards
    # closed-shell water -> restricted, no spin polarisation in the DFTB+ kwargs
    assert "Hamiltonian_SpinPolarisation" not in seen["geoopt_kwargs"]


def _mock_pipeline(monkeypatch):
    import ThermoScreening.thermo.api as api

    seen = {}

    class FakeGeoopt:
        def __init__(self, atoms, charge, **kwargs):
            seen["geoopt_kwargs"] = kwargs

        def potential_energy(self):
            return -1.0

        def read(self):
            return "optimized-atoms"

    class FakeHessian:
        def __init__(self, atoms, charge, **kwargs):
            seen["hessian_kwargs"] = kwargs

    class FakeModes:
        def __init__(self):
            self.wave_numbers = np.array([1.0, 2.0, 3.0])

    def fake_run_thermo(frequencies, atoms=None, spin=None, quasi_rrho=False, **kwargs):
        seen["run_thermo_spin"] = spin
        seen["run_thermo_quasi_rrho"] = quasi_rrho
        return "thermo-result"

    monkeypatch.setattr(api, "Geoopt", FakeGeoopt)
    monkeypatch.setattr(api, "Hessian", FakeHessian)
    monkeypatch.setattr(api, "Modes", FakeModes)
    monkeypatch.setattr(api, "run_thermo", fake_run_thermo)
    return api, seen


def test_dftbplus_thermo_spin_polarises_radical_automatically(monkeypatch, tmp_path):
    api, seen = _mock_pipeline(monkeypatch)
    # OH: 9 electrons (odd) -> auto doublet -> spin-polarised, no user input
    api.dftbplus_thermo(
        Atoms("OH", positions=[[0, 0, 0], [0.97, 0, 0]]), directory=str(tmp_path / "j")
    )

    assert seen["run_thermo_spin"] == 0.5  # analysis multiplicity matches
    assert seen["geoopt_kwargs"]["Hamiltonian_SpinPolarisation"] == "Colinear {"
    assert seen["geoopt_kwargs"]["Hamiltonian_SpinPolarisation_UnpairedElectrons"] == 1
    assert seen["hessian_kwargs"]["Hamiltonian_SpinConstants_O"] == "{ -0.02785 }"


def test_dftbplus_thermo_restricted_for_closed_shell(monkeypatch, tmp_path):
    api, seen = _mock_pipeline(monkeypatch)
    # water: 10 electrons (even) -> singlet -> restricted, no spin kwargs
    api.dftbplus_thermo(
        Atoms("OH2", positions=[[0, 0, 0.12], [0, 0.76, -0.48], [0, -0.76, -0.48]]),
        directory=str(tmp_path / "j"),
    )

    assert seen["run_thermo_spin"] == 0.0
    assert "Hamiltonian_SpinPolarisation" not in seen["geoopt_kwargs"]


def test_dftbplus_thermo_explicit_triplet(monkeypatch, tmp_path):
    api, seen = _mock_pipeline(monkeypatch)
    # O2 is even-electron but a triplet -> user declares spin=1 -> 2 unpaired
    api.dftbplus_thermo(
        Atoms("O2", positions=[[0, 0, 0], [0, 0, 1.2]]),
        directory=str(tmp_path / "j"),
        spin=1.0,
    )

    assert seen["run_thermo_spin"] == 1.0
    assert seen["geoopt_kwargs"]["Hamiltonian_SpinPolarisation_UnpairedElectrons"] == 2


def test_dftbplus_thermo_uses_given_spin_constants(monkeypatch, tmp_path):
    api, seen = _mock_pipeline(monkeypatch)
    # a caller-supplied spin-constant table flows into both DFTB+ steps
    custom = {"H": "{ -0.088 }", "O": "{ -0.099 }"}
    api.dftbplus_thermo(
        Atoms("OH", positions=[[0, 0, 0], [0.97, 0, 0]]),
        directory=str(tmp_path / "j"),
        spin_constants=custom,
    )

    assert seen["geoopt_kwargs"]["Hamiltonian_SpinConstants_O"] == "{ -0.099 }"
    assert seen["hessian_kwargs"]["Hamiltonian_SpinConstants_O"] == "{ -0.099 }"


def test_dftbplus_thermo_injects_solvation_into_both_steps(monkeypatch, tmp_path):
    api, seen = _mock_pipeline(monkeypatch)
    param = tmp_path / "param_gbsa_h2o.txt"
    param.write_text("data", encoding="utf-8")

    api.dftbplus_thermo(
        Atoms("OH2", positions=[[0, 0, 0.12], [0, 0.76, -0.48], [0, -0.76, -0.48]]),
        directory=str(tmp_path / "j"),
        solvation_param_file=str(param),
    )

    # the optimisation and the Hessian both run in solvent (consistent geometry
    # and frequencies)
    assert seen["geoopt_kwargs"]["Hamiltonian_Solvation"] == "GeneralizedBorn {"
    assert seen["hessian_kwargs"]["Hamiltonian_Solvation"] == "GeneralizedBorn {"


def test_dftbplus_thermo_injects_dispersion_into_both_steps(monkeypatch, tmp_path):
    api, seen = _mock_pipeline(monkeypatch)

    api.dftbplus_thermo(
        Atoms("OH2", positions=[[0, 0, 0.12], [0, 0.76, -0.48], [0, -0.76, -0.48]]),
        directory=str(tmp_path / "j"),
        dispersion="d3-bj",
    )

    for step in ("geoopt_kwargs", "hessian_kwargs"):
        assert seen[step]["Hamiltonian_Dispersion"] == "DftD3 {"
        assert seen[step]["Hamiltonian_Dispersion_Damping_a1"] == 0.5719


def test_dftbplus_thermo_no_dispersion_by_default(monkeypatch, tmp_path):
    api, seen = _mock_pipeline(monkeypatch)
    api.dftbplus_thermo(
        Atoms("OH2", positions=[[0, 0, 0.12], [0, 0.76, -0.48], [0, -0.76, -0.48]]),
        directory=str(tmp_path / "j"),
    )

    assert "Hamiltonian_Dispersion" not in seen["geoopt_kwargs"]


def test_dftbplus_thermo_gas_phase_has_no_solvation(monkeypatch, tmp_path):
    api, seen = _mock_pipeline(monkeypatch)
    api.dftbplus_thermo(
        Atoms("OH2", positions=[[0, 0, 0.12], [0, 0.76, -0.48], [0, -0.76, -0.48]]),
        directory=str(tmp_path / "j"),
    )

    assert "Hamiltonian_Solvation" not in seen["geoopt_kwargs"]


def test_dftbplus_thermo_forwards_quasi_rrho(monkeypatch, tmp_path):
    api, seen = _mock_pipeline(monkeypatch)
    api.dftbplus_thermo(
        Atoms("OH2", positions=[[0, 0, 0.12], [0, 0.76, -0.48], [0, -0.76, -0.48]]),
        directory=str(tmp_path / "j"),
        quasi_rrho=True,
    )

    assert seen["run_thermo_quasi_rrho"] is True


def test_xtb_thermo_pipeline(monkeypatch, tmp_path):
    import ThermoScreening.thermo.api as api

    seen = {}

    def fake_calculator(method):
        seen["method"] = method
        return "xtb-calc"

    def fake_optimise(atoms, calc, fmax=0.01):
        seen["charge"] = atoms.get_initial_charges().sum()
        seen["unpaired"] = atoms.get_initial_magnetic_moments().sum()
        seen["calc"] = calc
        return "optimized-atoms", -5.0, np.array([1500.0, 3600.0, 3700.0])

    def fake_run_thermo(frequencies, atoms=None, engine=None, spin=None,
                        quasi_rrho=False, **kwargs):
        seen["engine"] = engine
        seen["spin"] = spin
        seen["quasi_rrho"] = quasi_rrho
        seen["energy"] = kwargs.get("energy")
        return "thermo-result"

    monkeypatch.setattr(api, "xtb_calculator", fake_calculator)
    monkeypatch.setattr(api, "optimise_and_frequencies", fake_optimise)
    monkeypatch.setattr(api, "run_thermo", fake_run_thermo)

    # OH radical: 9 electrons -> auto doublet -> 1 unpaired electron for xTB
    result = api.xtb_thermo(
        Atoms("OH", positions=[[0, 0, 0], [0, 0, 0.97]]),
        directory=str(tmp_path / "j"),
        method="GFN1-xTB",
        quasi_rrho=True,
    )

    assert result == "thermo-result"
    assert seen["engine"] == "xtb"
    assert seen["spin"] == 0.5           # auto electron-count guess
    assert seen["unpaired"] == pytest.approx(1.0)
    assert seen["charge"] == pytest.approx(0.0)
    assert seen["method"] == "GFN1-xTB"
    assert seen["quasi_rrho"] is True
    assert seen["energy"] == -5.0


def test_xtb_cli_thermo_pipeline(monkeypatch, tmp_path):
    import ThermoScreening.thermo.api as api

    seen = {}

    def fake_run_xtb(atoms, charge=0.0, unpaired=0, method="GFN2-xTB", solvent=None):
        seen.update(charge=charge, unpaired=unpaired, method=method, solvent=solvent)
        return "optimized-atoms", -4.5, np.array([1500.0, 3600.0, 3700.0])

    def fake_run_thermo(frequencies, atoms=None, engine=None, spin=None,
                        quasi_rrho=False, **kwargs):
        seen.update(engine=engine, spin=spin, quasi_rrho=quasi_rrho,
                    energy=kwargs.get("energy"))
        return "thermo-result"

    monkeypatch.setattr(api, "run_xtb", fake_run_xtb)
    monkeypatch.setattr(api, "run_thermo", fake_run_thermo)

    # neutral OH radical -> auto doublet -> 1 unpaired electron; water solvation
    result = api.xtb_cli_thermo(
        Atoms("OH", positions=[[0, 0, 0], [0, 0, 0.97]]),
        directory=str(tmp_path / "j"),
        solvent="water",
        method="GFN1-xTB",
        quasi_rrho=True,
    )

    assert result == "thermo-result"
    assert seen["engine"] == "xtb"
    assert seen["spin"] == 0.5        # auto electron-count guess (9 electrons)
    assert seen["unpaired"] == 1      # round(2*S) passed to xtb --uhf
    assert seen["charge"] == 0.0
    assert seen["solvent"] == "water"
    assert seen["method"] == "GFN1-xTB"
    assert seen["quasi_rrho"] is True
    assert seen["energy"] == -4.5


def test_xtb_fukui_indices_pipeline(monkeypatch, tmp_path):
    import ThermoScreening.thermo.api as api

    seen = {}

    def fake_run_xtb_fukui(atoms, charge=0.0, unpaired=0, method="GFN2-xTB", solvent=None):
        seen.update(charge=charge, unpaired=unpaired, method=method, solvent=solvent)
        return [("O", 0.5, 0.5, 0.5), ("H", 0.25, 0.25, 0.25), ("H", 0.25, 0.25, 0.25)]

    monkeypatch.setattr(api, "run_xtb_fukui", fake_run_xtb_fukui)

    # neutral OH radical -> auto doublet -> 1 unpaired electron
    result = api.xtb_fukui_indices(
        Atoms("OH", positions=[[0, 0, 0], [0, 0, 0.97]]),
        directory=str(tmp_path / "j"),
        solvent="water",
        method="GFN1-xTB",
    )

    assert result == [("O", 0.5, 0.5, 0.5), ("H", 0.25, 0.25, 0.25), ("H", 0.25, 0.25, 0.25)]
    assert seen["unpaired"] == 1  # auto electron-count guess -> doublet -> 1 unpaired
    assert seen["charge"] == 0.0
    assert seen["solvent"] == "water"
    assert seen["method"] == "GFN1-xTB"


def test_xtb_fukui_indices_explicit_spin(monkeypatch, tmp_path):
    import ThermoScreening.thermo.api as api

    seen = {}

    def fake_run_xtb_fukui(atoms, charge=0.0, unpaired=0, method="GFN2-xTB", solvent=None):
        seen["unpaired"] = unpaired
        return []

    monkeypatch.setattr(api, "run_xtb_fukui", fake_run_xtb_fukui)

    api.xtb_fukui_indices(Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.74]]), charge=-2.0, spin=1.0)

    assert seen["unpaired"] == 2


if __name__ == "__main__":
    unittest.main()
