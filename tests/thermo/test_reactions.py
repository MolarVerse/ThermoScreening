import math
import os
import shutil

import pytest

from ThermoScreening.thermo.reactions import (
    SHE_ABSOLUTE_POTENTIAL,
    calibrate_reduction_reference,
    reaction_free_energy,
    reduction_potential,
)

_H_TO_EV = 27.211386245988034
_H_TO_KCAL = 627.5094740631  # 1 Hartree in kcal/mol
_H_TO_KJ = 2625.4996394798254  # 1 Hartree in kJ/mol


class _FakeThermo:
    """A stand-in exposing only what the reactions helpers use."""

    def __init__(self, eegtot):
        self._eegtot = eegtot

    def total_EeGtot(self):
        return self._eegtot


def test_reaction_free_energy_units():
    a, b, c = _FakeThermo(-1.0), _FakeThermo(-2.0), _FakeThermo(-3.5)
    # A + B -> C : dG = -3.5 - (-3.0) = -0.5 Hartree
    assert reaction_free_energy([a, b], [c], unit="H") == pytest.approx(-0.5)
    assert reaction_free_energy([a, b], [c], unit="eV") == pytest.approx(-0.5 * _H_TO_EV)
    assert reaction_free_energy([a, b], [c], unit="kcal") == pytest.approx(
        -0.5 * _H_TO_KCAL, rel=1e-4
    )
    assert reaction_free_energy([a, b], [c], unit="kJ") == pytest.approx(
        -0.5 * _H_TO_KJ, rel=1e-4
    )


def test_reaction_free_energy_stoichiometry():
    a, b = _FakeThermo(-1.0), _FakeThermo(-2.5)
    # 2 A -> B : dG = -2.5 - 2*(-1.0) = -0.5 Hartree
    assert reaction_free_energy([(2.0, a)], [b], unit="H") == pytest.approx(-0.5)


def test_reaction_free_energy_rejects_unknown_unit():
    with pytest.raises(ValueError, match="Unknown unit"):
        reaction_free_energy([_FakeThermo(0.0)], [_FakeThermo(0.0)], unit="furlong")


def test_reduction_potential_absolute_and_vs_reference():
    ox, red = _FakeThermo(0.0), _FakeThermo(-0.5)
    # dG_red = -0.5 Hartree -> E_abs = 0.5 * 27.2114 V
    assert reduction_potential(ox, red, reference_potential=0.0) == pytest.approx(
        0.5 * _H_TO_EV
    )
    assert reduction_potential(ox, red) == pytest.approx(
        0.5 * _H_TO_EV - SHE_ABSOLUTE_POTENTIAL
    )


def test_reduction_potential_n_electrons():
    ox, red = _FakeThermo(0.0), _FakeThermo(-1.0)
    # n = 2 : E_abs = 1.0 * 27.2114 / 2
    assert reduction_potential(ox, red, n_electrons=2, reference_potential=0.0) == pytest.approx(
        _H_TO_EV / 2
    )


@pytest.mark.parametrize("n_electrons", [0, -1])
def test_reduction_potential_rejects_non_positive_electrons(n_electrons):
    ox, red = _FakeThermo(0.0), _FakeThermo(-0.5)
    with pytest.raises(ValueError, match="positive"):
        reduction_potential(ox, red, n_electrons=n_electrons)


@pytest.mark.parametrize("n_electrons", [1, 2])
def test_calibrate_reduction_reference_reproduces_experiment(n_electrons):
    ox, red = _FakeThermo(0.0), _FakeThermo(-0.2)
    experimental = -0.75

    reference = calibrate_reduction_reference(
        ox, red, experimental, n_electrons=n_electrons
    )

    assert reduction_potential(
        ox,
        red,
        n_electrons=n_electrons,
        reference_potential=reference,
    ) == pytest.approx(experimental)


def test_calibrated_reference_transfers_relative_potential():
    reference_ox, reference_red = _FakeThermo(0.0), _FakeThermo(-0.2)
    target_ox, target_red = _FakeThermo(-1.0), _FakeThermo(-1.3)
    experimental = -0.75

    reference = calibrate_reduction_reference(
        reference_ox, reference_red, experimental
    )
    target = reduction_potential(
        target_ox, target_red, reference_potential=reference
    )

    assert target == pytest.approx(experimental + 0.1 * _H_TO_EV)


def test_two_step_reduction_uses_one_calibration_per_step():
    ref_neutral = _FakeThermo(-100.0)
    ref_anion = _FakeThermo(-100.13)
    ref_dianion = _FakeThermo(-100.23)
    target_neutral = _FakeThermo(-200.0)
    target_anion = _FakeThermo(-200.14)
    target_dianion = _FakeThermo(-200.25)
    experimental_e1 = -0.75
    experimental_e2 = -1.4

    e1_reference = calibrate_reduction_reference(
        ref_neutral, ref_anion, experimental_e1
    )
    e2_reference = calibrate_reduction_reference(
        ref_anion, ref_dianion, experimental_e2
    )
    e2e_reference = calibrate_reduction_reference(
        ref_neutral,
        ref_dianion,
        (experimental_e1 + experimental_e2) / 2,
        n_electrons=2,
    )

    target_e1 = reduction_potential(
        target_neutral, target_anion, reference_potential=e1_reference
    )
    target_e2 = reduction_potential(
        target_anion, target_dianion, reference_potential=e2_reference
    )
    target_e2e = reduction_potential(
        target_neutral,
        target_dianion,
        n_electrons=2,
        reference_potential=e2e_reference,
    )

    assert target_e1 == pytest.approx(experimental_e1 + 0.01 * _H_TO_EV)
    assert target_e2 == pytest.approx(experimental_e2 + 0.01 * _H_TO_EV)
    assert e2e_reference == pytest.approx((e1_reference + e2_reference) / 2)
    assert target_e2e == pytest.approx((target_e1 + target_e2) / 2)


def test_overall_two_electron_potential_is_mean_of_stepwise_potentials():
    neutral = _FakeThermo(-100.0)
    radical_anion = _FakeThermo(-100.13)
    dianion = _FakeThermo(-100.23)
    reference = 4.44

    first = reduction_potential(
        neutral, radical_anion, reference_potential=reference
    )
    second = reduction_potential(
        radical_anion, dianion, reference_potential=reference
    )
    overall = reduction_potential(
        neutral, dianion, n_electrons=2, reference_potential=reference
    )

    assert overall == pytest.approx((first + second) / 2)


@pytest.mark.parametrize("n_electrons", [0, -1])
def test_calibrate_reduction_reference_rejects_non_positive_electrons(
    n_electrons,
):
    ox, red = _FakeThermo(0.0), _FakeThermo(-0.5)
    with pytest.raises(ValueError, match="positive"):
        calibrate_reduction_reference(
            ox, red, -0.75, n_electrons=n_electrons
        )


def test_public_api_exported():
    from ThermoScreening.thermo import calibrate_reduction_reference as crr
    from ThermoScreening.thermo import reaction_free_energy as rfe
    from ThermoScreening.thermo import reduction_potential as rp
    from ThermoScreening.thermo import reactions

    assert crr is reactions.calibrate_reduction_reference
    assert rfe is reactions.reaction_free_energy
    assert rp is reactions.reduction_potential


xtb_available = shutil.which("xtb") is not None or "XTB_COMMAND" in os.environ


@pytest.mark.skipif(not xtb_available, reason="the native xtb binary is not available.")
def test_reduction_potential_end_to_end(tmp_path):
    # benzoquinone + e- -> radical anion, in water, via xtb-cli. GFN2 absolute
    # redox potentials are quantitatively poor, so this only checks that the
    # pipeline runs and the reduction is (correctly) favorable for a good acceptor.
    from ThermoScreening.thermo.api import xtb_cli_thermo
    from ThermoScreening.thermo.conformers import generate

    bq = generate("O=C1C=CC(=O)C=C1", max_conformers=3)[0]
    neutral = xtb_cli_thermo(bq, charge=0, solvent="water", directory=str(tmp_path / "n"))
    anion = xtb_cli_thermo(bq, charge=-1, solvent="water", directory=str(tmp_path / "a"))

    e_abs = reduction_potential(neutral, anion, reference_potential=0.0)
    assert math.isfinite(e_abs)
    assert e_abs > 0  # reduction of a quinone is favorable
