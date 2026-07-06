import math
import os
import shutil

import pytest

from ThermoScreening.thermo.reactions import (
    SHE_ABSOLUTE_POTENTIAL,
    reaction_free_energy,
    reduction_potential,
)

_H_TO_EV = 27.211386245988034
_H_TO_KCAL = 627.5094740631  # 1 Hartree in kcal/mol


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


def test_public_api_exported():
    from ThermoScreening.thermo import reaction_free_energy as rfe
    from ThermoScreening.thermo import reduction_potential as rp
    from ThermoScreening.thermo import reactions

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
