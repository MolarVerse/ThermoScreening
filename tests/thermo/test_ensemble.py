import math

import pytest

from ThermoScreening.thermo.ensemble import (
    boltzmann_weights,
    ensemble_free_energy,
    lowest_gibbs,
)

_KB_HARTREE_PER_K = 1.380649e-23 / 4.3597447222071e-18  # k_B in Hartree/K
_H_TO_KCAL = 627.5094740631


class _FakeThermo:
    """A stand-in exposing only the free energy the ensemble helpers use."""

    def __init__(self, eegtot):
        self._eegtot = eegtot

    def total_EeGtot(self):
        return self._eegtot


def test_boltzmann_weights_degenerate():
    weights = boltzmann_weights([_FakeThermo(-1.0), _FakeThermo(-1.0)])
    assert weights == pytest.approx([0.5, 0.5])


def test_boltzmann_weights_single_conformer():
    assert boltzmann_weights([_FakeThermo(-3.14)]) == pytest.approx([1.0])


def test_boltzmann_weights_normalised_and_ordered():
    weights = boltzmann_weights([_FakeThermo(-1.0), _FakeThermo(-1.001), _FakeThermo(-0.999)])
    assert sum(weights) == pytest.approx(1.0)
    # lower free energy -> larger population
    assert weights[1] > weights[0] > weights[2]


def test_boltzmann_weights_known_ratio():
    # gap of exactly k_B*T -> population ratio of e:1
    temperature = 300.0
    gap = _KB_HARTREE_PER_K * temperature
    weights = boltzmann_weights([_FakeThermo(0.0), _FakeThermo(gap)], temperature=temperature)
    assert weights[0] / weights[1] == pytest.approx(math.e)


def test_ensemble_free_energy_single_is_identity():
    assert ensemble_free_energy([_FakeThermo(-2.0)], unit="H") == pytest.approx(-2.0)


def test_ensemble_free_energy_degenerate_mixing_entropy():
    # two degenerate conformers -> G_ens = e - k_B*T*ln(2), below either one
    temperature = 298.15
    energy = -5.0
    kt = _KB_HARTREE_PER_K * temperature
    expected = energy - kt * math.log(2)
    result = ensemble_free_energy([_FakeThermo(energy), _FakeThermo(energy)], temperature=temperature)
    assert result == pytest.approx(expected)
    assert result < energy


def test_ensemble_free_energy_at_or_below_minimum():
    thermos = [_FakeThermo(-1.0), _FakeThermo(-1.2), _FakeThermo(-0.8)]
    g_ens = ensemble_free_energy(thermos, unit="H")
    assert g_ens <= min(t.total_EeGtot() for t in thermos)


def test_ensemble_free_energy_unit_conversion():
    h = ensemble_free_energy([_FakeThermo(-1.0), _FakeThermo(-1.2)], unit="H")
    kcal = ensemble_free_energy([_FakeThermo(-1.0), _FakeThermo(-1.2)], unit="kcal")
    assert kcal == pytest.approx(h * _H_TO_KCAL, rel=1e-4)


def test_ensemble_free_energy_rejects_unknown_unit():
    with pytest.raises(ValueError, match="Unknown unit"):
        ensemble_free_energy([_FakeThermo(0.0)], unit="furlong")


def test_lowest_gibbs_returns_min_conformer():
    a, b, c = _FakeThermo(-1.0), _FakeThermo(-1.5), _FakeThermo(-0.5)
    assert lowest_gibbs([a, b, c]) is b


@pytest.mark.parametrize("func", [boltzmann_weights, ensemble_free_energy, lowest_gibbs])
def test_empty_ensemble_raises(func):
    with pytest.raises(ValueError, match="at least one conformer"):
        func([])


@pytest.mark.parametrize("func", [boltzmann_weights, ensemble_free_energy])
def test_non_positive_temperature_raises(func):
    with pytest.raises(ValueError, match="temperature must be positive"):
        func([_FakeThermo(0.0)], temperature=0.0)


def test_public_api_exported():
    from ThermoScreening.thermo import boltzmann_weights as bw
    from ThermoScreening.thermo import ensemble_free_energy as efe
    from ThermoScreening.thermo import lowest_gibbs as lg
    from ThermoScreening.thermo import ensemble

    assert bw is ensemble.boltzmann_weights
    assert efe is ensemble.ensemble_free_energy
    assert lg is ensemble.lowest_gibbs
