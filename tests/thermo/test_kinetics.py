import math

import pytest

from ThermoScreening.thermo.kinetics import (
    eyring_rate_constant,
    wigner_tunneling_correction,
)

_KB_OVER_H = 1.380649e-23 / 6.62607015e-34  # kB/h, s^-1 K^-1


class _FakeThermo:
    """A stand-in exposing only what the kinetics helpers use."""

    def __init__(self, eegtot):
        self._eegtot = eegtot

    def total_EeGtot(self):
        return self._eegtot


def test_eyring_rate_constant_zero_barrier_is_kT_over_h():
    reactant = _FakeThermo(-10.0)
    ts = _FakeThermo(-10.0)  # dG-ddagger = 0
    temperature = 298.15

    k = eyring_rate_constant([reactant], ts, temperature=temperature)

    assert k == pytest.approx(_KB_OVER_H * temperature)


def test_eyring_rate_constant_known_barrier():
    # solve k = kB T / h * exp(-dG/RT) = 1 analytically for dG, independently
    # of eyring_rate_constant, and check it reproduces k = 1 s^-1
    R, N_A, H = 8.314462618, 6.02214076e23, 4.359744722207101e-18
    temperature = 298.15
    dG_j_per_mol = R * temperature * math.log(_KB_OVER_H * temperature)
    dG_hartree = dG_j_per_mol / (H * N_A)

    reactant = _FakeThermo(-10.0)
    ts = _FakeThermo(-10.0 + dG_hartree)

    k = eyring_rate_constant([reactant], ts, temperature=temperature)

    assert k == pytest.approx(1.0)


def test_eyring_rate_constant_stoichiometric_reactants():
    # bimolecular convention: reactants is a list of Thermo/(coeff, Thermo);
    # only the total G(reactants) matters, split however
    ts = _FakeThermo(-9.0)
    single = [_FakeThermo(-10.0)]
    split = [_FakeThermo(-6.0), _FakeThermo(-4.0)]

    k_single = eyring_rate_constant(single, ts, temperature=300.0)
    k_split = eyring_rate_constant(split, ts, temperature=300.0)

    assert k_single == pytest.approx(k_split)


def test_eyring_rate_constant_stoichiometry_tuple():
    ts = _FakeThermo(-9.0)
    reactant = _FakeThermo(-5.0)
    # 2 A -> TS, same total G as two separate -5.0 reactants
    k_tuple = eyring_rate_constant([(2.0, reactant)], ts, temperature=300.0)
    k_split = eyring_rate_constant(
        [_FakeThermo(-5.0), _FakeThermo(-5.0)], ts, temperature=300.0
    )

    assert k_tuple == pytest.approx(k_split)


def test_eyring_rate_constant_kappa_scales_linearly():
    reactant = _FakeThermo(-10.0)
    ts = _FakeThermo(-9.9)

    k1 = eyring_rate_constant([reactant], ts, temperature=298.15, kappa=1.0)
    k2 = eyring_rate_constant([reactant], ts, temperature=298.15, kappa=2.5)

    assert k2 == pytest.approx(2.5 * k1)


def test_eyring_rate_constant_higher_barrier_is_slower():
    reactant = _FakeThermo(-10.0)
    low_ts = _FakeThermo(-9.98)
    high_ts = _FakeThermo(-9.90)

    k_low = eyring_rate_constant([reactant], low_ts, temperature=298.15)
    k_high = eyring_rate_constant([reactant], high_ts, temperature=298.15)

    assert k_low > k_high > 0


def test_wigner_tunneling_correction_small_mode_near_one():
    kappa = wigner_tunneling_correction(-50.0, temperature=298.15)
    assert kappa == pytest.approx(1.0, abs=0.01)
    assert kappa > 1.0


def test_wigner_tunneling_correction_ignores_sign():
    assert wigner_tunneling_correction(-500.0, 298.15) == pytest.approx(
        wigner_tunneling_correction(500.0, 298.15)
    )


def test_wigner_tunneling_correction_increases_with_wavenumber():
    small = wigner_tunneling_correction(-100.0, 298.15)
    large = wigner_tunneling_correction(-1000.0, 298.15)
    assert large > small > 1.0


def test_wigner_tunneling_correction_matches_closed_form():
    # kappa = 1 + (1/24) * (h c |nu| / (kB T))^2, computed independently here
    h, c, kB = 6.62607015e-34, 299792458.0, 1.380649e-23
    nu, temperature = 1200.0, 310.0
    u = h * c * nu * 100.0 / (kB * temperature)
    expected = 1.0 + u**2 / 24.0

    assert wigner_tunneling_correction(-nu, temperature) == pytest.approx(expected)


def test_public_api_exported():
    from ThermoScreening.thermo import eyring_rate_constant as erc
    from ThermoScreening.thermo import wigner_tunneling_correction as wtc
    from ThermoScreening.thermo import kinetics

    assert erc is kinetics.eyring_rate_constant
    assert wtc is kinetics.wigner_tunneling_correction
