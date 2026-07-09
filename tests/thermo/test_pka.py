import math

import pytest

from ThermoScreening.thermo.pka import (
    pKa,
    calibrate_proton_reference,
    PROTON_AQUEOUS_FREE_ENERGY_KCAL,
)

_H_TO_KCAL = 627.5094740631  # 1 Hartree in kcal/mol
_R_KCAL = 8.314462618 / 4184.0  # kcal/(mol K)


class _FakeThermo:
    """A stand-in exposing only what the pKa helpers use."""

    def __init__(self, eegtot):
        self._eegtot = eegtot

    def total_EeGtot(self):
        return self._eegtot


def test_default_proton_reference_matches_the_documented_derivation():
    # G(H+, gas, 1 atm) = -6.28 kcal/mol (Bartmess 1994) + dG_solv(H+, aq,
    # 1 atm) = -264.0 kcal/mol (Tissandier 1998 / Kelly-Cramer-Truhlar 2006,
    # converted from their 1 mol/L reference state)
    assert PROTON_AQUEOUS_FREE_ENERGY_KCAL == pytest.approx(-270.28, abs=0.01)


def test_pKa_zero_delta_g_gives_pKa_from_reference_alone():
    acid = _FakeThermo(-10.0)
    base = _FakeThermo(-10.0)  # G(base) - G(acid) = 0

    result = pKa(acid, base, reference_free_energy=0.0)

    assert result == pytest.approx(0.0)


def test_pKa_known_value_acetic_acid_scale():
    # solve dG (Hartree) so that, with reference_free_energy=0, pKa lands at
    # the experimental acetic acid value (4.756) -- a numerical sign/scale
    # sanity check: a weak acid must give a POSITIVE pKa (unfavorable
    # dissociation), not a negative one
    temperature = 298.15
    target_pKa = 4.756
    dG_kcal = target_pKa * _R_KCAL * temperature * math.log(10)
    dG_hartree = dG_kcal / _H_TO_KCAL

    acid = _FakeThermo(-10.0)
    base = _FakeThermo(-10.0 + dG_hartree)

    result = pKa(acid, base, temperature=temperature, reference_free_energy=0.0)

    assert result == pytest.approx(target_pKa)
    assert result > 0  # weak acid -> positive pKa


def test_pKa_stronger_acid_has_lower_pKa():
    acid = _FakeThermo(-10.0)
    weak_base = _FakeThermo(-9.98)    # larger dG -> weaker acid -> higher pKa
    strong_base = _FakeThermo(-10.05)  # smaller (negative) dG -> stronger acid

    weak_pKa = pKa(acid, weak_base, reference_free_energy=0.0)
    strong_pKa = pKa(acid, strong_base, reference_free_energy=0.0)

    assert strong_pKa < weak_pKa


def test_pKa_reference_free_energy_shifts_result_linearly():
    acid, base = _FakeThermo(-10.0), _FakeThermo(-10.0)
    rt_ln10 = _R_KCAL * 298.15 * math.log(10)

    p0 = pKa(acid, base, reference_free_energy=0.0)
    p1 = pKa(acid, base, reference_free_energy=rt_ln10)  # +1 pKa unit worth

    assert p1 - p0 == pytest.approx(1.0)


def test_pKa_uses_default_reference_free_energy():
    acid, base = _FakeThermo(-10.0), _FakeThermo(-10.0)

    default_call = pKa(acid, base)
    explicit_call = pKa(acid, base, reference_free_energy=PROTON_AQUEOUS_FREE_ENERGY_KCAL)

    assert default_call == pytest.approx(explicit_call)


def test_calibrate_proton_reference_round_trips_with_pKa():
    acid = _FakeThermo(-10.0)
    base = _FakeThermo(-9.99)
    target_pKa = 7.2

    calibrated = calibrate_proton_reference(acid, base, target_pKa, temperature=310.0)
    recovered = pKa(acid, base, temperature=310.0, reference_free_energy=calibrated)

    assert recovered == pytest.approx(target_pKa)


def test_calibrate_proton_reference_matches_closed_form():
    acid, base = _FakeThermo(-10.0), _FakeThermo(-9.995)
    temperature = 298.15

    calibrated = calibrate_proton_reference(acid, base, 5.0, temperature=temperature)

    dG_kcal = (base.total_EeGtot() - acid.total_EeGtot()) * _H_TO_KCAL
    expected = 5.0 * _R_KCAL * temperature * math.log(10) - dG_kcal
    assert calibrated == pytest.approx(expected)


def test_public_api_exported():
    from ThermoScreening.thermo import pKa as pKa_pub
    from ThermoScreening.thermo import calibrate_proton_reference as cal_pub
    from ThermoScreening.thermo import pka as pka_module

    assert pKa_pub is pka_module.pKa
    assert cal_pub is pka_module.calibrate_proton_reference
