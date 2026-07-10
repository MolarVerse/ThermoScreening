import math
import os
import shutil

import pytest

from ThermoScreening.thermo.ensemble import (
    boltzmann_weights,
    ensemble_free_energy,
    lowest_gibbs,
    EnsembleThermo,
)

xtb_available = shutil.which("xtb") is not None or "XTB_COMMAND" in os.environ

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
@pytest.mark.parametrize("temperature", [0.0, -10.0])
def test_non_positive_temperature_raises(func, temperature):
    with pytest.raises(ValueError, match="temperature must be positive"):
        func([_FakeThermo(0.0)], temperature=temperature)


def test_ensemble_thermo_matches_ensemble_free_energy():
    thermos = [_FakeThermo(-1.0), _FakeThermo(-1.2), _FakeThermo(-0.8)]

    wrapped = EnsembleThermo(thermos)

    assert wrapped.total_EeGtot() == pytest.approx(ensemble_free_energy(thermos, unit="H"))


def test_ensemble_thermo_single_conformer_is_identity():
    assert EnsembleThermo([_FakeThermo(-2.0)]).total_EeGtot() == pytest.approx(-2.0)


def test_ensemble_thermo_forwards_temperature():
    thermos = [_FakeThermo(-1.0), _FakeThermo(-1.2)]

    wrapped = EnsembleThermo(thermos, temperature=350.0)

    assert wrapped.total_EeGtot() == pytest.approx(
        ensemble_free_energy(thermos, temperature=350.0, unit="H")
    )


def test_ensemble_thermo_empty_raises():
    with pytest.raises(ValueError, match="at least one conformer"):
        EnsembleThermo([])


def test_ensemble_thermo_drops_into_pKa():
    from ThermoScreening.thermo.pka import pKa

    acid_thermos = [_FakeThermo(-10.0), _FakeThermo(-9.999)]
    base_thermos = [_FakeThermo(-10.0 + 0.01), _FakeThermo(-9.995)]

    ensemble_result = pKa(
        EnsembleThermo(acid_thermos), EnsembleThermo(base_thermos), reference_free_energy=0.0
    )
    manual_result = pKa(
        _FakeThermo(ensemble_free_energy(acid_thermos, unit="H")),
        _FakeThermo(ensemble_free_energy(base_thermos, unit="H")),
        reference_free_energy=0.0,
    )

    assert ensemble_result == pytest.approx(manual_result)


def test_ensemble_thermo_drops_into_reduction_potential():
    from ThermoScreening.thermo.reactions import reduction_potential

    oxidized_thermos = [_FakeThermo(-10.0), _FakeThermo(-9.999)]
    reduced_thermos = [_FakeThermo(-10.01), _FakeThermo(-10.005)]

    ensemble_result = reduction_potential(
        EnsembleThermo(oxidized_thermos), EnsembleThermo(reduced_thermos)
    )
    manual_result = reduction_potential(
        _FakeThermo(ensemble_free_energy(oxidized_thermos, unit="H")),
        _FakeThermo(ensemble_free_energy(reduced_thermos, unit="H")),
    )

    assert ensemble_result == pytest.approx(manual_result)


def test_public_api_exported():
    from ThermoScreening.thermo import boltzmann_weights as bw
    from ThermoScreening.thermo import ensemble_free_energy as efe
    from ThermoScreening.thermo import lowest_gibbs as lg
    from ThermoScreening.thermo import EnsembleThermo as et
    from ThermoScreening.thermo import ensemble

    assert bw is ensemble.boltzmann_weights
    assert efe is ensemble.ensemble_free_energy
    assert lg is ensemble.lowest_gibbs
    assert et is ensemble.EnsembleThermo


@pytest.mark.skipif(not xtb_available, reason="the native xtb binary is not available.")
def test_ensemble_thermo_end_to_end_with_pKa(tmp_path):
    # 4-hydroxybutanoic acid has a flexible C-C-C-C backbone (unlike e.g.
    # glycolic acid, whose conformers RDKit's RMSD pruning collapses to one)
    # -- generate() gives it several genuinely distinct conformers, so this
    # exercises real Boltzmann averaging, not a degenerate single-conformer
    # ensemble. EnsembleThermo wraps the real xtb-cli energies over them and
    # duck-types as a Thermo, so it drops straight into pKa() with no other
    # changes.
    from ThermoScreening.thermo.api import xtb_cli_thermo
    from ThermoScreening.thermo.conformers import generate
    from ThermoScreening.thermo.pka import pKa

    def ensemble_thermo(smiles, charge, tag):
        conformers = generate(smiles, max_conformers=5)
        assert len(conformers) > 1  # otherwise this isn't testing an ensemble
        thermos = [
            xtb_cli_thermo(
                conf, charge=charge, solvent="water", directory=str(tmp_path / f"{tag}_{i}")
            )
            for i, conf in enumerate(conformers)
        ]
        return thermos, EnsembleThermo(thermos)

    acid_thermos, acid_ensemble = ensemble_thermo("OCCCC(=O)O", 0, "acid")
    base_thermos, base_ensemble = ensemble_thermo("OCCCC(=O)[O-]", -1, "base")

    # the ensemble free energy is at or below every individual conformer's
    assert acid_ensemble.total_EeGtot() <= min(t.total_EeGtot() for t in acid_thermos)
    assert base_ensemble.total_EeGtot() <= min(t.total_EeGtot() for t in base_thermos)

    ensemble_pKa = pKa(acid_ensemble, base_ensemble)

    assert math.isfinite(ensemble_pKa)
