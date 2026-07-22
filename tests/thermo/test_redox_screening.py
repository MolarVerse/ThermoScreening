import numpy as np
import pytest

from ThermoScreening.exceptions import TSValueError
from ThermoScreening.thermo.redox_screening import (
    DeltaRedoxModel,
    audit_redox_dataset,
    balanced_group_folds,
    canonical_structure_identity,
    grouped_delta_validation,
    pareto_front,
)


def _matched_series():
    approximate = np.array([-1.15, -0.82, -1.31, -0.73, -1.02, -0.91, -1.24, -0.67])
    groups = [
        {"amino": 1},
        {"nitro": 1},
        {"amino": 1, "hydroxy": 1},
        {"nitro": 1, "hydroxy": 1},
        {"hydroxy": 1},
        {"amino": 1, "nitro": 1},
        {"amino": 2},
        {"nitro": 2},
    ]
    correction = np.array(
        [
            0.04 + 0.12 * value + 0.03 * row.get("amino", 0) - 0.05 * row.get("nitro", 0)
            + 0.02 * row.get("hydroxy", 0)
            for value, row in zip(approximate, groups)
        ]
    )
    return approximate, approximate + correction, groups


def test_canonical_structure_identity_matches_equivalent_smiles():
    first = canonical_structure_identity("CCO")
    second = canonical_structure_identity("OCC")

    assert first.canonical_smiles == "CCO"
    assert first == second


def test_canonical_structure_identity_preserves_stereochemistry():
    clockwise = canonical_structure_identity("F[C@H](Cl)Br")
    anticlockwise = canonical_structure_identity("F[C@@H](Cl)Br")

    assert clockwise.structure_id != anticlockwise.structure_id


def test_canonical_structure_identity_distinguishes_fixed_h_tautomers():
    pyridone = canonical_structure_identity("O=c1cc[nH]cc1")
    hydroxypyridine = canonical_structure_identity("Oc1ccncc1")

    assert pyridone.structure_id != hydroxypyridine.structure_id


@pytest.mark.parametrize("smiles", ["", None, "not a smiles"])
def test_canonical_structure_identity_rejects_invalid_smiles(smiles):
    with pytest.raises(TSValueError, match="SMILES|Invalid"):
        canonical_structure_identity(smiles)


def test_canonical_structure_identity_rejects_missing_inchi_key(monkeypatch):
    monkeypatch.setattr(
        "ThermoScreening.thermo.redox_screening.inchi.InchiToInchiKey",
        lambda value: "",
    )

    with pytest.raises(TSValueError, match="derive an InChIKey"):
        canonical_structure_identity("CCO")


def test_audit_redox_dataset_reports_duplicates_and_conflicts():
    audit = audit_redox_dataset(
        ["CCO", "OCC", "CCN", "NCC"],
        {"E1": [-0.5, -0.49, -0.8, -0.72], "E2": [None, np.nan, -1.1, -1.1]},
        tolerance=0.02,
    )

    assert audit.unique_structures == 2
    assert audit.duplicate_groups == ((0, 1), (2, 3))
    assert len(audit.conflicts) == 1
    assert audit.conflicts[0].indices == (2, 3)
    assert audit.conflicts[0].spreads[0][0] == "E1"
    assert audit.conflicts[0].spreads[0][1] == pytest.approx(0.08)


def test_audit_redox_dataset_validates_columns_and_values():
    with pytest.raises(TSValueError, match="expected 1"):
        audit_redox_dataset(["CCO"], {"E1": [0.1, 0.2]})
    with pytest.raises(TSValueError, match="not numeric"):
        audit_redox_dataset(["CCO", "OCC"], {"E1": [0.1, "bad"]})
    with pytest.raises(TSValueError, match="non-negative"):
        audit_redox_dataset(["CCO"], tolerance=-1.0)
    with pytest.raises(TSValueError, match="infinite"):
        audit_redox_dataset(["CCO"], {"E1": [np.inf]})


def test_delta_model_recovers_a_transparent_linear_correction():
    approximate, reference, groups = _matched_series()

    model = DeltaRedoxModel.fit(approximate, reference, groups, ridge=0.0)
    prediction = model.predict(approximate, groups)

    assert prediction.corrected_potential == pytest.approx(reference, abs=1.0e-10)
    assert not prediction.extrapolated.any()
    assert np.all(np.isfinite(prediction.uncertainty))
    assert model.correction_coefficients["approximate_potential"] == pytest.approx(0.12)
    assert model.correction_coefficients["group:amino"] == pytest.approx(0.03)
    assert model.correction_coefficients["group:nitro"] == pytest.approx(-0.05)


def test_delta_model_accepts_repeated_group_names_as_counts():
    approximate, reference, groups = _matched_series()
    names = [[name for name, count in row.items() for _ in range(int(count))] for row in groups]

    model = DeltaRedoxModel.fit(approximate, reference, names, ridge=0.0)

    assert model.predict(approximate, groups).corrected_potential == pytest.approx(reference)


def test_delta_model_fits_supported_pair_interactions():
    approximate = np.array([-1.2, -0.7, -1.0, -0.8, -1.3, -0.6, -1.1, -0.9])
    groups = [
        {},
        {"a": 1},
        {"b": 1},
        {"a": 1, "b": 1},
        {"a": 2},
        {"b": 2},
        {"a": 2, "b": 1},
        {"a": 1, "b": 2},
    ]
    correction = []
    for value, row in zip(approximate, groups):
        a_count = row.get("a", 0)
        b_count = row.get("b", 0)
        correction.append(
            0.02
            + 0.05 * value
            + 0.01 * a_count
            - 0.03 * b_count
            + 0.08 * a_count * b_count
            + 0.04 * a_count * (a_count - 1) / 2
            - 0.02 * b_count * (b_count - 1) / 2
        )
    reference = approximate + correction

    model = DeltaRedoxModel.fit(
        approximate,
        reference,
        groups,
        include_interactions=True,
        min_interaction_count=1,
        ridge=0.0,
    )

    assert model.predict(approximate, groups).corrected_potential == pytest.approx(
        reference, abs=1.0e-10
    )
    assert "pair:a|b" in model.correction_coefficients


def test_delta_model_flags_unknown_groups_pairs_and_potential_range():
    approximate = [-1.1, -0.9, -0.7, -0.8]
    reference = [-1.0, -0.8, -0.6, -0.7]
    groups = [{"a": 1}, {"b": 1}, {"c": 1}, {"a": 1, "b": 1}]
    model = DeltaRedoxModel.fit(
        approximate,
        reference,
        groups,
        include_interactions=True,
        min_interaction_count=1,
    )

    prediction = model.predict([-0.8, -1.3], [{"a": 1, "c": 1}, {"unknown": 1}])

    assert prediction.extrapolated.tolist() == [True, True]
    assert "pair:a|c" in prediction.unknown_features[0]
    assert set(prediction.unknown_features[1]) == {
        "approximate_potential_range",
        "group:unknown",
    }
    assert np.isnan(prediction.uncertainty).all()


def test_grouped_delta_validation_keeps_labels_in_single_folds():
    approximate, reference, groups = _matched_series()
    labels = ["one", "one", "two", "two", "three", "three", "four", "four"]

    validation = grouped_delta_validation(
        approximate, reference, groups, labels, ridge=1.0e-6
    )
    fitted = DeltaRedoxModel.fit(
        approximate,
        reference,
        groups,
        validation_groups=labels,
        ridge=1.0e-6,
    )

    assert validation.predictions.shape == approximate.shape
    assert validation.folds == tuple(labels)
    assert validation.rmse >= 0.0
    assert fitted.validation is not None
    assert fitted.validation.rmse == pytest.approx(validation.rmse)
    assert np.all(np.isfinite(fitted.predict(approximate, groups).uncertainty))


@pytest.mark.parametrize(
    ("approximate", "reference", "groups", "message"),
    [
        ([1.0], [1.0], [{}], "At least two"),
        ([1.0, np.nan], [1.0, 2.0], [{}, {}], "finite"),
        ([1.0, 2.0], [1.0], [{}, {}], "same length"),
        ([1.0, 2.0], [1.0, 2.0], [{}], "expected 2"),
        ([1.0, 2.0], [1.0, 2.0], [{"nitro": -1}, {}], "non-negative"),
        ([1.0, 2.0], [1.0, 2.0], [{"nitro": 0.5}, {}], "integer"),
        ([[1.0, 2.0]], [1.0], [{}], "one-dimensional"),
        ([1.0, 2.0], [1.0, 2.0], ["nitro", []], "mapping or an iterable"),
        ([1.0, 2.0], [1.0, 2.0], [{1: 1}, {}], "non-empty strings"),
        ([1.0, 2.0], [1.0, 2.0], [{"nitro": "many"}, {}], "not numeric"),
    ],
)
def test_delta_model_rejects_invalid_training_data(
    approximate, reference, groups, message
):
    with pytest.raises(TSValueError, match=message):
        DeltaRedoxModel.fit(approximate, reference, groups)


def test_delta_model_rejects_invalid_regularization_options():
    with pytest.raises(TSValueError, match="ridge"):
        DeltaRedoxModel.fit([1, 2], [1, 2], [{}, {}], ridge=-1)
    with pytest.raises(TSValueError, match="min_interaction_count"):
        DeltaRedoxModel.fit(
            [1, 2],
            [1, 2],
            [{}, {}],
            include_interactions=True,
            min_interaction_count=0,
        )


def test_grouped_delta_validation_rejects_invalid_folds():
    with pytest.raises(TSValueError, match="at least two distinct"):
        grouped_delta_validation([1, 2, 3], [1, 2, 3], [{}, {}, {}], ["same"] * 3)
    with pytest.raises(TSValueError, match="fewer than two training"):
        grouped_delta_validation([1, 2, 3], [1, 2, 3], [{}, {}, {}], ["a", "a", "b"])
    with pytest.raises(TSValueError, match="same length"):
        grouped_delta_validation([1, 2], [1], [{}, {}], ["a", "b"])
    with pytest.raises(TSValueError, match="expected 2"):
        grouped_delta_validation([1, 2], [1, 2], [{}, {}], ["a"])
    with pytest.raises(TSValueError, match="hashable"):
        grouped_delta_validation([1, 2], [1, 2], [{}, {}], [["a"], ["b"]])


def test_balanced_group_folds_never_splits_a_group():
    groups = ["large"] * 5 + ["medium"] * 3 + ["small-a"] * 2 + ["small-b"]

    folds = balanced_group_folds(groups, n_folds=3)

    assert len(folds) == len(groups)
    assert all(len({fold for fold, value in zip(folds, groups) if value == group}) == 1 for group in set(groups))
    loads = [folds.count(index) for index in range(3)]
    assert max(loads) - min(loads) <= 2


def test_balanced_group_folds_validates_inputs():
    with pytest.raises(TSValueError, match="at least 2"):
        balanced_group_folds(["a", "b"], n_folds=1)
    with pytest.raises(TSValueError, match="only 2 distinct"):
        balanced_group_folds(["a", "b"], n_folds=3)
    with pytest.raises(TSValueError, match="hashable"):
        balanced_group_folds([["a"], ["b"]], n_folds=2)


def test_pareto_front_preserves_tradeoffs_and_equal_points():
    records = [
        {"name": "a", "potential": -0.8, "mass": 200, "uncertainty": 0.05},
        {"name": "b", "potential": -0.7, "mass": 210, "uncertainty": 0.04},
        {"name": "c", "potential": -0.9, "mass": 220, "uncertainty": 0.10},
        {"name": "d", "potential": -0.8, "mass": 200, "uncertainty": 0.05},
    ]

    front = pareto_front(
        records,
        {"potential": "max", "mass": "min", "uncertainty": "min"},
    )

    assert [record["name"] for record in front] == ["a", "b", "d"]


def test_pareto_front_handles_empty_input_and_validates_objectives():
    assert pareto_front([], {"potential": "max"}) == []
    with pytest.raises(TSValueError, match="At least one"):
        pareto_front([{"potential": 1.0}], {})
    with pytest.raises(TSValueError, match="'min' or 'max'"):
        pareto_front([{"potential": 1.0}], {"potential": "up"})
    with pytest.raises(TSValueError, match="finite"):
        pareto_front([{"potential": np.nan}], {"potential": "max"})
    with pytest.raises(TSValueError, match="missing or non-numeric"):
        pareto_front([{"mass": 10.0}], {"potential": "max"})
