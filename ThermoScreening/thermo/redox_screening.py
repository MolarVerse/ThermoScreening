"""Dataset-level tools for calibrated redox screening."""

from collections import Counter
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, replace
import math

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import inchi

from ThermoScreening.exceptions import TSValueError


@dataclass(frozen=True)
class MoleculeIdentity:
    """Canonical structure identity derived from a SMILES string."""

    canonical_smiles: str
    structure_id: str


@dataclass(frozen=True)
class DuplicateConflict:
    """Duplicate structure whose reported values disagree beyond tolerance."""

    structure_id: str
    indices: tuple[int, ...]
    spreads: tuple[tuple[str, float], ...]


@dataclass(frozen=True)
class DatasetAudit:
    """Canonical identities and duplicate findings for a redox dataset."""

    identities: tuple[MoleculeIdentity, ...]
    duplicate_groups: tuple[tuple[int, ...], ...]
    conflicts: tuple[DuplicateConflict, ...]

    @property
    def unique_structures(self):
        """Number of unique molecular structures in the dataset."""

        return len({identity.structure_id for identity in self.identities})


@dataclass(frozen=True)
class DeltaPrediction:
    """Corrected potentials and their model-domain diagnostics."""

    corrected_potential: np.ndarray
    correction: np.ndarray
    uncertainty: np.ndarray
    extrapolated: np.ndarray
    unknown_features: tuple[tuple[str, ...], ...]


@dataclass(frozen=True)
class CrossValidationResult:
    """Out-of-fold diagnostics from grouped delta-model validation."""

    predictions: np.ndarray
    errors: np.ndarray
    folds: tuple[object, ...]
    mae: float
    rmse: float
    bias: float


def canonical_structure_identity(smiles):
    """Return canonical isomeric SMILES and an InChIKey for ``smiles``."""

    if not isinstance(smiles, str) or not smiles.strip():
        raise TSValueError("SMILES must be a non-empty string.")
    molecule = Chem.MolFromSmiles(smiles)  # pylint: disable=no-member
    if molecule is None:
        raise TSValueError(f"Invalid SMILES: {smiles!r}.")

    canonical = Chem.MolToSmiles(  # pylint: disable=no-member
        molecule, canonical=True, isomericSmiles=True
    )
    with rdBase.BlockLogs():  # pylint: disable=c-extension-no-member
        fixed_h_inchi = inchi.MolToInchi(molecule, options="/FixedH")
    structure_id = inchi.InchiToInchiKey(fixed_h_inchi)
    if not structure_id:
        raise TSValueError(f"Could not derive an InChIKey from SMILES: {smiles!r}.")
    return MoleculeIdentity(canonical, structure_id)


def audit_redox_dataset(smiles, potentials=None, tolerance=0.02):
    """
    Canonicalize structures and identify duplicate measurement conflicts.

    ``potentials`` maps labels such as ``"E1"`` to sequences in volts. Missing
    values may be ``None`` or ``NaN``. Duplicate structures are always reported;
    a conflict is reported when any potential spans more than ``tolerance``.
    """

    smiles = tuple(smiles)
    if not math.isfinite(tolerance) or tolerance < 0.0:
        raise TSValueError("tolerance must be a finite, non-negative value.")

    values = {}
    for name, column in (potentials or {}).items():
        column = tuple(column)
        if len(column) != len(smiles):
            raise TSValueError(
                f"Potential column {name!r} has {len(column)} rows; expected {len(smiles)}."
            )
        parsed = []
        for index, value in enumerate(column):
            if value is None or value == "":
                parsed.append(np.nan)
                continue
            try:
                number = float(value)
            except (TypeError, ValueError) as exc:
                raise TSValueError(
                    f"Potential {name!r} at row {index} is not numeric: {value!r}."
                ) from exc
            if math.isinf(number):
                raise TSValueError(
                    f"Potential {name!r} at row {index} must not be infinite."
                )
            parsed.append(number)
        values[str(name)] = tuple(parsed)

    identities = tuple(canonical_structure_identity(value) for value in smiles)
    grouped = {}
    for index, identity in enumerate(identities):
        grouped.setdefault(identity.structure_id, []).append(index)

    duplicate_groups = tuple(
        tuple(indices) for indices in grouped.values() if len(indices) > 1
    )
    conflicts = []
    for indices in duplicate_groups:
        spreads = []
        for name, column in values.items():
            finite = [column[index] for index in indices if math.isfinite(column[index])]
            if len(finite) > 1:
                spread = max(finite) - min(finite)
                if spread > tolerance:
                    spreads.append((name, spread))
        if spreads:
            conflicts.append(
                DuplicateConflict(
                    identities[indices[0]].structure_id,
                    indices,
                    tuple(spreads),
                )
            )

    return DatasetAudit(identities, duplicate_groups, tuple(conflicts))


def _as_finite_vector(values, name):
    vector = np.asarray(values, dtype=float)
    if vector.ndim != 1:
        raise TSValueError(f"{name} must be one-dimensional.")
    if not np.all(np.isfinite(vector)):
        raise TSValueError(f"{name} must contain only finite values.")
    return vector


def _group_counts(groups):
    if isinstance(groups, Mapping):
        items = groups.items()
    else:
        if isinstance(groups, (str, bytes)) or not isinstance(groups, Iterable):
            raise TSValueError(
                "Each functional-group descriptor must be a mapping or an iterable of names."
            )
        items = Counter(groups).items()

    counts = {}
    for name, value in items:
        if not isinstance(name, str) or not name:
            raise TSValueError("Functional-group names must be non-empty strings.")
        try:
            count = float(value)
        except (TypeError, ValueError) as exc:
            raise TSValueError(f"Count for functional group {name!r} is not numeric.") from exc
        if not math.isfinite(count) or count < 0.0 or not count.is_integer():
            raise TSValueError(
                f"Count for functional group {name!r} must be a non-negative integer."
            )
        if count:
            counts[name] = int(count)
    return counts


def _normalise_descriptors(descriptors, expected_length):
    descriptors = tuple(_group_counts(groups) for groups in descriptors)
    if len(descriptors) != expected_length:
        raise TSValueError(
            f"functional_groups has {len(descriptors)} rows; expected {expected_length}."
        )
    return descriptors


def _pair_value(counts, first, second):
    if first == second:
        count = counts.get(first, 0.0)
        return count * (count - 1.0) / 2.0
    return counts.get(first, 0.0) * counts.get(second, 0.0)


def _available_interactions(descriptors, minimum_count):
    group_names = sorted({name for row in descriptors for name in row})
    interactions = []
    for first_index, first in enumerate(group_names):
        for second in group_names[first_index:]:
            support = sum(_pair_value(row, first, second) for row in descriptors)
            if support >= minimum_count:
                interactions.append((first, second))
    return tuple(interactions)


def _design_matrix(
    approximate,
    descriptors,
    approximate_center,
    group_names,
    interactions,
):
    matrix = np.ones(
        (len(approximate), 2 + len(group_names) + len(interactions)), dtype=float
    )
    matrix[:, 1] = approximate - approximate_center
    for column, name in enumerate(group_names, start=2):
        matrix[:, column] = [row.get(name, 0.0) for row in descriptors]
    offset = 2 + len(group_names)
    for column, (first, second) in enumerate(interactions, start=offset):
        matrix[:, column] = [
            _pair_value(row, first, second) for row in descriptors
        ]
    return matrix


@dataclass(frozen=True)
class DeltaRedoxModel:
    """Transparent linear correction from approximate to reference potentials."""

    approximate_center: float
    group_names: tuple[str, ...]
    interactions: tuple[tuple[str, str], ...]
    coefficients: np.ndarray
    covariance_factor: np.ndarray
    residual_std: float
    approximate_range: tuple[float, float]
    ridge: float
    include_interactions: bool
    validation: CrossValidationResult | None = None

    @classmethod
    def fit(
        cls,
        approximate,
        reference,
        functional_groups,
        *,
        include_interactions=False,
        min_interaction_count=5,
        ridge=1.0e-8,
        validation_groups=None,
    ):
        """
        Fit ``reference - approximate`` using group counts and optional pairs.

        ``validation_groups`` assigns each row to a held-out fold. Rows with the
        same label are never split between training and validation.
        """

        approximate = _as_finite_vector(approximate, "approximate")
        reference = _as_finite_vector(reference, "reference")
        if len(reference) != len(approximate):
            raise TSValueError("approximate and reference must have the same length.")
        descriptors = _normalise_descriptors(functional_groups, len(approximate))
        model = cls._fit(
            approximate,
            reference,
            descriptors,
            include_interactions=include_interactions,
            min_interaction_count=min_interaction_count,
            ridge=ridge,
        )
        if validation_groups is None:
            return model
        validation = grouped_delta_validation(
            approximate,
            reference,
            descriptors,
            validation_groups,
            include_interactions=include_interactions,
            min_interaction_count=min_interaction_count,
            ridge=ridge,
        )
        return replace(model, validation=validation)

    @classmethod
    def _fit(
        cls,
        approximate,
        reference,
        descriptors,
        *,
        include_interactions,
        min_interaction_count,
        ridge,
    ):
        if len(approximate) < 2:
            raise TSValueError("At least two matched rows are required to fit a delta model.")
        if not math.isfinite(ridge) or ridge < 0.0:
            raise TSValueError("ridge must be a finite, non-negative value.")
        if (
            isinstance(min_interaction_count, bool)
            or not isinstance(min_interaction_count, int)
            or min_interaction_count < 1
        ):
            raise TSValueError("min_interaction_count must be at least 1.")

        group_names = tuple(sorted({name for row in descriptors for name in row}))
        interactions = (
            _available_interactions(descriptors, min_interaction_count)
            if include_interactions
            else ()
        )
        center = float(np.mean(approximate))
        matrix = _design_matrix(
            approximate, descriptors, center, group_names, interactions
        )
        target = reference - approximate
        penalty = np.eye(matrix.shape[1]) * ridge
        penalty[0, 0] = 0.0
        cross_product = matrix.T @ matrix
        normal = cross_product + penalty
        inverse = np.linalg.pinv(normal)
        coefficients = inverse @ matrix.T @ target
        residual = target - matrix @ coefficients
        degrees_of_freedom = max(len(target) - np.linalg.matrix_rank(matrix), 1)
        residual_std = float(np.sqrt(np.sum(residual**2) / degrees_of_freedom))
        return cls(
            center,
            group_names,
            interactions,
            coefficients,
            inverse @ cross_product @ inverse,
            residual_std,
            (float(np.min(approximate)), float(np.max(approximate))),
            float(ridge),
            bool(include_interactions),
        )

    @property
    def correction_coefficients(self):
        """Model coefficients with an uncentered approximate-potential term."""

        result = {
            "intercept": float(
                self.coefficients[0]
                - self.coefficients[1] * self.approximate_center
            ),
            "approximate_potential": float(self.coefficients[1]),
        }
        offset = 2
        for name, value in zip(self.group_names, self.coefficients[offset:]):
            result[f"group:{name}"] = float(value)
        offset += len(self.group_names)
        for pair, value in zip(self.interactions, self.coefficients[offset:]):
            result[f"pair:{pair[0]}|{pair[1]}"] = float(value)
        return result

    def predict(self, approximate, functional_groups):
        """Correct approximate potentials and flag model-domain extrapolation."""

        approximate = _as_finite_vector(approximate, "approximate")
        descriptors = _normalise_descriptors(functional_groups, len(approximate))
        matrix = _design_matrix(
            approximate,
            descriptors,
            self.approximate_center,
            self.group_names,
            self.interactions,
        )
        correction = matrix @ self.coefficients
        leverage = np.einsum(
            "ij,jk,ik->i", matrix, self.covariance_factor, matrix
        )
        uncertainty_floor = self.residual_std
        if self.validation is not None:
            uncertainty_floor = max(uncertainty_floor, self.validation.rmse)
        uncertainty = uncertainty_floor * np.sqrt(np.maximum(1.0 + leverage, 1.0))

        known_groups = set(self.group_names)
        known_pairs = set(self.interactions)
        unknown = []
        for potential, row in zip(approximate, descriptors):
            missing = {f"group:{name}" for name in row if name not in known_groups}
            if self.include_interactions:
                row_pairs = _available_interactions((row,), 1)
                missing.update(
                    f"pair:{first}|{second}"
                    for first, second in row_pairs
                    if (first, second) not in known_pairs
                )
            if potential < self.approximate_range[0] or potential > self.approximate_range[1]:
                missing.add("approximate_potential_range")
            unknown.append(tuple(sorted(missing)))

        extrapolated = np.asarray([bool(items) for items in unknown], dtype=bool)
        uncertainty[extrapolated] = np.nan
        return DeltaPrediction(
            approximate + correction,
            correction,
            uncertainty,
            extrapolated,
            tuple(unknown),
        )


def grouped_delta_validation(
    approximate,
    reference,
    functional_groups,
    validation_groups,
    *,
    include_interactions=False,
    min_interaction_count=5,
    ridge=1.0e-8,
):
    """Validate a delta model while keeping each labelled group in one fold."""

    approximate = _as_finite_vector(approximate, "approximate")
    reference = _as_finite_vector(reference, "reference")
    if len(reference) != len(approximate):
        raise TSValueError("approximate and reference must have the same length.")
    descriptors = _normalise_descriptors(functional_groups, len(approximate))
    labels = tuple(validation_groups)
    if len(labels) != len(approximate):
        raise TSValueError(
            f"validation_groups has {len(labels)} rows; expected {len(approximate)}."
        )
    try:
        unique_labels = tuple(dict.fromkeys(labels))
    except TypeError as exc:
        raise TSValueError("validation_groups must contain hashable labels.") from exc
    if len(unique_labels) < 2:
        raise TSValueError("validation_groups must contain at least two distinct labels.")

    predictions = np.empty(len(approximate), dtype=float)
    for label in unique_labels:
        test = np.fromiter((candidate == label for candidate in labels), dtype=bool)
        train = ~test
        if np.count_nonzero(train) < 2:
            raise TSValueError(
                f"Validation fold {label!r} leaves fewer than two training rows."
            )
        model = DeltaRedoxModel._fit(
            approximate[train],
            reference[train],
            tuple(row for row, keep in zip(descriptors, train) if keep),
            include_interactions=include_interactions,
            min_interaction_count=min_interaction_count,
            ridge=ridge,
        )
        prediction = model.predict(
            approximate[test],
            tuple(row for row, keep in zip(descriptors, test) if keep),
        )
        predictions[test] = prediction.corrected_potential

    errors = predictions - reference
    return CrossValidationResult(
        predictions,
        errors,
        labels,
        float(np.mean(np.abs(errors))),
        float(np.sqrt(np.mean(errors**2))),
        float(np.mean(errors)),
    )


def balanced_group_folds(group_ids, n_folds=5):
    """Assign complete groups to deterministic folds with similar row counts."""

    group_ids = tuple(group_ids)
    if isinstance(n_folds, bool) or not isinstance(n_folds, int) or n_folds < 2:
        raise TSValueError("n_folds must be an integer of at least 2.")
    try:
        counts = Counter(group_ids)
    except TypeError as exc:
        raise TSValueError("group_ids must contain hashable values.") from exc
    if len(counts) < n_folds:
        raise TSValueError(
            f"Cannot build {n_folds} folds from only {len(counts)} distinct groups."
        )

    loads = [0] * n_folds
    assignment = {}
    ordered = sorted(counts.items(), key=lambda item: (-item[1], repr(item[0])))
    for group_id, count in ordered:
        fold = min(range(n_folds), key=lambda index: (loads[index], index))
        assignment[group_id] = fold
        loads[fold] += count
    return tuple(assignment[group_id] for group_id in group_ids)


def pareto_front(records, objectives):
    """
    Return the non-dominated records for named minimization/maximization goals.

    ``objectives`` maps record keys to ``"min"`` or ``"max"``. Equal objective
    vectors are retained, and the input order is preserved.
    """

    records = list(records)
    if not objectives:
        raise TSValueError("At least one Pareto objective is required.")
    if not records:
        return []

    columns = []
    for name, direction in objectives.items():
        if direction not in {"min", "max"}:
            raise TSValueError(
                f"Objective {name!r} direction must be 'min' or 'max', got {direction!r}."
            )
        try:
            column = np.asarray([record[name] for record in records], dtype=float)
        except (KeyError, TypeError, ValueError) as exc:
            raise TSValueError(f"Objective {name!r} is missing or non-numeric.") from exc
        if not np.all(np.isfinite(column)):
            raise TSValueError(f"Objective {name!r} must contain only finite values.")
        columns.append(column if direction == "min" else -column)

    points = np.column_stack(columns)
    unique_points, inverse = np.unique(points, axis=0, return_inverse=True)
    frontier = []
    for index, point in enumerate(unique_points):
        if not frontier:
            frontier.append(index)
            continue
        current = unique_points[frontier]
        dominated = np.any(
            np.all(current <= point, axis=1) & np.any(current < point, axis=1)
        )
        if dominated:
            continue
        removes = np.all(point <= current, axis=1) & np.any(point < current, axis=1)
        frontier = [
            kept for kept, remove in zip(frontier, removes) if not remove
        ]
        frontier.append(index)

    selected = set(frontier)
    return [record for record, point_index in zip(records, inverse) if point_index in selected]
