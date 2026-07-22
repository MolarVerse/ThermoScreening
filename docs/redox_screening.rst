Calibrated redox screening
==========================

A single reference compound removes a constant offset from computed reduction
potentials. Across a chemically diverse series, the remaining error can still
depend on the approximate potential, functional groups, substitution positions,
and interactions between substituents. ``DeltaRedoxModel`` fits that residual:

.. math::

   \Delta E = E_{\mathrm{reference}} - E_{\mathrm{approximate}}

The model is a regularized linear regression. Its coefficients remain available
through ``correction_coefficients``. It can include integer functional-group
counts and selected pair terms without requiring a machine-learning dependency.

Data audit
----------

Canonicalize structures before splitting or fitting. The identity function uses
canonical isomeric SMILES and a fixed-hydrogen InChIKey, so stereoisomers,
protonation states, and fixed-hydrogen tautomers remain distinct.

.. code-block:: python

   from ThermoScreening.thermo import audit_redox_dataset

   audit = audit_redox_dataset(
       smiles,
       {"E1_DFT": dft_e1, "E1_DFTB": dftb_e1},
       tolerance=0.02,  # volts
   )

   print(audit.unique_structures)
   print(audit.duplicate_groups)
   print(audit.conflicts)

Every duplicate is reported. ``conflicts`` identifies structures whose values
differ by more than the tolerance. Inspect or exclude conflicting calculations;
do not silently average them. Keep the original calculation provenance,
convergence state, geometry, charge, spin, solvent, and method alongside each
record.

Grouped validation
------------------

Random row splits can place symmetry-equivalent structures or closely related
positional isomers in both training and validation. Define a group identifier
that represents the generalization claim, then assign complete groups to folds.
For a substituted scaffold, a sorted substituent-composition tuple is a useful
default because all positional isomers remain together.

.. code-block:: python

   from ThermoScreening.thermo import balanced_group_folds, DeltaRedoxModel

   # One mapping per molecule. Include position labels when their effects matter.
   descriptors = [
       {"nitro": 1, "nitro@beta": 1},
       {"hydroxy": 2, "hydroxy@alpha": 1, "hydroxy@beta": 1},
       {"amino": 1, "amino@alpha": 1},
   ]
   composition_groups = [
       ("nitro",),
       ("hydroxy", "hydroxy"),
       ("amino",),
   ]
   folds = balanced_group_folds(composition_groups, n_folds=3)

   model = DeltaRedoxModel.fit(
       dftb_e1,
       dft_e1,
       descriptors,
       ridge=1.0e-4,
       validation_groups=folds,
   )

   print(model.validation.mae)
   print(model.validation.rmse)
   print(model.correction_coefficients)

Fit separate models for the first and second reductions. They describe different
charge-state transitions and generally have different residual errors. Select
the ridge strength and descriptor scheme without consulting the final
experimental test set.

Prediction and domain checks
----------------------------

.. code-block:: python

   prediction = model.predict(candidate_dftb_e1, candidate_descriptors)

   corrected = prediction.corrected_potential
   uncertainty = prediction.uncertainty
   extrapolated = prediction.extrapolated
   reasons = prediction.unknown_features

The uncertainty combines the regression residual with leverage and, when
grouped validation was supplied, uses its RMSE as a lower bound. It is a model
diagnostic rather than a calibrated probabilistic confidence interval. An
unknown functional group, unsupported pair interaction, or potential outside
the training range sets ``extrapolated`` and leaves the uncertainty as ``NaN``.
Such candidates require higher-level calculations instead of automatic ranking.

Pair interactions are optional:

.. code-block:: python

   interaction_model = DeltaRedoxModel.fit(
       dftb_e1,
       dft_e1,
       descriptors,
       include_interactions=True,
       min_interaction_count=10,
       ridge=1.0e-3,
       validation_groups=folds,
   )

Use pair terms only when they improve grouped validation. A minimum support
avoids assigning unconstrained coefficients to rare combinations.

Candidate selection
-------------------

Do not replace electrochemical objectives with potential divided by molecular
mass. That ratio is not specific energy. Keep competing properties separate and
select the non-dominated candidates:

.. code-block:: python

   from ThermoScreening.thermo import pareto_front

   candidates = [
       {"name": "A", "E1": -0.65, "mass": 240.2, "uncertainty": 0.04},
       {"name": "B", "E1": -0.58, "mass": 278.3, "uncertainty": 0.03},
   ]
   selected = pareto_front(
       candidates,
       {"E1": "max", "mass": "min", "uncertainty": "min"},
   )

Add solubility, stability, reversibility, or a physically defined cell voltage
when those data are available. Pareto selection exposes the trade-offs; it does
not decide which trade-off is appropriate for an application.

Final validation
----------------

The corrected DFTB+ result remains a screening estimate. Recalculate Pareto-front
candidates, extrapolations, and any case where ``E2 >= E1`` with the higher-level
method. Check multiple conformers, confirm that optimized structures are minima,
and inspect spin state and charge localization for the radical anion and dianion.
Keep solvent, temperature, reference electrode, and standard-state conventions
consistent across computed and experimental values.

Reserve experimental measurements as a final test set. Do not use the same
measurements to select descriptors, tune regularization, calibrate the reference,
and report final accuracy.

Anthraquinone benchmark
-----------------------

The workflow was checked against the matched DFTB+/DFT table from the
`Anthraquinone-screening repository
<https://github.com/jolint1/Anthraquinone-screening>`_. Every row belonging to
one of the 22 duplicated canonical structures was excluded, leaving 5,035
singleton structures before filtering missing potentials. Five balanced folds
kept all positional isomers of a substituent composition together. The models
used a ridge value of ``1e-4`` and no pair interactions.

.. list-table:: Grouped five-fold errors against DFT
   :header-rows: 1
   :widths: 22 13 13 13 13

   * - Model
     - E1 MAE
     - E1 RMSE
     - E2 MAE
     - E2 RMSE
   * - Uncorrected DFTB+
     - 82.9 mV
     - 117.7 mV
     - 92.6 mV
     - 135.2 mV
   * - Global intercept and slope
     - 58.7 mV
     - 78.2 mV
     - 64.4 mV
     - 92.6 mV
   * - Substituent counts
     - 31.7 mV
     - 46.5 mV
     - 44.7 mV
     - 71.3 mV
   * - Substituent and alpha/beta position
     - 29.2 mV
     - 42.7 mV
     - 35.4 mV
     - 59.6 mV

The E1 and E2 comparisons contain 5,022 and 5,013 structures, respectively.
These numbers measure agreement with the available DFT calculations, not with
experiment. Experimental potentials measured under consistent conditions must
remain an independent final validation set.
