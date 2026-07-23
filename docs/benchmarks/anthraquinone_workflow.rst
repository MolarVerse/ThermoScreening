Anthraquinone workflow regression
=================================

The test suite contains a compact three-molecule anthraquinone fixture covering
the parent compound and two hydroxy substitution positions. It validates:

- deterministic SMILES embedding and molecular formulas;
- construction of the oxidized, singly reduced and doubly reduced states;
- separate reference calibration of both reduction steps;
- potential-inversion classification; and
- reproducible input-set provenance.

The fixture tests workflow invariants with controlled state energies. It is not
an absolute-potential accuracy benchmark and does not copy results from an
external dataset. Quantitative comparisons require the same structures,
electronic-structure method, solvent model, reference electrode and standard
states.
