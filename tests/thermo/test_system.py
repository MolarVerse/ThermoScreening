import unittest
import warnings
import math
import numpy as np
from ase.atoms import Atoms
import ThermoScreening.thermo.system as system_module
from ThermoScreening.thermo.system import System, dim, dof, linearity, rotational_symmetry_number, default_spin, frequency_dof, check_frequency_length
from ThermoScreening.thermo.atoms import Atom
from ThermoScreening.thermo.cell import Cell
from ThermoScreening.exceptions import TSValueError
import pytest


def _atoms(spec):
    return [
        Atom(symbol=symbol, position=np.array(position, dtype=float))
        for symbol, position in spec
    ]


def _ring(symbol, radius, z_coord, count, phase=0.0):
    return [
        (
            symbol,
            (
                radius * math.cos(phase + 2 * math.pi * index / count),
                radius * math.sin(phase + 2 * math.pi * index / count),
                z_coord,
            ),
        )
        for index in range(count)
    ]


SYMMETRY_EDGE_CASES = [
    (
        "linear_h2",
        [("H", (0, 0, -0.37)), ("H", (0, 0, 0.37))],
        2,
        "D*h",
    ),
    (
        "linear_co",
        [("C", (0, 0, -0.64)), ("O", (0, 0, 0.49))],
        1,
        "C*v",
    ),
    (
        "linear_co2",
        [("O", (0, 0, -1.16)), ("C", (0, 0, 0)), ("O", (0, 0, 1.16))],
        2,
        "D*h",
    ),
    (
        "water",
        [
            ("O", (0.000000, 0.000000, 0.000000)),
            ("H", (0.957200, 0.000000, 0.000000)),
            ("H", (-0.239987, 0.927297, 0.000000)),
        ],
        2,
        "C2v",
    ),
    (
        "ammonia",
        [
            ("N", (0, 0, 0.116)),
            ("H", (0, 0.939, -0.272)),
            ("H", (0.813, -0.4695, -0.272)),
            ("H", (-0.813, -0.4695, -0.272)),
        ],
        3,
        "C3v",
    ),
    (
        "methane",
        [
            ("C", (0, 0, 0)),
            ("H", (0.629, 0.629, 0.629)),
            ("H", (-0.629, -0.629, 0.629)),
            ("H", (-0.629, 0.629, -0.629)),
            ("H", (0.629, -0.629, -0.629)),
        ],
        12,
        "Td",
    ),
    ("bf3", [("B", (0, 0, 0))] + _ring("F", 1.3, 0, 3), 6, "D3h"),
    (
        "square_planar_sf4_geometry",
        [
            ("S", (0, 0, 0)),
            ("F", (1, 0, 0)),
            ("F", (-1, 0, 0)),
            ("F", (0, 1, 0)),
            ("F", (0, -1, 0)),
        ],
        8,
        "D4h",
    ),
    (
        "allene",
        [
            ("C", (0, 0, 0)),
            ("C", (0, 0, -1.3)),
            ("C", (0, 0, 1.3)),
            ("H", (0.9, 0, -1.3)),
            ("H", (-0.9, 0, -1.3)),
            ("H", (0, 0.9, 1.3)),
            ("H", (0, -0.9, 1.3)),
        ],
        4,
        "D2d",
    ),
    (
        "ethene",
        [
            ("C", (-0.67, 0, 0)),
            ("C", (0.67, 0, 0)),
            ("H", (-1.232, 0.928, 0)),
            ("H", (-1.232, -0.928, 0)),
            ("H", (1.232, 0.928, 0)),
            ("H", (1.232, -0.928, 0)),
        ],
        4,
        "D2h",
    ),
    (
        "ethane_staggered",
        [("C", (0, 0, -0.77)), ("C", (0, 0, 0.77))]
        + _ring("H", 1.0, -1.27, 3)
        + _ring("H", 1.0, 1.27, 3, math.pi / 3),
        6,
        "D3d",
    ),
    (
        "ethane_eclipsed",
        [("C", (0, 0, -0.77)), ("C", (0, 0, 0.77))]
        + _ring("H", 1.0, -1.27, 3)
        + _ring("H", 1.0, 1.27, 3),
        6,
        "D3h",
    ),
    (
        "benzene",
        [
            ("C", (1.397, 0, 0)),
            ("C", (0.6985, 1.2098, 0)),
            ("C", (-0.6985, 1.2098, 0)),
            ("C", (-1.397, 0, 0)),
            ("C", (-0.6985, -1.2098, 0)),
            ("C", (0.6985, -1.2098, 0)),
            ("H", (2.479, 0, 0)),
            ("H", (1.2395, 2.1468, 0)),
            ("H", (-1.2395, 2.1468, 0)),
            ("H", (-2.479, 0, 0)),
            ("H", (-1.2395, -2.1468, 0)),
            ("H", (1.2395, -2.1468, 0)),
        ],
        12,
        "D6h",
    ),
    (
        "sf6",
        [
            ("S", (0, 0, 0)),
            ("F", (1.56, 0, 0)),
            ("F", (-1.56, 0, 0)),
            ("F", (0, 1.56, 0)),
            ("F", (0, -1.56, 0)),
            ("F", (0, 0, 1.56)),
            ("F", (0, 0, -1.56)),
        ],
        24,
        "Oh",
    ),
    (
        "hof",
        [("O", (0, 0, 0)), ("H", (0.96, 0.02, 0.01)), ("F", (-0.25, 1.39, 0.13))],
        1,
        "Cs",
    ),
    (
        "substituted_tetrahedron",
        [
            ("C", (0, 0, 0)),
            ("H", (1, 1, 1)),
            ("N", (-1, -1, 1)),
            ("O", (-1, 1, -1)),
            ("F", (1, -1, -1)),
        ],
        1,
        "C1",
    ),
]


def _transformed_spec(spec, mode):
    symbols = [symbol for symbol, _ in spec]
    coordinates = np.array([position for _, position in spec], dtype=float)

    if mode == "translated":
        coordinates = coordinates + np.array([3.2, -1.7, 0.9])
    elif mode == "rotated":
        axis = np.array([0.2, -0.7, 0.68])
        axis = axis / np.linalg.norm(axis)
        angle = 1.137
        cross_product_matrix = np.array(
            [
                [0, -axis[2], axis[1]],
                [axis[2], 0, -axis[0]],
                [-axis[1], axis[0], 0],
            ]
        )
        rotation = (
            np.eye(3)
            + math.sin(angle) * cross_product_matrix
            + (1 - math.cos(angle)) * (cross_product_matrix @ cross_product_matrix)
        )
        coordinates = coordinates @ rotation.T
    elif mode == "permuted":
        order = list(reversed(range(len(spec))))
        symbols = [symbols[index] for index in order]
        coordinates = coordinates[order]
    else:
        raise ValueError(mode)

    return list(zip(symbols, map(tuple, coordinates)))


def _symmetry_case_ids(value):
    if isinstance(value, str):
        return value
    return None



def test_dim():
    atoms = []
    with pytest.raises(TSValueError) as e:
        dim(atoms)
    assert str(e.value) == "The number of atoms must be greater than 0."
    
def test_dof():
    atoms = []
    with pytest.raises(TSValueError) as e:
        dof(atoms)
    assert str(e.value) == "The number of atoms must be greater than 0."


def test_frequency_dof_keeps_highest_modes():
    # nine raw modes (ascending), keep the top dof=3 vibrations
    freqs = np.arange(1.0, 10.0)
    kept = frequency_dof(freqs, 3)
    assert list(kept) == [7.0, 8.0, 9.0]


def test_frequency_dof_identity_when_length_matches():
    freqs = np.array([100.0, 200.0, 300.0])
    assert list(frequency_dof(freqs, 3)) == [100.0, 200.0, 300.0]


def test_frequency_dof_raises_on_underflow():
    # too few frequencies must raise instead of wrapping/duplicating modes
    with pytest.raises(TSValueError, match="expected at least 4 frequencies"):
        frequency_dof(np.array([1.0, 2.0, 3.0]), 4)


def test_check_frequency_length_validates_input_count():
    assert check_frequency_length(np.arange(9.0), 3) is True
    assert check_frequency_length(np.array([1.0, 2.0, 3.0]), 3) is True
    assert check_frequency_length(np.array([1.0, 2.0]), 3) is False

def test_linearity():
    atoms = [Atom(symbol='H', position=np.array([0, 0, 0]))]

    with pytest.raises(TSValueError) as e:
        linearity(atoms)
    assert str(e.value) == "Number of atoms must be greater than 1. The system is monoatomic."


def _atoms_at(spec, offset):
    return [
        Atom(symbol=symbol, position=np.array(position, dtype=float) + offset)
        for symbol, position in spec
    ]


@pytest.mark.parametrize("offset", [np.zeros(3), np.array([10.0, -5.0, 3.0])])
def test_linearity_is_translation_invariant(offset):
    # Inertia tensor must be built from center-of-mass-relocated coordinates,
    # so classification does not depend on where the molecule sits in space.
    co2 = [("C", [0, 0, 0]), ("O", [0, 0, 1.16]), ("O", [0, 0, -1.16])]
    # water is planar; when it lies in a coordinate plane the old second-moment
    # tensor produced a spurious zero eigenvalue and called it linear.
    water = [("O", [0, 0, 0]), ("H", [0.757, 0.586, 0]), ("H", [-0.757, 0.586, 0])]

    assert linearity(_atoms_at(co2, offset)) is True
    assert linearity(_atoms_at(water, offset)) is False


def test_linearity_detects_off_axis_linear_molecule():
    # A linear molecule not aligned to a Cartesian axis must still be linear.
    hcn = [("H", [0, 0, 0]), ("C", [0.6, 0.6, 0.6]), ("N", [1.2, 1.2, 1.2])]
    assert linearity(_atoms_at(hcn, np.zeros(3))) is True


def test_linearity_rejects_slightly_bent_near_linear_geometry():
    atoms = _atoms(
        [
            ("N", (0.0174573422, -1.1613421749, -0.0041534236)),
            ("C", (0.0025321236, -0.0034427793, 0.0017993915)),
            ("C", (-0.0161141932, 1.3722098653, 0.0093487294)),
            ("N", (-0.0326430226, 2.5300827589, 0.0160909027)),
        ]
    )

    assert linearity(atoms) is False
    assert dof(atoms) == 6


@pytest.mark.parametrize(
    "symbols,charge,expected",
    [
        (["O", "H", "H"], 0, 0.0),       # water, 10 electrons (even) -> singlet
        (["C", "H", "H", "H"], 0, 0.5),  # methyl radical, 9 electrons (odd) -> doublet
        (["O", "O"], 0, 0.0),            # O2: even electrons -> minimum-spin guesses singlet
        (["O", "H"], -1, 0.0),           # hydroxide anion, 10 electrons (even) -> singlet
    ],
)
def test_default_spin_from_electron_count(symbols, charge, expected):
    atoms = [Atom(symbol=s, position=np.array([i * 1.2, 0.0, 0.0]))
             for i, s in enumerate(symbols)]
    assert default_spin(atoms, charge) == expected


def test_system_accepts_explicit_spin():
    # an explicit spin overrides the guess, e.g. triplet O2 (even electrons)
    atoms = [Atom(symbol="O", position=np.array([0.0, 0, 0])),
             Atom(symbol="O", position=np.array([1.2, 0, 0]))]
    system = System(
        atoms, periodicity=False, cell=None, charge=0, spin=1.0,
        electronic_energy=0.0, vibrational_frequencies=np.array([1580.0]),
    )
    assert system.spin == 1.0


def test_system_rejects_negative_spin():
    atoms = [Atom(symbol="H", position=np.array([0.0, 0, 0])),
             Atom(symbol="H", position=np.array([1.0, 0, 0]))]
    with pytest.raises(TSValueError, match="spin must be non-negative"):
        System(
            atoms, periodicity=False, cell=None, charge=0, spin=-1.0,
            electronic_energy=0.0, vibrational_frequencies=np.array([1580.0]),
        )
             

def test_system_accepts_explicit_symmetry_number():
    atoms = [
        Atom(symbol="H", position=np.array([0.0, 0.0, 0.0])),
        Atom(symbol="H", position=np.array([0.0, 0.0, 0.74])),
    ]

    system = System(
        atoms,
        charge=0,
        electronic_energy=-1.0,
        vibrational_frequencies=np.array([4400.0]),
        symmetry_number=1,
    )

    assert system.rotational_symmetry_number == 1


@pytest.mark.parametrize("symmetry_number", [0, -1, True])
def test_system_rejects_invalid_symmetry_number(symmetry_number):
    atoms = [
        Atom(symbol="H", position=np.array([0.0, 0.0, 0.0])),
        Atom(symbol="H", position=np.array([0.0, 0.0, 0.74])),
    ]

    with pytest.raises(TSValueError, match="positive integer"):
        System(
            atoms,
            charge=0,
            electronic_energy=-1.0,
            vibrational_frequencies=np.array([4400.0]),
            symmetry_number=symmetry_number,
        )


def test_rotational_symmetry_number_accepts_property(monkeypatch):
    class FakePointGroupAnalyzer:
        get_rotational_symmetry_number = 7

        def __init__(self, molecule):
            self.molecule = molecule

    atoms = [Atom(symbol='H', position=np.array([0, 0, 0])), Atom(symbol='H', position=np.array([0, 0, 1]))]
    monkeypatch.setattr(system_module, "PointGroupAnalyzer", FakePointGroupAnalyzer)

    assert rotational_symmetry_number(atoms) == 7


def test_symmetry_analysis_converts_zero_imaginary_positions(monkeypatch):
    captured_coordinates = []

    class FakePointGroupAnalyzer:
        sch_symbol = "D*h"

        def __init__(self, molecule):
            self.molecule = molecule
            self.get_rotational_symmetry_number = 2
            captured_coordinates.append(np.array(molecule.cart_coords))

    atoms = [
        Atom(symbol="H", position=np.array([0.0 + 0.0j, 0.0, 0.0])),
        Atom(symbol="H", position=np.array([1.0 + 0.0j, 0.0, 0.0])),
    ]
    monkeypatch.setattr(system_module, "PointGroupAnalyzer", FakePointGroupAnalyzer)

    assert rotational_symmetry_number(atoms) == 2
    assert system_module.rotational_group_calc(atoms) == "D*h"
    assert not np.iscomplexobj(captured_coordinates[0])
    assert not np.iscomplexobj(captured_coordinates[1])


def test_real_position_converts_near_zero_imaginary_noise():
    position = system_module._real_position(np.array([1.0 + 1e-12j, 0.0, 0.0]))

    np.testing.assert_array_equal(position, np.array([1.0, 0.0, 0.0]))
    assert not np.iscomplexobj(position)


def test_point_group_analyzer_suppresses_pymatgen_complex_warning(monkeypatch):
    class FakePointGroupAnalyzer:
        sch_symbol = "C1"

        def __init__(self, molecule):
            self.molecule = molecule
            warnings.warn_explicit(
                "Casting complex values to real discards the imaginary part",
                system_module.ComplexWarning,
                filename="operations.py",
                lineno=1,
                module="pymatgen.core.operations",
            )

    monkeypatch.setattr(system_module, "PointGroupAnalyzer", FakePointGroupAnalyzer)

    with warnings.catch_warnings():
        warnings.simplefilter("error", system_module.ComplexWarning)
        analyzer = system_module._point_group_analyzer(object())

    assert analyzer.sch_symbol == "C1"


def test_point_group_analyzer_keeps_unexpected_complex_warning(monkeypatch):
    class FakePointGroupAnalyzer:
        def __init__(self, molecule):
            self.molecule = molecule
            warnings.warn_explicit(
                "unexpected complex warning",
                system_module.ComplexWarning,
                filename="other.py",
                lineno=1,
                module="other.module",
            )

    monkeypatch.setattr(system_module, "PointGroupAnalyzer", FakePointGroupAnalyzer)

    with warnings.catch_warnings():
        warnings.simplefilter("error", system_module.ComplexWarning)
        with pytest.raises(system_module.ComplexWarning, match="unexpected"):
            system_module._point_group_analyzer(object())


def test_symmetry_analysis_rejects_imaginary_positions():
    atoms = [
        Atom(symbol="H", position=np.array([0.0 + 1.0j, 0.0, 0.0])),
        Atom(symbol="H", position=np.array([1.0, 0.0, 0.0])),
    ]

    with pytest.raises(TSValueError, match="Atom positions must be real-valued"):
        rotational_symmetry_number(atoms)


def test_rotational_group_rejects_imaginary_positions():
    atoms = [
        Atom(symbol="H", position=np.array([0.0, 0.0, 0.0])),
        Atom(symbol="H", position=np.array([1.0, 1.0e-6j, 0.0])),
    ]

    with pytest.raises(TSValueError, match="Atom positions must be real-valued"):
        system_module.rotational_group_calc(atoms)


def test_water_symmetry_regression():
    atoms = [
        Atom(symbol="O", position=np.array([0.000000, 0.000000, 0.000000])),
        Atom(symbol="H", position=np.array([0.957200, 0.000000, 0.000000])),
        Atom(symbol="H", position=np.array([-0.239987, 0.927297, 0.000000])),
    ]

    assert rotational_symmetry_number(atoms) == 2
    assert system_module.rotational_group_calc(atoms) == "C2v"


def test_monoatomic_symmetry_regression():
    atoms = [Atom(symbol="H", position=np.array([0.0, 0.0, 0.0]))]

    assert rotational_symmetry_number(atoms) == 1
    assert system_module.rotational_group_calc(atoms) == "Kh"


def test_monoatomic_system_initializes():
    atoms = [Atom(symbol="H", position=np.array([0.0, 0.0, 0.0]))]

    system = System(
        atoms,
        periodicity=False,
        cell=None,
        solvation=None,
        solvent=None,
        charge=0,
        electronic_energy=-0.5,
        vibrational_frequencies=np.array([1.0, 2.0, 3.0]),
    )

    assert system.rotational_symmetry_number == 1
    assert system.rotational_group == "Kh"


def test_symmetry_analysis_rejects_overlapping_positions():
    atoms = _atoms([("H", (0.0, 0.0, 0.0)), ("H", (0.0, 0.0, 0.0))])

    with pytest.raises(TSValueError, match="must not overlap"):
        rotational_symmetry_number(atoms)
    with pytest.raises(TSValueError, match="must not overlap"):
        system_module.rotational_group_calc(atoms)


@pytest.mark.parametrize(
    "name,spec,expected_number,expected_group",
    SYMMETRY_EDGE_CASES,
    ids=_symmetry_case_ids,
)
def test_rotational_symmetry_edge_case_regressions(
    name,
    spec,
    expected_number,
    expected_group,
):
    del name
    atoms = _atoms(spec)

    assert rotational_symmetry_number(atoms) == expected_number
    assert system_module.rotational_group_calc(atoms) == expected_group


@pytest.mark.parametrize("mode", ["translated", "rotated", "permuted"])
@pytest.mark.parametrize(
    "name,spec,expected_number,expected_group",
    [
        case
        for case in SYMMETRY_EDGE_CASES
        if case[0] in {"linear_co2", "methane", "benzene", "sf6", "hof"}
    ],
    ids=_symmetry_case_ids,
)
def test_rotational_symmetry_is_coordinate_frame_invariant(
    name,
    spec,
    expected_number,
    expected_group,
    mode,
):
    del name
    atoms = _atoms(_transformed_spec(spec, mode))

    assert rotational_symmetry_number(atoms) == expected_number
    assert system_module.rotational_group_calc(atoms) == expected_group


@pytest.mark.parametrize(
    "name,spec,expected_number,expected_group",
    [
        case
        for case in SYMMETRY_EDGE_CASES
        if case[0] in {"water", "methane", "benzene", "sf6"}
    ],
    ids=_symmetry_case_ids,
)
def test_rotational_symmetry_tolerates_tiny_coordinate_noise(
    name,
    spec,
    expected_number,
    expected_group,
):
    del name
    rng = np.random.default_rng(314159)
    noisy_spec = [
        (symbol, tuple(np.asarray(position, dtype=float) + rng.normal(0, 1e-8, 3)))
        for symbol, position in spec
    ]
    atoms = _atoms(noisy_spec)

    assert rotational_symmetry_number(atoms) == expected_number
    assert system_module.rotational_group_calc(atoms) == expected_group


class TestSystem(unittest.TestCase):
    atoms = [
                Atom(symbol='H', position=np.array([0, 0, 0])),
                Atom(symbol='H', position=np.array([1, 0, 0])),
                Atom(symbol='H', position=np.array([2, 0, 0]))
            ]
    
    def test_invalid_init(self):
        with pytest.raises(TSValueError) as e:
            System(atoms=None,periodicity=False,cell=None,solvation=None,solvent=None,charge=0,electronic_energy=-33.6052447996,vibrational_frequencies=np.array([-1,-2,0.1,1,2,3,4,5,6]))
        assert str(e.value) == "Atoms must be provided."
        
        with pytest.raises(TSValueError) as e:
            System(atoms=self.atoms,periodicity=False,cell=None,solvation=None,solvent=None,charge=0,electronic_energy=-33.6052447996,vibrational_frequencies=None)
        assert str(e.value) == "Vibrational frequencies must be provided."
        
        with pytest.raises(TSValueError) as e:
            System(atoms=[],periodicity=False,cell=None,solvation=None,solvent=None,charge=0,electronic_energy=-33.6052447996,vibrational_frequencies=np.array([-1,-2,0.1,1,2,3,4,5,6]))
        assert str(e.value) == "The number of atoms must be greater than 0."
        
        
    def test_system(self):
        
            
        system = System(self.atoms,periodicity=False,cell=None,solvation=None,solvent=None,charge=0,electronic_energy=-33.6052447996,vibrational_frequencies=np.array([-1,-2,0.1,1,2,3,4,5,6]))

        np.testing.assert_array_equal(system.atom_names(), ['H', 'H', 'H'])

        np.testing.assert_array_equal(
            system.coord(), np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0]]))

        np.testing.assert_array_equal(system.atomic_masses()[0], [
                                      1.00794, 1.00794, 1.00794])

        assert len(system.atoms) == 3

        assert system.charge == 0

        assert system.spin == 0.5  # 3 H atoms -> 3 electrons (doublet, minimum-spin default)

        assert system.number_of_atoms == 3

        assert system.mass == 3.02382

        assert system.dim == 1

        assert system.dof == 3 * 3 - 5

        assert system.rotational_symmetry_number == 2

        assert system.rotational_group == 'D*h'

        assert system._spacegroup_number is None

        assert system._spacegroup is None

        assert system.periodicity is False

    def test_system_accepts_periodic_ndarray_cell(self):
        # read_xyz/read_gen hand back a raw ndarray cell for periodic input;
        # System must construct instead of raising on the spacegroup hint, and
        # degrade the (unused) spacegroup to None.
        cell = np.array([[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]])

        system = System(
            self.atoms,
            periodicity=True,
            cell=cell,
            solvation=None,
            solvent=None,
            charge=0,
            electronic_energy=-33.6052447996,
            vibrational_frequencies=np.array([-1, -2, 0.1, 1, 2, 3, 4, 5, 6]),
        )

        assert system.spacegroup_number is None
        assert system.spacegroup is None

    def test_spacegroup_properties_map_to_their_attributes(self):
        system = System(
            self.atoms,
            periodicity=False,
            cell=None,
            solvation=None,
            solvent=None,
            charge=0,
            electronic_energy=-33.6052447996,
            vibrational_frequencies=np.array([-1, -2, 0.1, 1, 2, 3, 4, 5, 6]),
        )
        system._spacegroup_number = 1
        system._spacegroup = "P1"

        assert system.spacegroup_number == 1
        assert system.spacegroup == "P1"

        np.testing.assert_array_equal(system.center_of_mass, [1, 0, 0])

        assert system.pbc == [False, False, False]

        assert system.cell == None

        assert system.solvation == False

        assert system.solvent == ''

        assert system.electronic_energy == -33.6052447996

        np.testing.assert_array_equal(
            system.vibrational_frequencies, np.array([-1, -2, 0.1, 1, 2, 3, 4, 5, 6]))

        np.testing.assert_array_equal(
            system.imaginary_frequencies, np.array([-1, -2]))

        np.testing.assert_array_equal(
            system._real_vibrational_frequencies, np.array([3, 4, 5, 6]))

        assert system._has_imaginary_frequencies == True

        assert system._check_frequency_length == True

        np.testing.assert_array_equal(system.x(), np.array([0, 1, 2]))

        np.testing.assert_array_equal(system.y(), np.array([0, 0, 0]))

        np.testing.assert_array_equal(system.z(), np.array([0, 0, 0]))

    def test_system_aq(self):

        atoms = [
            Atom(symbol='O', position=np.array(
                [0.00008112, 0.00000188, -2.05402286])),
            Atom(symbol='O', position=np.array(
                [-0.00009254, 0.00000502, -7.65886167])),
            Atom(symbol='C', position=np.array(
                [-3.68652488, -0.00000123, -5.56937006])),
            Atom(symbol='C', position=np.array(
                [-2.49045658, 0.00000104, -6.25271675])),
            Atom(symbol='C', position=np.array(
                [-1.22985769, 0.00000200, -5.58346285])),
            Atom(symbol='C', position=np.array(
                [-1.22981936, 0.00000095, -4.12932133])),
            Atom(symbol='C', position=np.array(
                [-2.49034390, -0.00000124, -3.45992559])),
            Atom(symbol='C', position=np.array(
                [-3.68647387, -0.00000247, -4.14315959])),
            Atom(symbol='C', position=np.array(
                [-0.00004268, 0.00000431, -6.32159897])),
            Atom(symbol='C', position=np.array(
                [0.00004092, 0.00000241, -3.39127855])),
            Atom(symbol='C', position=np.array(
                [1.22985982, 0.00000519, -4.12939595])),
            Atom(symbol='C', position=np.array(
                [1.22981336, 0.00000624, -5.58353992])),
            Atom(symbol='C', position=np.array(
                [2.49035536, 0.00001035, -6.25289256])),
            Atom(symbol='H', position=np.array(
                [2.48994553, 0.00001106, -7.34077798])),
            Atom(symbol='C', position=np.array(
                [3.68646897, 0.00001431, -5.56962708])),
            Atom(symbol='C', position=np.array(
                [3.68653147, 0.00001307, -4.14341821])),
            Atom(symbol='C', position=np.array(
                [2.49044502, 0.00000805, -3.46009695])),
            Atom(symbol='H', position=np.array(
                [-4.63230122, -0.00000198, -6.10807371])),
            Atom(symbol='H', position=np.array(
                [-2.49018198, 0.00000208, -7.34059729])),
            Atom(symbol='H', position=np.array(
                [-2.48990606, -0.00000205, -2.37203685])),
            Atom(symbol='H', position=np.array(
                [-4.63218409, -0.00000429, -3.60434535])),
            Atom(symbol='H', position=np.array(
                [4.63217171, 0.00001845, -6.10844479])),
            Atom(symbol='H', position=np.array(
                [4.63230600, 0.00001612, -3.60469800])),
            Atom(symbol='H', position=np.array(
                [2.49016578, 0.00000689, -2.37221336]))
        ]

        charge = -2

        vibrational_frequencies = np.array([-13.800,-2.54,6.36,25.33,46.98,53.790,90.3, 143.96, 207.59, 222.69, 229.17, 233.44, 281.5, 364.2, 388.97, 393.46, 418.36, 439.41, 449.72, 469.27, 513.78, 597.7, 606.64, 624.86, 629.79, 678.06, 686.75, 701.51, 717.21, 729.74, 734.03, 802.32, 807.81, 825.31, 896.02, 896.94, 898.48, 898.59, 959.91,
                                           960.79, 1043.79, 1047.42, 1056.91, 1120, 1127.16, 1130.62, 1131.3, 1192.85, 1202.44, 1275.69, 1321.66, 1323.16, 1333.16, 1397.04, 1433.18, 1504.85, 1514.35, 1518.2, 1543.43, 1596.63, 1625.67, 1677.6, 1699.25, 1701.29, 2979.33, 2980.88, 2986.11, 2986.64, 2995.14, 2995.94, 3000.44, 3001.33])

        elecronic_energy = -33.8387075830

        system = System(atoms, periodicity=False, cell=None, solvation=True, solvent="DMF", charge=charge,
                        electronic_energy=elecronic_energy, vibrational_frequencies=vibrational_frequencies)

        np.testing.assert_array_equal(system.atom_names(), [
                                      'O', 'O', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'H', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H'])

        np.testing.assert_array_equal(
            system.coord(), np.array([[0.00008112, 0.00000188, -2.05402286], [-0.00009254, 0.00000502, -7.65886167], [-3.68652488, -0.00000123, -5.56937006], [-2.49045658, 0.00000104, -6.25271675], [-1.22985769, 0.00000200, -5.58346285], [-1.22981936, 0.00000095, -4.12932133], [-2.49034390, -0.00000124, -3.45992559], [-3.68647387, -0.00000247, -4.14315959], [-0.00004268, 0.00000431, -6.32159897], [0.00004092, 0.00000241, -3.39127855], [1.22985982, 0.00000519, -4.12939595], [1.22981336, 0.00000624, -5.58353992], [2.49035536, 0.00001035, -6.25289256], [2.48994553, 0.00001106, -7.34077798], [3.68646897, 0.00001431, -5.56962708], [3.68653147, 0.00001307, -4.14341821], [2.49044502, 0.00000805, -3.46009695], [-4.63230122, -0.00000198, -6.10807371], [-2.49018198, 0.00000208, -7.34059729], [-2.48990606, -0.00000205, -2.37203685], [-4.63218409, -0.00000429, -3.60434535], [4.63217171, 0.00001845, -6.10844479], [4.63230600, 0.00001612, -3.60469800], [2.49016578, 0.00000689, -2.37221336]]))

        x = np.array([8.112000e-05 ,-9.254000e-05 ,-3.686525e+00 ,-2.490457e+00 ,-1.229858e+00 ,-1.229819e+00 ,-2.490344e+00 ,-3.686474e+00 ,-4.268000e-05 ,4.092000e-05 ,1.229860e+00 ,1.229813e+00 ,2.490355e+00 ,2.489946e+00 ,3.686469e+00 ,3.686531e+00 ,2.490445e+00 ,-4.632301e+00 ,-2.490182e+00 ,-2.489906e+00 ,-4.632184e+00 ,4.632172e+00 ,4.632306e+00 ,2.490166e+00])
        y= np.array([1.880000e-06 ,5.020000e-06 ,-1.230000e-06 ,1.040000e-06 ,2.000000e-06 ,9.500000e-07 ,-1.240000e-06 ,-2.470000e-06 ,4.310000e-06 ,2.410000e-06 ,5.190000e-06 ,6.240000e-06 ,1.035000e-05 ,1.106000e-05 ,1.431000e-05 ,1.307000e-05 ,8.050000e-06 ,-1.980000e-06 ,2.080000e-06 ,-2.050000e-06 ,-4.290000e-06 ,1.845000e-05 ,1.612000e-05 ,6.890000e-06])
        z = np.array([-2.054023e+00 ,-7.658862e+00 ,-5.569370e+00 ,-6.252717e+00 ,-5.583463e+00 ,-4.129321e+00 ,-3.459926e+00 ,-4.143160e+00 ,-6.321599e+00 ,-3.391279e+00 ,-4.129396e+00 ,-5.583540e+00 ,-6.252893e+00 ,-7.340778e+00 ,-5.569627e+00 ,-4.143418e+00 ,-3.460097e+00 ,-6.108074e+00 ,-7.340597e+00 ,-2.372037e+00 ,-3.604345e+00 ,-6.108445e+00 ,-3.604698e+00 ,-2.372213e+00])
        np.testing.assert_allclose(system.coord(), np.array([x, y, z]).T, atol=1e-6)
        np.testing.assert_allclose(system.x(), x, atol=1e-6)
        np.testing.assert_allclose(system.y(), y, atol=1e-6)
        np.testing.assert_allclose(system.z(), z, atol=1e-6)
        np.testing.assert_array_equal(system.atomic_masses()[0], [15.9994])
        np.testing.assert_array_equal(system.atomic_masses()[1], [15.9994])
        np.testing.assert_array_equal(system.atomic_masses()[2], [12.0107])
        np.testing.assert_array_equal(system.atomic_masses()[23], [1.00794])

        assert len(system.atoms) == 24

        assert system.charge == -2

        assert system.spin == 0

        assert system.number_of_atoms == 24

        assert system.mass == 208.21212

        assert system.dim == 3

        assert system.dof == 3 * 24 - 6

        assert system.rotational_symmetry_number == 4

        assert system.rotational_group == 'D2h'

        assert system._spacegroup_number == None

        assert system._spacegroup == None

        assert system.periodicity == False

        # use ase to get the center of mass
        ase_atoms = Atoms(symbols=system.atom_names(), positions=system.coord())
        print(ase_atoms.get_center_of_mass())
        def center_of_mass(system):
                total_mass = np.sum(system.atomic_masses())
                center_of_mass = np.zeros(3)
                for i in range(system.number_of_atoms):
                        center_of_mass += system.coord()[i] * system.atomic_masses()[i]
                center_of_mass /= total_mass
                return center_of_mass
        np.testing.assert_allclose(ase_atoms.get_center_of_mass(), center_of_mass(system), atol=1e-10)
        np.testing.assert_array_equal(system.center_of_mass, center_of_mass(system))



        assert system.pbc == [False, False, False]

        assert system.cell == None

        assert system.solvation == True

        assert system.solvent == 'DMF'

        assert system.electronic_energy == -33.8387075830

        np.testing.assert_array_equal(
            system.vibrational_frequencies, np.array([-13.800,-2.54,6.36,25.33,46.98,53.790,90.3, 143.96, 207.59, 222.69, 229.17, 233.44, 281.5, 364.2, 388.97, 393.46, 418.36, 439.41, 449.72, 469.27, 513.78, 597.7, 606.64, 624.86, 629.79, 678.06, 686.75, 701.51, 717.21, 729.74, 734.03, 802.32, 807.81, 825.31, 896.02, 896.94, 898.48, 898.59, 959.91, 960.79, 1043.79, 1047.42, 1056.91, 1120, 1127.16, 1130.62, 1131.3, 1192.85, 1202.44, 1275.69, 1321.66, 1323.16, 1333.16, 1397.04, 1433.18, 1504.85, 1514.35, 1518.2, 1543.43, 1596.63, 1625.67, 1677.6, 1699.25, 1701.29, 2979.33, 2980.88, 2986.11, 2986.64, 2995.14, 2995.94, 3000.44, 3001.33]))

        np.testing.assert_array_equal(
            system.imaginary_frequencies, np.array([-13.800,-2.54]))

        np.testing.assert_array_equal(
            system._real_vibrational_frequencies, np.array([90.3, 143.96, 207.59, 222.69, 229.17, 233.44, 281.5, 364.2, 388.97, 393.46, 418.36, 439.41, 449.72, 469.27, 513.78, 597.7, 606.64, 624.86, 629.79, 678.06, 686.75, 701.51, 717.21, 729.74, 734.03, 802.32, 807.81, 825.31, 896.02, 896.94, 898.48, 898.59, 959.91, 960.79, 1043.79, 1047.42, 1056.91, 1120, 1127.16, 1130.62, 1131.3, 1192.85, 1202.44, 1275.69, 1321.66, 1323.16, 1333.16, 1397.04, 1433.18, 1504.85, 1514.35, 1518.2, 1543.43, 1596.63, 1625.67, 1677.6, 1699.25, 1701.29, 2979.33, 2980.88, 2986.11, 2986.64, 2995.14, 2995.94, 3000.44, 3001.33]))

        assert system._has_imaginary_frequencies == True

        assert system._check_frequency_length == True



if __name__ == '__main__':
    unittest.main()
