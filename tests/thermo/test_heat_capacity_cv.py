"""The reported heat capacity is the ideal-gas constant-volume Cv.

Cross-checked non-circularly against the closed-form ideal-gas Cv (re-derived
here from the partition functions) and pinned to reference values. ASE's
IdealGasThermo exposes no Cv accessor, so this is validated analytically rather
than against ASE.
"""

import numpy as np
import pytest
from ase.build import molecule

from ThermoScreening.thermo.api import run_thermo

_h = 6.62607015e-34
_c = 299792458.0
_kB = 1.380649e-23
_R = 8.314462618 / 4.184  # cal/(mol*K)
_T = 298.15


def _analytic_cv(wavenumbers, n_rot):
    """Closed-form ideal-gas Cv (cal/mol/K) from translation, rotation, vibration."""
    x = _h * _c * np.asarray(wavenumbers, float) * 100.0 / (_kB * _T)
    vib = _R * np.sum(x**2 * np.exp(x) / np.expm1(x) ** 2)
    return 1.5 * _R + 0.5 * n_rot * _R + vib


# molecule, vibrational wavenumbers (cm^-1), rotational dof, reference Cv
_CASES = [
    ("H2O", [1538.86, 3642.96, 3651.45], 3, 6.026981),
    ("CH4", [1384.87, 1384.87, 1384.87, 1556.89, 1556.89,
             3096.01, 3109.71, 3109.71, 3109.71], 3, 6.418981),
    ("CO2", [599.92, 599.92, 1425.30, 2594.02], 2, 7.130121),
    ("N2", [2437.38], 2, 4.970154),
]


@pytest.mark.parametrize("name,wavenumbers,n_rot,expected_cv", _CASES)
def test_heat_capacity_is_constant_volume_cv(name, wavenumbers, n_rot, expected_cv):
    thermo = run_thermo(
        np.array(wavenumbers),
        atoms=molecule(name),
        temperature=_T,
        pressure=101325,
        energy=0.0,
        engine="dftb+",
    )
    cv = thermo.total_heat_capacity("cal/(mol*K)")

    assert thermo._n_rot == n_rot
    # matches the independently re-derived analytic Cv, and the pinned reference
    assert cv == pytest.approx(_analytic_cv(wavenumbers, n_rot), abs=1e-3)
    assert cv == pytest.approx(expected_cv, abs=1e-3)
