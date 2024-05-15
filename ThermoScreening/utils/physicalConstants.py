import numpy as np
"""
Physical constants and conversion factors that is used in the program.
The values are taken from the CODATA 2018 recommended values.

Attributes:
-----------
PhysicalConstants : dict
    A dictionary containing physical constants.

Examples:
---------
>>> from ThermoScreening.thermo import PhysicalConstants
>>> print(PhysicalConstants["hbar"])
1.0545718001391127e-34
"""
PhysicalConstants = {
    "h": 6.62607015 * 10**(-34), # J s
    "hbar": 1.054571817 * 10**(-34), # J s
    "c": 299792458, # m/s
    "kB": 1.380649 * 10**(-23), # J K^-1
    "R": 8.314462618, # J mol^-1 K^-1
    "N_A": 6.02214076 * 10**(23), # mol^-1 
    "u": 1.66053906660 * 10**(-27), # kg
    "H": 4.3597447222071 * 10**(-18), # J
    "cal": 4.184, # J
    "A": 10**(-10), # m
    "HztoGHz": 10**(-9), # GHz
    "eV": 1.602176634 * 10**(-19), # J 
    # NOT USED:
    # "pi": 3.14159265358979323846,  
    # "kcal": 4184, # J  
    # "angstrom": 10**(-10), # m
    # "bohr_radius": 5.29177210903 * 10**(-11), # m
}
