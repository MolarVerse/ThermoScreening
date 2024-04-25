import numpy as np
import scipy.constants as const


"""
This is a class for storing physical constants.
"""

PhysicalConstants = {
        "hbar" : np.double(const.hbar),
        "h" : np.double(const.h),
        "pi" : np.double(const.pi),
        "c" : np.double(const.c),
        "kB" : np.double(const.k),
        "R" : np.double(const.gas_constant),
        "h" : np.double(const.h),
        "u" : np.double(const.u),
        "kcal" : const.calorie/1000,
        "angstrom" : np.double(const.angstrom),
        "eV" : np.double(const.electron_volt),
        "H" : const.physical_constants["Hartree energy"][0],
        "amu" : const.physical_constants["atomic unit of mass"][0],
        "bohr_radius" : 5.29177210903*10**(-1),
        "Na" : const.N_A,
        "cal" : const.calorie,
        "A"   : 10**(-10),
        "angstrom_to_meter" : 10**(-10),
        "u_to_kg" : 1.66053906660*10**(-27),  
        "kcal_to_kJ" : 4.184,
        "kcal_to_kJ_mol" : 0.239005736,
        "kcal_to_J" : 4184,
        "kcal_to_J_mol" : 239.005736,
        "HztoGHz" : 10**(-9),
        "JinHartree" : 4.3597447222071*10**(-18),
}   
