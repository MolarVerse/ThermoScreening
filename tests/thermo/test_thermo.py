import unittest
import numpy as np
from ThermoScreening.thermo.system import System
from ThermoScreening.thermo.atoms import Atom
from ThermoScreening.thermo.cell import Cell
from ThermoScreening.thermo.thermo import Thermo
import pytest
import sys, os

class TestThermo(unittest.TestCase):
    def test_thermo(self):
        # test all parameters in the class thermo   
        
        path = os.getcwd()
        #contains tests delete from path
        if path.endswith("tests"):
            path = path[:-5]
        frequencies = np.loadtxt(path+'/tests/data/thermo/frequency.txt', usecols=1)
        cell = None
        atoms = np. loadtxt(path+'/tests/data/thermo/geo_opt.xyz', usecols=0, dtype=str, skiprows=2)
        coord = np. loadtxt(path+'/tests/data/thermo/geo_opt.xyz', usecols=(1,2,3), skiprows=2)
        atom_list = []
        for i in range(len(atoms)):
            atom_list.append(Atom(symbol=atoms[i], position=coord[i]))

        
        system = System(atoms=atom_list,cell=cell,periodicity=False,solvation=None,solvent=None,charge=0,electronic_energy=-33.6052447996,vibrational_frequencies=frequencies)

        thermo = Thermo(system=system,temperature=298.15,pressure=101325,engine="DFTB+")

        assert thermo._temperature == 298.15

        assert thermo._pressure == 101325

        assert thermo._engine == "DFTB+"

        thermo.run()

        assert thermo.electronic_energy() == -33.6052447996

        assert thermo.total_energy("H") ==  0.18747059799192353

        assert thermo.total_enthalpy("H") ==  0.1884147828595504

        assert thermo.total_entropy("cal/(mol*K)") ==  107.5750500838071

        assert thermo.total_gibbs_free_energy("H") ==  0.13635822039734372 

        assert thermo.total_heat_capacity("cal/(mol*K)") ==   46.853366199557115

        assert thermo.total_EeZPE() ==  -33.42938538587388 

        assert thermo.total_EeEtot() ==  -33.417774201608076 

        assert thermo.total_EeHtot() == -33.41683001674045 

        assert thermo.total_EeGtot() == -33.46888657920266  # -33.46794021246214 
      

        





if __name__ == '__main__':
    unittest.main()