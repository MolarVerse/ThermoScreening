import unittest
import numpy as np
from unittest.mock import patch, mock_open
from  ThermoScreening.thermo.atoms import Atom 
import pytest

class TestAtoms(unittest.TestCase):

    def test_atom(self):
        symbol =  'C'
        number = 6
        mass = 12.0107
        position = np.array([0.1,0.1,0.1])
        atom = Atom(symbol=symbol,position=position)
        assert atom.symbol == symbol
        assert atom.number == number
        assert atom.mass == mass
        np.array_equal(atom.position,position)
        
        atom = Atom(number=number,position=position)
        assert atom.symbol.lower() == symbol.lower()
        assert atom.number == number
        assert atom.mass == mass
        np.testing.assert_array_equal(atom.position,position)

        atom = Atom(symbol=symbol,number=number,position=position)
        assert atom.symbol == symbol
        assert atom.number == number
        assert atom.mass == mass
        np.testing.assert_array_equal(atom.position,position)

        atom = Atom(symbol=symbol,number=number,position=position)
        assert atom.symbol == symbol
        assert atom.number == number
        assert atom.mass == mass
        np.testing.assert_array_equal(atom.position,position)

        with pytest.raises(ValueError) as exception:
            atom = Atom(position=position)
        assert str(exception.value) == "Either symbol or number has to be given to initialize the atom."

        with pytest.raises(ValueError) as exception:
            atom = Atom(symbol=symbol,number=7)
        assert str(exception.value) == "The symbol and atomic number are not consistent."

        with pytest.raises(ValueError) as exception:
            atom = Atom(symbol=symbol)
        assert str(exception.value) == "The position of the atom has to be given to initialize the atom."



        atom.change_atom(symbol="O")
        assert atom.symbol == 'O'
        assert atom.number == 8
        assert atom.mass == 15.9994

        atom.change_atom(number=6)
        assert atom.symbol == 'C'
        assert atom.number == 6
        assert atom.mass == 12.0107

        atom.change_atom(position=np.array([0.2,0.2,0.2]))
        assert atom.symbol == 'C'
        assert atom.number == 6
        assert atom.mass == 12.0107


        with pytest.raises(ValueError) as exception:
            atom.change_atom(symbol='ZZZ')
        assert str(exception.value) == "The chemical symbol ZZZ is not known."

        with pytest.raises(ValueError) as exception:
            atom.change_atom(number=999999)
        assert str(exception.value) == "The atomic number 999999 is not known."

        with pytest.raises(ValueError) as exception:
            atom = Atom(symbol=symbol,position=np.array([0.1,0.1,0.1,0.1])) 
        assert str(exception.value) == "The position of the atom has to be a 3D vector."


if __name__ == '__main__':
    unittest.main()

