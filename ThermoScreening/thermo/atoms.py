import logging

import numpy as np
from PQAnalysis.core.atom.element import (
    atomicMasses,
    atomicNumbers,
    atomicNumbersReverse,
)

from ThermoScreening.exceptions import TSValueError
from ThermoScreening.utils.custom_logging import setup_logger
from ThermoScreening import __package_name__

class Atom:
    """
    This is a class for atom objects.

    Attributes
    ----------
    symbol : str
        The chemical symbol of the atom.
    number : int
        The atomic number of the atom.
    mass : float
        The mass of the atom.
    position : np.ndarray
        The position of the atom.
    charge : float
        The charge of the atom.

    Methods
    -------
    change_atom(symbol=None, number=None, position=None)
        Changes the atom.

    Examples
    --------
    >>> atom = Atom(symbol="H", position=[0.0, 0.0, 0.0])
    >>> atom.symbol
    'H'
    >>> atom.number
    1
    >>> atom.mass
    1.00794
    >>> atom.position
    array([0., 0., 0.])
    """
    
    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    def __init__(
        self,
        symbol: str | None = None,
        number: int | None = None,
        position: np.ndarray | None = None,
    ) -> None:
        """
        Initializes the Atom class.

        Parameters
        ----------
        symbol : str, optional, default=None
            The chemical symbol of the atom.
        number : int, optional, default=None
            The atomic number of the atom.
        position : np.ndarray, optinal, default=[0.0,0.0,0.0]
            The position of the atom.

        Raises
        ------
        TSValueError
            If neither symbol nor number is given.
            If the symbol and atomic number are not consistent.
            If the position of the atom is not given.
            If the position of the atom is not a 3D vector.
            If the atomic number is not known.
            If the chemical symbol is not known.
        """

        if symbol is None and number is None:
            self.logger.error(
                "Either symbol or number has to be given to initialize the atom.",
                exception=TSValueError
            )

        if symbol is not None and number is not None:
            if number != atomicNumbers[symbol.lower()]:
                self.logger.error(
                    "The symbol and atomic number are not consistent.",
                    exception=TSValueError
                )
        if number is not None:
            self._number = number
            try:
                self._symbol = atomic_Symbol[int(number)].capitalize()
            except:
                self.logger.error(
                    "The atomic number %s is not known." % number,
                    exception=TSValueError
                )
            self._mass = atomicMasses[self._symbol.lower()]
            self._configuration = atomicElectronConfigurations[self._symbol.lower()]
        else:
            self._symbol = symbol
            try:
                self._number = atomicNumbers[symbol.lower()]
            except:
                self.logger.error(
                    "The chemical symbol %s is not known." % symbol,
                    exception=TSValueError
                )
            self._mass = atomicMasses[symbol.lower()]
            self._configuration = atomicElectronConfigurations[symbol.lower()]

        if position is None:
            self.logger.error(
                "The position of the atom has to be given to initialize the atom.",
                exception=TSValueError
            )

        else:
            self._position = position

        if len(position) != 3:
            self.logger.error(
                "The position of the atom has to be a 3D vector.",
                exception=TSValueError
            )

    @property
    def symbol(self):
        """
        Chemical symbol of the atom.

        Returns
        -------
        str
            The chemical symbol of the atom.
        """
        return self._symbol

    @property
    def mass(self):
        """
        Mass of the atom.

        Returns
        -------
        float
            The mass of the atom.
        """
        return self._mass

    @property
    def number(self):
        """
        Atomic number of the atom.

        Returns
        -------
        int
            The atomic number of the atom.
        """
        return self._number

    @property
    def position(self):
        """
        Position of the atom.

        Returns
        -------
        np.ndarray
            The position of the atom.
        """
        return self._position

    @position.setter
    def position(self, array):
        """
        Sets the position of the atom.

        Parameters
        ----------
        value : np.ndarray
            The position of the atom.
            
        Raises
        ------
        TSValueError
            If the length of the position is not 3.
        """
        if len(array) != 3:
            self.logger.error(
                "The position of the atom has to be a 3D vector.",
                exception=TSValueError
            )
        self._position = array

    def change_atom(
        self,
        symbol: str | None = None,
        number: int | None = None,
        position: np.ndarray | None = None,
    ):
        """
        Changes the atom.
        Either the symbol or the number has to be given.
        The position can also be changed.

        Parameters
        ----------
        symbol : str, optional, default=None
            The chemical symbol of the atom.
        number : int, optional, default=None
            The atomic number of the atom.
        position : np.ndarray, optinal, default=None
            The position of the atom.

        """
        if symbol is not None:
            self._symbol = symbol
            try:
                self._number = atomicNumbers[self._symbol.lower()]
            except:
                self.logger.error(
                    "The chemical symbol %s is not known." % symbol,
                    exception=TSValueError
                )

            self._mass = atomicMasses[symbol.lower()]
            self._configuration = atomicElectronConfigurations[symbol.lower()]

        if number is not None:
            self._number = number
            try:
                self._symbol = atomic_Symbol[int(number)].capitalize()
            except:
                self.logger.error(
                    "The atomic number %s is not known." % number,
                    exception=TSValueError
                )
            self._mass = atomicMasses[self._symbol.lower()]
            self._configuration = atomicElectronConfigurations[self._symbol.lower()]
        if position is not None:
            
            self.position = position        


atomic_Symbol = atomicNumbersReverse

atomicElectronConfigurations = {
    "h": "1s1",
    "d": "1s1",
    "t": "1s1",
    "he": "1s2",
    "li": "1s2 2s1",
    "be": "1s2 2s2",
    "b": "1s2 2s2 2p1",
    "c": "1s2 2s2 2p2",
    "n": "1s2 2s2 2p3",
    "o": "1s2 2s2 2p4",
    "f": "1s2 2s2 2p5",
    "ne": "1s2 2s2 2p6",
    "na": "1s2 2s2 2p6 3s1",
    "mg": "1s2 2s2 2p6 3s2",
    "al": "1s2 2s2 2p6 3s2 3p1",
    "si": "1s2 2s2 2p6 3s2 3p2",
    "p": "1s2 2s2 2p6 3s2 3p3",
    "s": "1s2 2s2 2p6 3s2 3p4",
    "cl": "1s2 2s2 2p6 3s2 3p5",
    "ar": "1s2 2s2 2p6 3s2 3p6",
    "k": "1s2 2s2 2p6 3s2 3p6 4s1",
    "ca": "1s2 2s2 2p6 3s2 3p6 4s2",
    "sc": "1s2 2s2 2p6 3s2 3p6 4s2 3d1",
    "ti": "1s2 2s2 2p6 3s2 3p6 4s2 3d2",
    "v": "1s2 2s2 2p6 3s2 3p6 4s2 3d3",
    "cr": "1s2 2s2 2p6 3s2 3p6 4s1 3d5",
    "mn": "1s2 2s2 2p6 3s2 3p6 4s2 3d5",
    "fe": "1s2 2s2 2p6 3s2 3p6 4s2 3d6",
    "co": "1s2 2s2 2p6 3s2 3p6 4s2 3d7",
    "ni": "1s2 2s2 2p6 3s2 3p6 4s2 3d8",
    "cu": "1s2 2s2 2p6 3s2 3p6 4s1 3d10",
    "zn": "1s2 2s2 2p6 3s2 3p6 4s2 3d10",
    "ga": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p1",
    "ge": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p2",
    "as": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p3",
    "se": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p4",
    "br": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p5",
    "kr": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6",
    "rb": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s1",
    "sr": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s2",
    "y": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s2 4d1",
    "zr": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s2 4d2",
    "nb": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s1 4d4",
    "mo": "1s2 2s2 2p6 3s2 3p6 4s1 3d5",
    "tc": "1s2 2s2 2p6 3s2 3p6 4s2 3d5",
    "ru": "1s2 2s2 2p6 3s2 3p6 4s2 3d7",
    "rh": "1s2 2s2 2p6 3s2 3p6 4s2 3d8",
    "pd": "1s2 2s2 2p6 3s2 3p6 4s2 3d10",
    "ag": "1s2 2s2 2p6 3s2 3p6 4s1 3d10",
    "cd": "1s2 2s2 2p6 3s2 3p6 4s2 3d10",
    "in": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p1",
    "sn": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p2",
    "sb": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p3",
    "te": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p4",
    "i": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p5",
    "xe": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6",
    "cs": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s1",
    "ba": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s2",
    "la": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s2 4d1",
    "ce": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s2 4d2",
    "pr": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s1 4d4",
    "nd": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s1 4d5",
    "pm": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s1 4d6",
    "sm": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s1 4d7",
    "eu": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s1 4d8",
    "gd": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s1 4d10",
    "tb": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s1 4d10 5f1",
    "dy": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s1 4d10 5f2",
    "ho": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s1 4d10 5f3",
    "er": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s2 4d10 5f3",
    "tm": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s2 4d10 5f4",
    "yb": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s2 4d10 5f5",
    "lu": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s2 4d10 5f6",
    "hf": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s2 4d2",
    "ta": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s1 4d4",
    "w": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p5",
    "re": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6",
    "os": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6",
    "ir": "1s2 2s2 2p6 3s2 3p6 4s2 3d7",
    "pt": "1s2 2s2 2p6 3s2 3p6 4s2 3d9",
    "au": "1s2 2s2 2p6 3s2 3p6 4s2 3d10",
    "hg": "1s2 2s2 2p6 3s2 3p6 4s1 3d10",
    "tl": "1s2 2s2 2p6 3s2 3p6 4s2 3d10",
    "pb": "1s2 2s2 2p6 3s2 3p6 4s2 3d10",
    "bi": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p1",
    "po": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p2",
    "at": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p3",
    "rn": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p4",
    "fr": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s1",
    "ra": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s2",
    "ac": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s2 4f1",
    "th": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s2 4f2",
    "pa": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s2 4f3",
    "u": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s1 4f4",
    "np": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6",
    "pu": "1s2 2s2 2p6 3s2 3p6 4s2 3d10",
    "am": "1s2 2s2 2p6 3s2 3p6 4s2 3d10",
    "cm": "1s2 2s2 2p6 3s2 3p6 4s2",
    "bk": "1s2 2s2 2p6 3s2 3p6 4s2",
    "cf": "1s2 2s2 2p6 3s2 3p6 4s1",
    "es": "1s2 2s2 2p6 3s2 3p6",
    "fm": "1s2 2s2 2p6 3s2",
    "md": "1s2 2s2 2p6 3s2",
    "no": "1s2 2s2 2p6",
    "lr": "1s2 2s2",
    "q": "1s1",
    "x": "1s1",
    "cav": "1s1",
    "sup": "1s1",
    "dum": "1s1",
}
