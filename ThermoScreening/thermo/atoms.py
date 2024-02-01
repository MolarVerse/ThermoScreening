import numpy as np  
import scipy.constants as const
    
"""
This is a class for atom
"""
class Atom:
 

    def __init__(self,
                symbol    :  str | None = None,
                number    : int | None = None,
                configuration    : str | None = None,
                position  : np.ndarray | None = None,
                ) -> None:

        """
        Initializes the Atom class.
        
        Parameters
        ----------
        symbol : str, optional, default=None
            The chemical symbol of the atom.
        mass : float, optional, default=None
            The mass of the atom.
        number : int, optional, default=None
            The atomic number of the atom.
        configuration : str, optional, default=None
            The electron configuration of the atom.
        position : np.ndarray, optinal, default=[0.0,0.0,0.0]
            The position of the atom.

        """
       
        if symbol is None and number is None:
            raise ValueError("Either symbol or number has to be given to initialize the atom.")
        
        if symbol is not None and number is not None: 
            if number != atomicNumbers[symbol]:
                raise ValueError("The symbol and number are not consistent.")  
        if number is not None:
            self._number = number
            self._symbol = atomicNumbers[int(number)]
            self._mass = atomicMasses[symbol.lower()]
            self._configuration = atomicElectronConfigurations[symbol.lower()]
        else:
            self._symbol = symbol
            self._number = atomicNumbers[symbol.lower()]
            self._mass = atomicMasses[symbol.lower()]
            self._configuration = atomicElectronConfigurations[symbol.lower()]


        
        if position is None:
            self._position = np.zeros(3)
        else:
            self._position = position

        if len(position) != 3:
            raise ValueError("The length of position must be 3.")
       
    
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
    def charge(self):
        """
        Charge of the atom.

        Returns
        -------
        float
            The charge of the atom.
        """
        return self._charge
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
        """
        if len(array) != 3:
            raise ValueError("The length of position must be 3.")
        self._position = array
    
    @charge.setter
    def charge(self, value):
        """
        Sets the charge of the atom.

        Parameters
        ----------
        value : float
            The charge of the atom.
        """
        self._charge = value
    
    def change_atom(self, _symbol=None, _mass=None, _number=None, _charge=None, _position=None):
        """
        Changes the atom to another atom.
        """
        self = Atom(_symbol=_symbol, _mass=_mass, _number=_number, _charge=_charge, _position=_position)






atomicMasses = {"h":    1.00794,    "d":    2.014101778,   "t":    3.0160492675,
                "he":   4.002602,   "li":   6.941,         "be":   9.012182,
                "b":   10.811,      "c":   12.0107,        "n":   14.0067,
                "o":   15.9994,     "f":   18.9984032,     "ne":  20.1797,
                "na":  22.989770,   "mg":  24.3050,        "al":  26.981538,
                "si":  28.0855,     "p":   30.973761,      "s":   32.065,
                "cl":  35.453,      "ar":  39.948,         "k":   39.0983,
                "ca":  40.078,      "sc":  44.955910,      "ti":  47.880,
                "v":   50.9415,     "cr":  51.9961,        "mn":  54.938049,
                "fe":  55.845,      "co":  58.933200,      "ni":  58.6934,
                "cu":  63.546,      "zn":  65.399,         "ga":  69.723,
                "ge":  72.64,       "as":  74.92160,       "se":  78.96,
                "br":  79.904,      "kr":  83.798,         "rb":  85.4678,
                "sr":  87.62,       "y":   88.90585,       "zr":  91.224,
                "nb":  92.90638,    "mo":  95.94,          "tc":  98.9063,
                "ru": 101.07,       "rh": 102.9055,        "pd": 106.42,
                "ag": 107.8682,     "cd": 112.411,         "in": 114.818,
                "sn": 118.71,       "sb": 121.76,          "te": 127.6,
                "i":  126.90447,    "xe": 131.293,         "cs": 132.90546,
                "ba": 137.327,      "la": 138.9055,        "ce": 140.116,
                "pr": 140.90765,    "nd": 144.24,          "pm": 146.9151,
                "sm": 150.36,       "eu": 151.964,         "gd": 157.25,
                "tb": 158.92534,    "dy": 162.5,           "ho": 164.93032,    
                "er": 167.259,      "tm": 168.93421,       "yb": 173.04,       
                "lu": 174.967,      "hf": 178.49,          "ta": 180.9479,     
                "w":  183.84,       "re": 186.207,         "os": 190.23,       
                "ir": 192.217,      "pt": 195.078,         "au": 196.96655,    
                "hg": 200.59,       "tl": 204.3833,        "pb": 207.2,        
                "bi": 208.98038,    "po": 208.9824,        "at": 209.9871,     
                "rn": 222.0176,     "fr": 223.0197,        "ra": 226.0254,     
                "ac": 227.0278,     "th": 232.0381,        "pa": 231.03588,    
                "u":  238.0289,     "np": 237.0482,        "pu": 244.0642,     
                "am": 243.0614,     "cm": 247.0703,        "bk": 247.0703,     
                "cf": 251.0796,     "es": 252.0829,        "fm": 257.0951,     
                "md": 258.0986,     "no": 259.1009,        "lr": 260.1053,
                "q":  999.00000,    "x":  999.00000,       "cav": 1000.00000, 
                "sup": 1000000.0, 
                "dum": 1.0}

atomicNumbers = {"h":     1,  "d":     1,  "t":    1,
                 "he":    2,  "li":    3,  "be":   4,
                 "b":     5,  "c":     6,  "n":    7,
                 "o":     8,  "f":     9,  "ne":  10,
                 "na":   11,  "mg":   12,  "al":  13,
                 "si":   14,  "p":    15,  "s":   16,
                 "cl":   17,  "ar":   18,  "k":   19,
                 "ca":   20,  "sc":   21,  "ti":  22,
                 "v":    23,  "cr":   24,  "mn":  25,
                 "fe":   26,  "co":   27,  "ni":  28,
                 "cu":   29,  "zn":   30,  "ga":  31,
                 "ge":   32,  "as":   33,  "se":  34,
                 "br":   35,  "kr":   36,  "rb":  37,
                 "sr":   38,  "y":    39,  "zr":  40,
                 "nb":   41,  "mo":   42,  "tc":  43,
                 "ru":   44,  "rh":   45,  "pd":  46,
                 "ag":   47,  "cd":   48,  "in":  49,
                 "sn":   50,  "sb":   51,  "te":  52,
                 "i":    53,  "xe":   54,  "cs":  55,
                 "ba":   56,  "la":   57,  "ce":  58,
                 "pr":   59,  "nd":   60,  "pm":  61,
                 "sm":   62,  "eu":   63,  "gd":  64,
                 "tb":   65,  "dy":   66,  "ho":  67,
                 "er":   68,  "tm":   69,  "yb":  70,
                 "lu":   71,  "hf":   72,  "ta":  73,
                 "w":    74,  "re":   75,  "os":  76,
                 "ir":   77,  "pt":   78,  "au":  79,
                 "hg":   80,  "tl":   81,  "pb":  82,
                 "bi":   83,  "po":   84,  "at":  85,
                 "rn":   86,  "fr":   87,  "ra":  88,
                 "ac":   89,  "th":   90,  "pa":  91,
                 "u":    92,  "np":   93,  "pu":  94,
                 "am":   95,  "cm":   96,  "bk":  97,
                 "cf":   98,  "es":   99,  "fm": 100,
                 "md":  101,  "no":  102,  "lr": 103, 
                 "q":   999,  "x":   999,  "cav": 1000, 
                 "sup": 1000000, 
                 "dum": 1}

atomicElectronConfigurations = {"h":     "1s1",  "d":     "1s1",  "t":    "1s1",
                                "he":    "1s2",  "li":    "1s2 2s1",  "be":   "1s2 2s2",
                                "b":     "1s2 2s2 2p1",  "c":     "1s2 2s2 2p2",  "n":    "1s2 2s2 2p3",
                                "o":     "1s2 2s2 2p4",  "f":     "1s2 2s2 2p5",  "ne":  "1s2 2s2 2p6",
                                "na":   "1s2 2s2 2p6 3s1",  "mg":   "1s2 2s2 2p6 3s2",  "al":  "1s2 2s2 2p6 3s2 3p1",
                                "si":   "1s2 2s2 2p6 3s2 3p2",  "p":    "1s2 2s2 2p6 3s2 3p3",  "s":   "1s2 2s2 2p6 3s2 3p4",
                                "cl":   "1s2 2s2 2p6 3s2 3p5",  "ar":   "1s2 2s2 2p6 3s2 3p6",  "k":   "1s2 2s2 2p6 3s2 3p6 4s1",
                                "ca":   "1s2 2s2 2p6 3s2 3p6 4s2",  "sc":   "1s2 2s2 2p6 3s2 3p6 4s2 3d1",  "ti":  "1s2 2s2 2p6 3s2 3p6 4s2 3d2",
                                "v":    "1s2 2s2 2p6 3s2 3p6 4s2 3d3",  "cr":   "1s2 2s2 2p6 3s2 3p6 4s1 3d5",  "mn":  "1s2 2s2 2p6 3s2 3p6 4s2 3d5",
                                "fe":   "1s2 2s2 2p6 3s2 3p6 4s2 3d6",  "co":   "1s2 2s2 2p6 3s2 3p6 4s2 3d7",  "ni":  "1s2 2s2 2p6 3s2 3p6 4s2 3d8",
                                "cu":   "1s2 2s2 2p6 3s2 3p6 4s1 3d10",  "zn":   "1s2 2s2 2p6 3s2 3p6 4s2 3d10",  "ga":  "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p1",
                                "ge":   "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p2",  "as":   "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p3",  "se":  "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p4",
                                "br":   "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p5",  "kr":   "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6",  "rb":  "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s1",
                                "sr":   "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s2",  "y":   "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s2 4d1",  "zr":  "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s2 4d2",
                                "nb":   "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s1 4d4",  "mo":   "1s2 2s2 2p6 3s2 3p6 4s1 3d5",  "tc":  "1s2 2s2 2p6 3s2 3p6 4s2 3d5",
                                "ru":   "1s2 2s2 2p6 3s2 3p6 4s2 3d7",  "rh":   "1s2 2s2 2p6 3s2 3p6 4s2 3d8",  "pd":  "1s2 2s2 2p6 3s2 3p6 4s2 3d10",
                                "ag":   "1s2 2s2 2p6 3s2 3p6 4s1 3d10",  "cd":   "1s2 2s2 2p6 3s2 3p6 4s2 3d10",  "in":  "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p1",
                                "sn":   "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p2",  "sb":   "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p3",  "te":  "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p4",
                                "i":    "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p5",  "xe":   "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6",  "cs":  "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s1",
                                "ba":   "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s2",  "la":   "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s2 4d1",  "ce":  "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s2 4d2",
                                "pr":   "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s1 4d4",  "nd":   "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s1 4d5",  "pm":  "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s1 4d6",
                                "sm":   "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s1 4d7",  "eu":   "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s1 4d8",  "gd":  "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s1 4d10",
                                "tb":   "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s1 4d10 5f1",  "dy":   "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s1 4d10 5f2",  "ho":  "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s1 4d10 5f3",
                                "er":   "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s2 4d10 5f3",  "tm":   "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s2 4d10 5f4",  "yb":  "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s2 4d10 5f5",
                                "lu":   "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s2 4d10 5f6",  "hf":   "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s2 4d2",  "ta":  "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s1 4d4",
                                "w":    "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p5",  "re":   "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6",  "os":  "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6",
                                "ir":   "1s2 2s2 2p6 3s2 3p6 4s2 3d7",  "pt":   "1s2 2s2 2p6 3s2 3p6 4s2 3d9",  "au":  "1s2 2s2 2p6 3s2 3p6 4s2 3d10",
                                "hg":   "1s2 2s2 2p6 3s2 3p6 4s1 3d10",  "tl":   "1s2 2s2 2p6 3s2 3p6 4s2 3d10",  "pb":  "1s2 2s2 2p6 3s2 3p6 4s2 3d10",
                                "bi":   "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p1",  "po":   "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p2",  "at":  "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p3",
                                "rn":   "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p4",  "fr":   "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s1",  "ra":  "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s2",
                                "ac":   "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s2 4f1",  "th":   "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s2 4f2",  "pa":  "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s2 4f3",
                                "u":    "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s1 4f4",  "np":   "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6",  "pu":  "1s2 2s2 2p6 3s2 3p6 4s2 3d10",
                                "am":   "1s2 2s2 2p6 3s2 3p6 4s2 3d10",  "cm":   "1s2 2s2 2p6 3s2 3p6 4s2",  "bk":  "1s2 2s2 2p6 3s2 3p6 4s2",
                                "cf":   "1s2 2s2 2p6 3s2 3p6 4s1",  "es":   "1s2 2s2 2p6 3s2 3p6",  "fm":  "1s2 2s2 2p6 3s2",
                                "md":   "1s2 2s2 2p6 3s2",  "no":   "1s2 2s2 2p6",  "lr":  "1s2 2s2",
                                "q":   "1s1",  "x":   "1s1",  "cav": "1s1",
                                "sup": "1s1",
                                "dum": "1s1"}
