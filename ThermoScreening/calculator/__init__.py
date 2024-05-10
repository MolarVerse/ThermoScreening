import os

# without dftb+ and modes, the module is useless
if os.system("which dftb+ > /dev/null") != 0:
    raise FileNotFoundError("The dftb+ executable is not found.")

if os.system("which modes > /dev/null") != 0:
    raise FileNotFoundError("The modes executable is not found.")

# import the classes
from .dftbplus import Geoopt, Hessian, Modes