"""Program header/banner."""

import sys

from beartype.typing import Any

from ThermoScreening.version import __version__


def print_header(file: Any = None) -> None:
    """
    A function to print the header
    
    The header is printed to standard error stream. 
    """

    header = rf"""
 ________   ______  
|        \ /      \  
 \$$$$$$$$|  $$$$$$\ 
   | $$   | $$___\$$ 
   | $$    \$$    \  
   | $$    _\$$$$$$\ 
   | $$   |  \__| $$ 
   | $$    \$$    $$ 
    \$$     \$$$$$$  

ThermoScreening - v{__version__}
    
"""
    if file is None:
        print(header, file=sys.stderr)
    else:
        print(header, file=file)
