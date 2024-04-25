from ..__version__ import __version__
import sys

def print_header() -> None:
    """
    A function to print the header
    
    The header is printed to standard error stream. 
    """

    header = f"""
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
    print(header, file=sys.stderr)