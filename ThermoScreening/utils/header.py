from ..__version__ import __version__
import sys

def print_header(file: str | None = None) -> None:
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
    if file is None:
        print(header, file=sys.stderr)
    else:
        print(header, file=file)