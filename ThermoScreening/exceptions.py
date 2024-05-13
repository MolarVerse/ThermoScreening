"""
A module defining exceptions of the ThermoScreening package.
"""

class ThermoScreeningException(Exception):
    """Base class for exceptions in this package."""
    
class TSValueError(ThermoScreeningException):
    """Exception raised for errors in the input value."""
