"""
This module is the main entry point for the application.

It is responsible for calling the main function in the thermo module.

Example:
    $ python -m ThermoScreening
"""
import sys

from .cli import main


if __name__ == "__main__":
    sys.exit(main())
