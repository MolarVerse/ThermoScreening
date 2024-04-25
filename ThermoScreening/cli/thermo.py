import numpy as np
from argparse import ArgumentParser
from ThermoScreening.thermo.api import execute
import time
from ..__version__ import __version__


def parse_args():
    """
    Parse command line arguments

    Returns:
    args: argparse.Namespace: command line arguments
    """
    parser = ArgumentParser(description="ThermoScreening")
    parser.add_argument("input_file", type=str, help="Input file")
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Verbose output", default=True
    )
    args = parser.parse_args()
    return args


def main():
    """
    Main function to run the thermo cli. It parses the command line arguments
    and executes the program. Execution is timed and the time elapsed is printed.

    Example:
    python thermo.py input.in output.out -p plot.png -v -t

    Returns:
    None
    """
    parser_args = parse_args()
    input_file = parser_args.input
    verbose = parser_args.verbose

    if verbose:
        print("Input file: ", input_file)
        print("Verbose: ", verbose)

    # start timer
    start = time.time()
    # Execute
    execute(input_file)

    # end timer
    end = time.time()
    print("#####################################\n")
    print("Time elapsed: ", end - start, " s")


if __name__ == "__main__":
    main()
