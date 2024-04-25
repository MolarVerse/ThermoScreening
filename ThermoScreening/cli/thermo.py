import numpy as np
from argparse import ArgumentParser
from ThermoScreening.thermo.api import execute
import logging
import time
import sys
from ..__version__ import __version__

logger = logging.getLogger(__name__)

def parse_args():
    """
    Parse command line arguments

    Returns:
    args: argparse.Namespace: command line arguments
    """
    parser = ArgumentParser(description="ThermoScreening")
    parser.add_argument("input", type=str, help="Input file")
    # output file is optional and defaults to STDOUT
    parser.add_argument(
        "-o", "--output", type=str, default=sys.stdout, help="Output file"
    )
    parser.add_argument("-p", "--plot", type=str, default="plot.png", help="Plot file")
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Verbose output", default=True
    )
    parser.add_argument("-t", "--test", action="store_true", help="Test run")
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
    output_file = parser_args.output
    plot_file = parser_args.plot
    verbose = parser_args.verbose
    test = parser_args.test


    # use logging to print verbose output if verbose is True
    # if verbose is False, do not print verbose output
    if verbose:
        logging.basicConfig(level=logging.INFO)
        logging.info(f"Input file: {input_file}")
        logging.info(f"Output file: {output_file}")
        logging.info(f"Plot file: {plot_file}")
        logging.info(f"Verbose: {verbose}")
        logging.info(f"Test: {test}")
        logging.info(f"Version: {__version__}")
        start = time.time()

    # Execute main program
    execute(input_file)
    
    # Print time elapsed if verbose is True using logging
    if verbose:
        end = time.time()
        logging.info(f"Time elapsed: {end - start:.2f} seconds")


if __name__ == "__main__":
    main()
