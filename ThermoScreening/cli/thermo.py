import numpy as np
from argparse import ArgumentParser
from ThermoScreening.thermo.api import execute
import time

"""
Main function
"""
"""
Parser for command line arguments

"""
def parse_args():
    parser = ArgumentParser(description="ThermoScreening")
    parser.add_argument("-i", "--input", type=str, default="input.in", help="Input file")
    parser.add_argument("-o", "--output", type=str, default="output.out", help="Output file")
    parser.add_argument("-p", "--plot", type=str, default="plot.png", help="Plot file")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output", default=True)
    parser.add_argument("-t", "--test", action="store_true", help="Test run")
    args = parser.parse_args()
    return args

def main():
    print("#####################################")
    print("############ThermoScreening##########")
    print("#####################################")

    parser_args = parse_args()
    input_file = parser_args.input
    output_file = parser_args.output
    plot_file = parser_args.plot
    verbose = parser_args.verbose
    test = parser_args.test

    if verbose:
        print("Input file: ", input_file)
        print("Output file: ", output_file)
        print("Plot file: ", plot_file)
        print("Verbose: ", verbose)
        print("Test: ", test)
    
    #timer 
    start = time.time()
    # Execute
    execute(input_file)
    end = time.time()
    print("#####################################\n")
    print("Time elapsed: ", end-start, " s")


if __name__ == "__main__":
    main()