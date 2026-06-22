from argparse import ArgumentParser
import sys
import time

from ThermoScreening.cli.dftb_setup import (
    DEFAULT_SLAKO_URL,
    check_dftb_setup,
    dftb_prefix_export,
    format_diagnostics,
    install_slakos,
)
from ThermoScreening.thermo.api import execute
from ThermoScreening.version import __version__


DFTB_COMMANDS = {"setup-dftb", "doctor"}


def _run_parser():
    parser = ArgumentParser(description="ThermoScreening")
    parser.add_argument("input_file", type=str, help="Input file")
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Verbose output", default=True
    )
    return parser


def _command_parser():
    parser = ArgumentParser(description="ThermoScreening")
    subparsers = parser.add_subparsers(dest="command", required=True)

    setup_parser = subparsers.add_parser(
        "setup-dftb",
        help="Download the default DFTB+ Slater-Koster parameter set.",
    )
    setup_parser.add_argument(
        "--install-root",
        default=None,
        help="Directory where Slater-Koster parameter sets are installed.",
    )
    setup_parser.add_argument(
        "--url",
        default=DEFAULT_SLAKO_URL,
        help="Archive URL for the default Slater-Koster parameter set.",
    )
    setup_parser.add_argument(
        "--force",
        action="store_true",
        help="Download and extract even when the parameter set already exists.",
    )

    subparsers.add_parser(
        "doctor",
        help="Check DFTB+ executables and Slater-Koster parameter configuration.",
    )

    return parser


def parse_args(argv=None):
    """
    Parse command line arguments

    Returns:
    --------
    args : argparse.Namespace
        Command line arguments
    """
    if argv is None:
        argv = sys.argv[1:]

    if argv and argv[0] in DFTB_COMMANDS:
        parser = _command_parser()
        args = parser.parse_args(argv)
    else:
        parser = _run_parser()
        args = parser.parse_args(argv)
        args.command = "run"

    return args


def run_setup_dftb(parser_args):
    """
    Download the default DFTB+ Slater-Koster parameter set.
    """

    parameter_dir = install_slakos(
        install_root=parser_args.install_root,
        url=parser_args.url,
        force=parser_args.force,
    )

    print("Slater-Koster files: ", parameter_dir)
    print("Shell configuration:")
    print(dftb_prefix_export(parameter_dir))

    return 0


def run_doctor():
    """
    Check whether DFTB+ executables and parameters are available.
    """

    diagnostics = check_dftb_setup()
    print(format_diagnostics(diagnostics))

    return 0 if all(item.ok for item in diagnostics) else 1


def main():
    """
    Main function to run the thermo cli. It parses the command line arguments
    and executes the program. Execution is timed and the time elapsed is printed.

    Returns:
    --------
    None

    Examples:
    ---------
    >>> python thermo.py input_file.txt -v
    """

    parser_args = parse_args()

    command = getattr(parser_args, "command", "run")

    if command == "setup-dftb":
        return run_setup_dftb(parser_args)

    if command == "doctor":
        return run_doctor()

    input_file = parser_args.input_file
    verbose = parser_args.verbose

    if verbose:
        print("Input file: ", input_file)
        print("Verbose: ", verbose)

        # start timer
        start = time.time()

    #################################
    # Main execution of the program #
    #################################
    execute(input_file)

    if verbose:
        # end timer
        end = time.time()

        print("#####################################\n")
        print("Time elapsed: ", end - start, " s")

    return None


if __name__ == "__main__":
    main()
