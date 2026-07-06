"""The ``thermo`` command-line entry point."""

from argparse import ArgumentParser
import sys
import time

from ThermoScreening.cli.dftb_setup import (
    check_dftb_setup,
    dftb_prefix_export,
    format_diagnostics,
    install_gbsa_param,
    install_slakos,
)
from ThermoScreening.thermo.api import execute
from ThermoScreening.thermo.screening import screen
from ThermoScreening.thermo.conformers import generate as generate_conformers, write_conformers
from ThermoScreening.version import __version__


SUBCOMMANDS = {"setup-dftb", "doctor", "screen", "conformers"}


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
        help="Download a DFTB+ Slater-Koster parameter set.",
    )
    setup_parser.add_argument(
        "--install-root",
        default=None,
        help="Directory where Slater-Koster parameter sets are installed.",
    )
    setup_parser.add_argument(
        "--parameter-set",
        default="3ob",
        choices=["3ob", "mio"],
        help="Parameter set to download (default '3ob').",
    )
    setup_parser.add_argument(
        "--url",
        default=None,
        help="Archive URL override (defaults to the release URL for the set).",
    )
    setup_parser.add_argument(
        "--solvent",
        default=None,
        help="Instead of a parameter set, download the GBSA implicit-solvation "
        "parameter file for this solvent (e.g. 'water').",
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

    screen_parser = subparsers.add_parser(
        "screen",
        help="Run thermochemistry screening over a set of molecules.",
    )
    screen_parser.add_argument(
        "source",
        help="Directory of .xyz/.gen structures, or a .csv manifest "
        "(columns: path, optional name/charge).",
    )
    screen_parser.add_argument(
        "-o", "--out", default="results",
        help="Output stem; writes <out>.csv and <out>.json. Default 'results'.",
    )
    screen_parser.add_argument(
        "--charge", type=float, default=0.0,
        help="Charge for directory input and manifest rows without a charge.",
    )
    screen_parser.add_argument(
        "--temperature", type=float, default=298.15, help="Temperature in K.",
    )
    screen_parser.add_argument(
        "--pressure", type=float, default=101325.0, help="Pressure in Pa.",
    )
    screen_parser.add_argument(
        "--directory", default="screening",
        help="Root working directory; each molecule runs in <directory>/<name>.",
    )
    screen_parser.add_argument(
        "--parameter-set", default="3ob", choices=["3ob", "mio"],
        help="Slater-Koster parameter set (default '3ob'). Selects the Hamiltonian "
        "parameters and matching spin constants.",
    )
    screen_parser.add_argument(
        "--solvent", default=None,
        help="GBSA/ALPB implicit-solvation solvent (e.g. 'water') applied to every "
        "molecule. Default gas phase. Install with 'setup-dftb --solvent <name>'.",
    )
    screen_parser.add_argument(
        "--quasi-rrho", action="store_true",
        help="Use Grimme's quasi-RRHO vibrational entropy (better for low-frequency "
        "modes) instead of the pure harmonic oscillator.",
    )
    screen_parser.add_argument(
        "--engine", default="dftb+", choices=["dftb+", "xtb", "xtb-cli"],
        help="Calculation engine (default 'dftb+'). 'xtb' uses GFN-xTB via tblite "
        "(gas-phase); 'xtb-cli' uses the native xtb binary, which also does "
        "implicit solvation of charged radicals. parameter-set applies to dftb+; "
        "solvent to dftb+ and xtb-cli.",
    )
    screen_parser.add_argument(
        "--method", default="GFN2-xTB", choices=["GFN2-xTB", "GFN1-xTB"],
        help="GFN-xTB parametrisation for the xtb engines (default 'GFN2-xTB').",
    )
    screen_parser.add_argument(
        "--resume", action="store_true",
        help="Reuse successful results from a prior <out>.json and only (re)run "
        "the missing/failed molecules.",
    )

    conf_parser = subparsers.add_parser(
        "conformers",
        help="Generate conformers from a SMILES (RDKit ETKDG) and write them as "
        ".xyz files ready to screen.",
    )
    conf_parser.add_argument("smiles", help="Molecule as a SMILES string.")
    conf_parser.add_argument(
        "-o", "--out-dir", default="conformers",
        help="Directory to write <prefix>_<i>.xyz files (default 'conformers').",
    )
    conf_parser.add_argument(
        "--max-conformers", type=int, default=10,
        help="Maximum number of conformers to embed (default 10).",
    )
    conf_parser.add_argument(
        "--energy-window", type=float, default=None,
        help="Keep only conformers within this many kcal/mol of the lowest.",
    )
    conf_parser.add_argument(
        "--no-optimize", action="store_true",
        help="Skip MMFF force-field optimisation of the embedded conformers.",
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

    if argv and argv[0] in SUBCOMMANDS:
        parser = _command_parser()
        args = parser.parse_args(argv)
    else:
        parser = _run_parser()
        args = parser.parse_args(argv)
        args.command = "run"

    return args


def run_setup_dftb(parser_args):
    """
    Download a DFTB+ Slater-Koster parameter set, or a GBSA solvent parameter.
    """

    if parser_args.solvent is not None:
        param_file = install_gbsa_param(
            parser_args.solvent,
            install_root=parser_args.install_root,
            url=parser_args.url,
            force=parser_args.force,
        )
        print("GBSA solvation parameters:", param_file)
        return 0

    parameter_dir = install_slakos(
        install_root=parser_args.install_root,
        url=parser_args.url,
        force=parser_args.force,
        parameter_set=parser_args.parameter_set,
    )

    print("Slater-Koster files: ", parameter_dir)
    print("Shell configuration:")
    print(dftb_prefix_export(parameter_dir))

    return 0


def run_doctor():
    """
    Check whether the calculation backends are available.
    """

    diagnostics = check_dftb_setup()
    print(format_diagnostics(diagnostics))

    # optional backends (xtb, tblite) do not fail the check
    return 0 if all(item.ok for item in diagnostics if not item.optional) else 1


def run_screen(parser_args):
    """
    Run thermochemistry screening over a set of molecules.
    """

    results = screen(
        parser_args.source,
        out=parser_args.out,
        charge=parser_args.charge,
        temperature=parser_args.temperature,
        pressure=parser_args.pressure,
        directory=parser_args.directory,
        parameter_set=parser_args.parameter_set,
        solvent=parser_args.solvent,
        quasi_rrho=parser_args.quasi_rrho,
        engine=parser_args.engine,
        method=parser_args.method,
        resume=parser_args.resume,
    )

    failed = sum(1 for record in results if record["status"] != "ok")
    print(f"Screened {len(results)} molecules ({failed} failed).")
    print(f"Results: {parser_args.out}.csv, {parser_args.out}.json")

    return 1 if failed else 0


def run_conformers(parser_args):
    """
    Generate conformers from a SMILES and write them as .xyz files.
    """

    try:
        conformers = generate_conformers(
            parser_args.smiles,
            max_conformers=parser_args.max_conformers,
            optimize=not parser_args.no_optimize,
            energy_window=parser_args.energy_window,
        )
    except (ValueError, ImportError) as exc:
        print(f"Conformer generation failed: {exc}", file=sys.stderr)
        return 1

    paths = write_conformers(conformers, parser_args.out_dir)
    print(f"Generated {len(paths)} conformer(s) in {parser_args.out_dir}/")
    print(f"Screen them with: thermo screen {parser_args.out_dir}")
    return 0


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

    if command == "screen":
        return run_screen(parser_args)

    if command == "conformers":
        return run_conformers(parser_args)

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
