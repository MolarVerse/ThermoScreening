"""The ``thermo`` command-line entry point."""

from argparse import ArgumentParser, REMAINDER
import sys
import time

from ThermoScreening.cli.dftb_setup import (
    check_dftb_setup,
    format_diagnostics,
    install_gbsa_param,
    install_slakos,
)
from ThermoScreening.exceptions import TSValueError
from ThermoScreening.cli.slurm import write_slurm_array_script, submit_slurm_array
from ThermoScreening.thermo.api import execute
from ThermoScreening.thermo.screening import (
    collect_shards,
    redox_screen,
    screen,
    screen_shard_directory,
    rank_by_gibbs,
)
from ThermoScreening.thermo.conformers import generate as generate_conformers, write_conformers
from ThermoScreening.version import __version__


SUBCOMMANDS = {
    "run",
    "setup-dftb",
    "doctor",
    "screen",
    "collect",
    "slurm",
    "redox",
    "conformers",
}


def _add_run_arguments(parser):
    """Add arguments for the legacy input-file workflow."""
    parser.add_argument("input_file", type=str, help="Input file")
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Print timing information."
    )


def _command_parser():
    parser = ArgumentParser(description="ThermoScreening")
    parser.add_argument("--version", action="version", version=__version__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser(
        "run",
        help="Run thermochemistry from a ThermoScreening input file.",
    )
    _add_run_arguments(run_parser)

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

    doctor_parser = subparsers.add_parser(
        "doctor",
        help="Check calculation backends and parameter configuration.",
    )
    doctor_parser.add_argument(
        "--engine",
        choices=["all", "dftb+", "xtb", "xtb-cli"],
        default="all",
        help="Backend to require. 'all' succeeds when at least one backend is ready.",
    )
    doctor_parser.add_argument(
        "--parameter-set",
        choices=["3ob", "mio"],
        default="3ob",
        help="DFTB+ parameter set to check (default '3ob').",
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
        "--dispersion", default=None, choices=["d3-bj"],
        help="Add a Grimme dispersion correction to the DFTB+ Hamiltonian "
        "('d3-bj', 3ob-recommended parameters). Default none.",
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
    screen_parser.add_argument(
        "-j", "--jobs", type=int, default=1,
        help="Number of molecules to run concurrently (default 1). Each job runs "
        "in its own process and directory.",
    )
    screen_parser.add_argument(
        "--shard-index",
        type=int,
        default=None,
        help="Zero-based cluster shard to execute; requires --shard-count.",
    )
    screen_parser.add_argument(
        "--shard-count",
        type=int,
        default=None,
        help="Total number of deterministic cluster shards.",
    )

    collect_parser = subparsers.add_parser(
        "collect",
        help="Validate and combine distributed screening shards.",
    )
    collect_parser.add_argument(
        "shard_directory",
        help="Directory containing shard-*.json and shard-*-run.json files.",
    )
    collect_parser.add_argument(
        "-o",
        "--out",
        default="results",
        help="Output stem for the combined CSV, JSON, and run metadata.",
    )

    slurm_parser = subparsers.add_parser(
        "slurm",
        help="Generate or submit a Slurm array for screen or redox.",
    )
    slurm_parser.add_argument("--tasks", type=int, required=True)
    slurm_parser.add_argument("--script", default="thermoscreening.slurm")
    slurm_parser.add_argument("--job-name", default="thermoscreening")
    slurm_parser.add_argument("--cpus-per-task", type=int, default=1)
    slurm_parser.add_argument("--time", dest="walltime", default=None)
    slurm_parser.add_argument("--mem", dest="memory", default=None)
    slurm_parser.add_argument("--partition", default=None)
    slurm_parser.add_argument("--account", default=None)
    slurm_parser.add_argument(
        "--preamble",
        default=None,
        help="Shell file inserted before the worker command for modules or variables.",
    )
    slurm_parser.add_argument(
        "--submit",
        action="store_true",
        help="Submit the array and a dependent collection job with sbatch.",
    )
    slurm_parser.add_argument(
        "command_args",
        nargs=REMAINDER,
        help="Screen or redox command after '--'.",
    )

    redox_parser = subparsers.add_parser(
        "redox",
        help="Run a three-state redox screen from one starting geometry.",
    )
    redox_parser.add_argument(
        "source",
        help="SMILES, one .xyz/.gen file, a directory, or a CSV containing "
        "a path or SMILES per row.",
    )
    redox_parser.add_argument(
        "-o",
        "--out",
        default="redox-results",
        help="Output stem for aggregate and per-state CSV/JSON results.",
    )
    redox_parser.add_argument(
        "--charge", type=int, default=None,
        help="Oxidized-state charge. Ionic SMILES are inferred; structures default to 0.",
    )
    redox_parser.add_argument(
        "--reference",
        default=None,
        help="Reference molecule name, SMILES, or structure path.",
    )
    redox_parser.add_argument(
        "--reference-e1",
        type=float,
        default=None,
        help="Measured first reduction potential of the reference in V.",
    )
    redox_parser.add_argument(
        "--reference-e2",
        type=float,
        default=None,
        help="Measured second reduction potential of the reference in V.",
    )
    redox_parser.add_argument(
        "--reference-charge",
        type=int,
        default=None,
        help="Oxidized-state charge of a separate reference, if not inferable.",
    )
    redox_parser.add_argument(
        "--potential-scale",
        default=None,
        help="Label for calibrated potentials, for example 'Fc/Fc+'.",
    )
    redox_parser.add_argument(
        "--max-conformers",
        type=int,
        default=20,
        help="SMILES conformers to embed before selecting the lowest MMFF structure.",
    )
    redox_parser.add_argument(
        "--temperature", type=float, default=298.15, help="Temperature in K.",
    )
    redox_parser.add_argument(
        "--pressure", type=float, default=101325.0, help="Pressure in Pa.",
    )
    redox_parser.add_argument(
        "--directory",
        default="redox-screening",
        help="Working directory for generated inputs and state calculations.",
    )
    redox_parser.add_argument(
        "--parameter-set", default=None, choices=["3ob", "mio"],
        help="DFTB+ Slater-Koster set; DFTB+ defaults to '3ob'.",
    )
    redox_parser.add_argument(
        "--solvent",
        default=None,
        help="Implicit-solvation solvent applied consistently to every state.",
    )
    redox_parser.add_argument(
        "--dispersion", default=None, choices=["d3-bj"],
        help="DFTB+ dispersion correction. Default none.",
    )
    redox_parser.add_argument(
        "--quasi-rrho",
        action="store_true",
        help="Use quasi-RRHO vibrational entropy for every state.",
    )
    redox_parser.add_argument(
        "--engine", default="dftb+", choices=["dftb+", "xtb", "xtb-cli"],
        help="Calculation engine (default 'dftb+').",
    )
    redox_parser.add_argument(
        "--method", default=None, choices=["GFN2-xTB", "GFN1-xTB"],
        help="GFN parametrisation; xTB engines default to 'GFN2-xTB'.",
    )
    redox_parser.add_argument(
        "--resume",
        action="store_true",
        help="Reuse only state calculations with matching input fingerprints.",
    )
    redox_parser.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=1,
        help="Number of charge-state calculations to run concurrently.",
    )
    redox_parser.add_argument(
        "--shard-index",
        type=int,
        default=None,
        help="Zero-based cluster shard to execute; requires --shard-count.",
    )
    redox_parser.add_argument(
        "--shard-count",
        type=int,
        default=None,
        help="Total number of deterministic cluster shards.",
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

    parser = _command_parser()
    if argv and argv[0] not in SUBCOMMANDS and not argv[0].startswith("-"):
        argv = ["run", *argv]
    return parser.parse_args(argv)


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
    print("The parameter set is ready for automatic discovery.")

    return 0


def run_doctor(parser_args):
    """
    Check whether the calculation backends are available.
    """

    diagnostics = check_dftb_setup(parameter_set=parser_args.parameter_set)
    print(format_diagnostics(diagnostics))

    status = {item.name: item.ok for item in diagnostics}
    ready = {
        "dftb+": all(
            status.get(name, False)
            for name in ("dftb+", "modes", "parameters", "C-C.skf")
        ),
        "xtb": status.get("tblite", False),
        "xtb-cli": status.get("xtb", False),
    }
    if parser_args.engine == "all":
        available = ", ".join(name for name, ok in ready.items() if ok) or "none"
        print(f"Usable engines: {available}")
        return 0 if any(ready.values()) else 1
    return 0 if ready[parser_args.engine] else 1


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
        dispersion=parser_args.dispersion,
        quasi_rrho=parser_args.quasi_rrho,
        engine=parser_args.engine,
        method=parser_args.method,
        resume=parser_args.resume,
        jobs=parser_args.jobs,
        shard_index=getattr(parser_args, "shard_index", None),
        shard_count=getattr(parser_args, "shard_count", None),
    )

    ranked = rank_by_gibbs(results)
    if ranked:
        print("Comparable structures ranked by Gibbs free energy:")
        for position, record in enumerate(ranked, start=1):
            print(f"  {position}. {record['name']}  G = {record['G_total_hartree']:.6f} Ha")

    failures = [record for record in results if record["status"] != "ok"]
    if failures:
        print(f"Failed ({len(failures)}):")
        for record in failures:
            print(f"  {record['name']}: {record['error']}")

    print(f"Screened {len(results)} molecules ({len(failures)} failed).")
    shard_index = getattr(parser_args, "shard_index", None)
    if shard_index is None:
        result_stem = parser_args.out
    else:
        result_stem = screen_shard_directory(parser_args.out) / f"shard-{shard_index:05d}"
    print(f"Results: {result_stem}.csv, {result_stem}.json")

    return 1 if failures else 0


def run_collect(parser_args):
    """Combine and validate distributed screening results."""
    try:
        results = collect_shards(parser_args.shard_directory, out=parser_args.out)
    except TSValueError as exc:
        print(f"Collection failed: {exc}", file=sys.stderr)
        return 1

    failures = [record for record in results if record.get("status") != "ok"]
    print(f"Collected {len(results)} molecules ({len(failures)} failed).")
    print(f"Results: {parser_args.out}.csv, {parser_args.out}.json")
    return 1 if failures else 0


def run_slurm(parser_args):
    """Generate or submit a Slurm screening array."""
    command_args = list(parser_args.command_args)
    if command_args and command_args[0] == "--":
        command_args.pop(0)
    try:
        if not command_args or command_args[0] not in {"screen", "redox"}:
            raise TSValueError("Slurm arrays support the screen and redox commands.")
        nested = parse_args(command_args)
        script = write_slurm_array_script(
            command_args,
            tasks=parser_args.tasks,
            script=parser_args.script,
            local_jobs=nested.jobs,
            cpus_per_task=parser_args.cpus_per_task,
            job_name=parser_args.job_name,
            walltime=parser_args.walltime,
            memory=parser_args.memory,
            partition=parser_args.partition,
            account=parser_args.account,
            preamble=parser_args.preamble,
        )
        print(f"Slurm script: {script}")
        if parser_args.submit:
            array_id, collector_id = submit_slurm_array(
                script,
                shard_directory=screen_shard_directory(nested.out),
                out=nested.out,
                job_name=parser_args.job_name,
                partition=parser_args.partition,
                account=parser_args.account,
            )
            print(f"Submitted array job {array_id} and collector job {collector_id}.")
        else:
            print(f"Submit: sbatch {script}")
            print(
                "Collect after completion: thermo collect "
                f"{screen_shard_directory(nested.out)} -o {nested.out}"
            )
    except (OSError, TSValueError, ValueError) as exc:
        print(f"Slurm setup failed: {exc}", file=sys.stderr)
        return 1
    return 0


def run_redox(parser_args):
    """Run the three-state redox-screening workflow."""

    try:
        results = redox_screen(
            parser_args.source,
            out=parser_args.out,
            charge=parser_args.charge,
            temperature=parser_args.temperature,
            pressure=parser_args.pressure,
            directory=parser_args.directory,
            parameter_set=parser_args.parameter_set,
            solvent=parser_args.solvent,
            dispersion=parser_args.dispersion,
            quasi_rrho=parser_args.quasi_rrho,
            engine=parser_args.engine,
            method=parser_args.method,
            resume=parser_args.resume,
            jobs=parser_args.jobs,
            reference=parser_args.reference,
            reference_e1=parser_args.reference_e1,
            reference_e2=parser_args.reference_e2,
            reference_charge=parser_args.reference_charge,
            potential_scale=parser_args.potential_scale,
            max_conformers=parser_args.max_conformers,
            shard_index=parser_args.shard_index,
            shard_count=parser_args.shard_count,
        )
    except (TSValueError, ValueError) as exc:
        print(f"Redox screen failed: {exc}", file=sys.stderr)
        return 1

    successful = [result for result in results if result["status"] == "ok"]
    if successful:
        print("Reduction potentials:")
        for result in successful:
            print(
                f"  {result['name']}: E1={result['E1_V']:.4f} V, "
                f"E2={result['E2_V']:.4f} V, E2e={result['E2e_V']:.4f} V "
                f"({result['potential_scale']})"
            )

    failures = [result for result in results if result["status"] != "ok"]
    if failures:
        print(f"Failed ({len(failures)}):")
        for result in failures:
            print(f"  {result['name']}: {result['error']}")
    print(f"Processed {len(results)} molecules ({len(failures)} failed).")
    if parser_args.shard_index is None:
        result_stem = parser_args.out
    else:
        result_stem = (
            screen_shard_directory(parser_args.out)
            / f"shard-{parser_args.shard_index:05d}"
        )
    print(f"Results: {result_stem}.csv, {result_stem}.json")
    print(f"State results: {result_stem}-states.csv, {result_stem}-states.json")
    print(f"Run metadata: {result_stem}-run.json")
    return 1 if failures else 0


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
        return run_doctor(parser_args)

    if command == "screen":
        return run_screen(parser_args)

    if command == "collect":
        return run_collect(parser_args)

    if command == "slurm":
        return run_slurm(parser_args)

    if command == "redox":
        return run_redox(parser_args)

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
