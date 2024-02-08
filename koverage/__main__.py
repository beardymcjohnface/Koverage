"""
Entrypoint for Koverage
"""

import os
import click

from snaketool_utils.cli_utils import (
    OrderedCommands,
    run_snakemake,
    initialise_config,
    echo_click,
)


def snake_base(rel_path):
    """Get the filepath to a Snaketool system file (relative to __main__.py)

    Args:
        rel_path (str): Filepath relative to __main__.py

    Returns (str): Resolved filepath
    """
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), rel_path)


def get_version():
    with open(snake_base("VERSION"), "r") as f:
        version = f.readline()
    return version


def print_citation():
    with open(snake_base("CITATION"), "r") as f:
        for line in f:
            echo_click(line)


def default_to_output(ctx, param, value):
    """Callback for click options; places value in output directory unless specified"""
    if param.default == value:
        return os.path.join(ctx.params["output"], value)
    return value


def common_options(func):
    """Common command line args
    Define common command line args here, and include them with the @common_options decorator below.
    """
    options = [
        click.option(
            "--output",
            help="Output directory",
            type=click.Path(dir_okay=True, writable=True, readable=True),
            default="koverage.out",
            show_default=True,
        ),
        click.option(
            "--configfile",
            default="koverage.config.yaml",
            show_default=False,
            callback=default_to_output,
            help="Custom config file [default: (outputDir)/koverage.config.yaml]",
        ),
        click.option(
            "--system-config",
            default=snake_base(os.path.join("config", "config.yaml")),
            hidden=True,
        ),
        click.option(
            "--threads", help="Number of threads to use", default=8, show_default=True
        ),
        click.option(
            "--minimap",
            help="Minimap preset",
            default="sr",
            show_default=True,
            type=click.Choice(["map-pb", "map-ont", "map-hifi", "sr"]),
        ),
        click.option(
            "--pafs",
            is_flag=True,
            show_default=True,
            default=False,
            help="Save the (compressed) PAF files",
        ),
        click.option(
            "--bin-width",
            help="Bin width for estimating read depth variance",
            show_default=True,
            default=100,
        ),
        click.option(
            "--kmer-size", help="Size of kmers to use", show_default=True, default=25
        ),
        click.option(
            "--kmer-sample",
            help="Sample every [INT]th kmer",
            show_default=True,
            default=100,
        ),
        click.option(
            "--kmer-min",
            help="Min kmers to try to sample per contig",
            show_default=True,
            default=50,
        ),
        click.option(
            "--kmer-max",
            help="Max kmers to sample per contig",
            show_default=True,
            default=5000,
        ),
        click.option(
            "--report-max-ctg",
            help="Only include the top N contigs by coverage in the summary HMTL report (use -1 for all contigs)",
            show_default=True,
            default=1000,
        ),
        click.option(
            "--report/--no-report",
            default=True,
            help="Generate HTML summary report of coverage",
            show_default=True,
        ),
        click.option(
            "--use-conda/--no-use-conda",
            default=True,
            help="Use conda for Snakemake rules",
            show_default=True,
        ),
        click.option(
            "--conda-prefix",
            default=snake_base(os.path.join("workflow", "conda")),
            help="Custom conda env directory",
            type=click.Path(),
            show_default=False,
        ),
        click.option(
            "--workflow-profile",
            default="koverage.profile",
            show_default=False,
            callback=default_to_output,
            help="Custom config file [default: (outputDir)/koverage.profile/]",
        ),
        click.option(
            "--system-workflow-profile",
            default=snake_base(os.path.join("config", "profile", "config.yaml")),
            hidden=True,
        ),
        click.option(
            "--log",
            default="koverage.log",
            callback=default_to_output,
            hidden=True,
        ),
        click.argument("snake_args", nargs=-1),
    ]
    for option in reversed(options):
        func = option(func)
    return func


@click.group(
    cls=OrderedCommands, context_settings=dict(help_option_names=["-h", "--help"])
)
@click.version_option(get_version(), "-v", "--version", is_flag=True)
def cli():
    """Quickly get coverage statistics given reads and an assembly
    \b
    For more options, run:
    koverage command --help"""
    pass


def print_splash():
    click.echo(
        """
\b
██╗  ██╗ ██████╗ ██╗   ██╗███████╗██████╗  █████╗  ██████╗ ███████╗
██║ ██╔╝██╔═══██╗██║   ██║██╔════╝██╔══██╗██╔══██╗██╔════╝ ██╔════╝
█████╔╝ ██║   ██║██║   ██║█████╗  ██████╔╝███████║██║  ███╗█████╗  
██╔═██╗ ██║   ██║╚██╗ ██╔╝██╔══╝  ██╔══██╗██╔══██║██║   ██║██╔══╝  
██║  ██╗╚██████╔╝ ╚████╔╝ ███████╗██║  ██║██║  ██║╚██████╔╝███████╗
╚═╝  ╚═╝ ╚═════╝   ╚═══╝  ╚══════╝╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝ ╚══════╝
"""
    )


help_msg_extra = """
\b
CLUSTER EXECUTION:
koverage run ... --profile [profile]
For information on Snakemake profiles see:
https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles
\b
RUN EXAMPLES:
Required:           koverage run --ref [file] --reads [dir]
Specify threads:    koverage run ... --threads [threads]
Disable conda:      koverage run ... --no-use-conda 
Change defaults:    koverage run ... --snake-default="-k --nolock"
Add Snakemake args: koverage run ... --dry-run --keep-going --touch
Specify targets:    koverage run ... map kmer
                    koverage run ... print_targets
Available targets:
    map             Mapping-based coverage (default).
    kmer            Kmer-based coverage. Faster for large reference FASTAs 
                    but only provides depth metrics.
    coverm          Wrapper for CoverM (minimap2 -> samtools sort -> coverm).
    print_targets   List available targets.
"""


@click.command(
    epilog=help_msg_extra,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@click.option("--reads", help="Input file/directory", type=str, required=True)
@click.option("--ref", help="Input reference fasta file", type=str, required=True)
@common_options
def run(**kwargs):
    """Run Koverage"""
    # Config to add or update in configfile
    merge_config = {
        "koverage": {
            "args": kwargs
        }
    }

    # run!
    run_snakemake(
        # Full path to Snakefile
        snakefile_path=snake_base(os.path.join("workflow", "Snakefile")),
        merge_config=merge_config,
        **kwargs
    )


@click.command(
    epilog=help_msg_extra,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@common_options
def test(**kwargs):
    """Run test dataset for Koverage"""
    # Config to add or update in configfile
    kwargs["reads"] = snake_base(os.path.join("test_data", "reads"))
    kwargs["ref"] = snake_base(os.path.join("test_data", "ref.fa"))

    merge_config = {
        "koverage": {
            "args": kwargs
        }
    }

    # run!
    run_snakemake(
        # Full path to Snakefile
        snakefile_path=snake_base(os.path.join("workflow", "Snakefile")),
        merge_config=merge_config,
        **kwargs
    )


@click.command()
@common_options
def config(**kwargs):
    """Copy the system default config file"""
    initialise_config(**kwargs)


@click.command()
def citation(**kwargs):
    """Print the citation(s) for this tool"""
    print_citation()


cli.add_command(run)
cli.add_command(test)
cli.add_command(config)
cli.add_command(citation)


def main():
    print_splash()
    cli()


if __name__ == "__main__":
    main()
