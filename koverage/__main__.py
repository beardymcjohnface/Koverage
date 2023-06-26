"""
Entrypoint for Koverage
"""

import os
import click

from .util import (
    snake_base,
    get_version,
    default_to_output,
    copy_config,
    run_snakemake,
    OrderedCommands,
    print_citation,
)


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
            default=50,
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
            "--snake-default",
            multiple=True,
            default=[
                "--rerun-incomplete",
                "--printshellcmds",
                "--nolock",
                "--show-failed-logs",
            ],
            help="Customise Snakemake runtime args",
            show_default=True,
        ),
        click.option(
            "--pyspy",
            is_flag=True,
            show_default=True,
            default=False,
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
Required:           koverage run --assembly [file] --reads [dir]
Specify threads:    koverage run ... --threads [threads]
Disable conda:      koverage run ... --no-use-conda 
Change defaults:    koverage run ... --snake-default="-k --nolock"
Add Snakemake args: koverage run ... --dry-run --keep-going --touch
Specify targets:    koverage run ... map kmer
                    koverage run ... print_targets
Available targets:
    map             Mapping-based coverage (default)
    kmer            Kmer-based coverage. This is faster than mapping for 
                    large reference FASTAs but does not provide read 
                    counts, RPKM values etc.
    bench           A more typical approach--included for benchmarking
                    (minimap2 -> samtools sort -> coverm)
    print_targets   List available targets
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
        "args": {
            "reads": kwargs["reads"],
            "ref": kwargs["ref"],
            "minimap": kwargs["minimap"],
            "pafs": kwargs["pafs"],
            "bin_width": kwargs["bin_width"],
            "output": kwargs["output"],
            "kmer_size": kwargs["kmer_size"],
            "kmer_sample": kwargs["kmer_sample"],
            "kmer_min": kwargs["kmer_min"],
            "kmer_max": kwargs["kmer_max"],
            "log": kwargs["log"],
            "pyspy": kwargs["pyspy"],
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
    merge_config = {
        "args": {
            "reads": snake_base(os.path.join("test_data", "reads")),
            "ref": snake_base(os.path.join("test_data", "ref.fa")),
            "minimap": kwargs["minimap"],
            "pafs": kwargs["pafs"],
            "bin_width": kwargs["bin_width"],
            "output": kwargs["output"],
            "kmer_size": kwargs["kmer_size"],
            "kmer_sample": kwargs["kmer_sample"],
            "kmer_min": kwargs["kmer_min"],
            "kmer_max": kwargs["kmer_max"],
            "log": kwargs["log"],
            "pyspy": kwargs["pyspy"],
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
def config(configfile, **kwargs):
    """Copy the system default config file"""
    copy_config(configfile)


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
