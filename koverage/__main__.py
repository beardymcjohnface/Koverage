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
            default="config.yaml",
            show_default=False,
            callback=default_to_output,
            help="Custom config file [default: (outputDir)/config.yaml]",
        ),
        click.option(
            "--threads", help="Number of threads to use", default=8, show_default=True
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
Specify targets:    koverage run ... all print_targets
Available targets:
    all             Run everything (default)
    print_targets   List available targets
"""


@click.command(
    epilog=help_msg_extra,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@click.option("--reads", help="Input file/directory", type=str, required=True)
@click.option("--assembly", help="Input fasta file", type=str, required=True)
@click.option('--library', help='Library type', default='paired', show_default=True,
              type=click.Choice(['paired', 'single', 'longread']))
@click.option("--bams", is_flag=True, show_default=True, default=False, help="Save BAM files")
@click.option("--bin-width", help="Bin width for estimating read depth variance", default=50)
@common_options
def run(reads, assembly, library, bams, bin_width, output, log, **kwargs):
    """Run Koverage"""
    # Config to add or update in configfile
    merge_config = {
        "args": {
            "reads": reads,
            "assembly": assembly,
            "library": library,
            "bams": bams,
            "bin_width": bin_width,
            "output": output,
            "log": log
        }
    }

    # run!
    run_snakemake(
        # Full path to Snakefile
        snakefile_path=snake_base(os.path.join("workflow", "Snakefile")),
        merge_config=merge_config,
        log=log,
        **kwargs
    )


@click.command(
    epilog=help_msg_extra,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@click.option('--library', help='Library type', default='paired', show_default=True,
              type=click.Choice(['paired', 'single', 'longread']))
@click.option("--bams", is_flag=True, show_default=True, default=False, help="Save BAM files")
@click.option("--bin-width", help="Bin width for estimating read depth variance", default=50)
@common_options
def test(library, bams, bin_width, output, log, **kwargs):
    """Run test dataset for Koverage"""
    # Config to add or update in configfile
    merge_config = {
        "args": {
            "reads": snake_base(os.path.join("test", "reads")),
            "assembly": snake_base(os.path.join("test", "ref.fa")),
            "library": library,
            "bams": bams,
            "bin_width": bin_width,
            "output": output,
            "log": log
        }
    }

    # run!
    run_snakemake(
        # Full path to Snakefile
        snakefile_path=snake_base(os.path.join("workflow", "Snakefile")),
        merge_config=merge_config,
        log=log,
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
    cli()


if __name__ == "__main__":
    main()
