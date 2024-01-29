#!/usr/bin/env python3


"""Combine the coverage statistics from all samples

This script will take the coverage information from each samples' coverage file and output the population-wide counts for each contig.

- `collect_coverage_stats` - Read and add the counts from all samples
- `print_sample_coverage` - Print the combined coverage statistics for all contigs
"""


import logging
import os
import subprocess


def collect_coverage_stats(input_file):
    """Combine the mapped coverage stats for all samples.

    Args:
        input_file (str): Text TSV file (Sample\tContig\tCount\tRPM\tRPKM\tRPK\tTPM\tHitrate\tVariance)

    Returns:
        all_coverage (dict):
            - key (str): contig ID
            - value (dict):
                - count (int): number of reads
                - rpm (float): reads per million
                - rpkm (float): reads per kilobase million
                - tpm (float): transcripts per million
    """
    all_coverage = {}
    with open(input_file, "r") as infh:
        infh.readline()
        for line in infh:
            l = line.strip().split()
            try:
                assert type(all_coverage[l[1]]) is dict
            except (AssertionError, KeyError):
                all_coverage[l[1]] = {
                    "count": 0,
                    "rpm": 0,
                    "rpkm": 0,
                    "rpk": 0,
                    "tpm": 0,
                }
            all_coverage[l[1]]["count"] += int(l[2])
            all_coverage[l[1]]["rpm"] += float(l[3])
            all_coverage[l[1]]["rpkm"] += float(l[4])
            all_coverage[l[1]]["rpk"] += float(l[5])
            all_coverage[l[1]]["tpm"] += float(l[6])
    return all_coverage


def print_sample_coverage(output_file, all_coverage):
    """Print the combined coverage statistics from collect_coverage_stats().

    Args:
        output_file (str): Text TSV filepath for writing
        all_coverage (dict):
            - key (str): contig ID
            - value (dict):
                - count (int): number of reads
                - rpm (float): reads per million
                - rpkm (float): reads per kilobase million
                - tpm (float): transcripts per million
    """
    with open(output_file, "w") as outCov:
        outCov.write("Contig\tCount\tRPM\tRPKM\tRPK\tTPM\n")
        for contig in sorted(all_coverage.keys()):
            outCov.write(
                "\t".join(
                    [
                        contig,
                        str(all_coverage[contig]["count"]),
                        "{:.{}g}".format(all_coverage[contig]["rpm"], 4),
                        "{:.{}g}".format(all_coverage[contig]["rpkm"], 4),
                        "{:.{}g}".format(all_coverage[contig]["rpk"], 4),
                        "{:.{}g}".format(all_coverage[contig]["tpm"], 4) + "\n",
                    ]
                )
            )


def main(input_file, output_file, log_file, **kwargs):
    # if kwargs["pyspy"]:
    #     subprocess.Popen(
    #         [
    #             "py-spy",
    #             "record",
    #             "-s",
    #             "-o",
    #             kwargs["pyspy_svg"],
    #             "--pid",
    #             str(os.getpid()),
    #         ]
    #     )
    logging.basicConfig(filename=log_file, filemode="w", level=logging.DEBUG)
    logging.debug("Collecting combined coverage stats")
    all_coverage = collect_coverage_stats(input_file)
    logging.debug("Printing all sample coverage")
    print_sample_coverage(output_file, all_coverage)


if __name__ == "__main__":
    main(
        snakemake.input[0],
        snakemake.output.all_cov,
        snakemake.log[0],
        # pyspy=snakemake.params.pyspy,
        # pyspy_svg=snakemake.log.pyspy,
    )
