#!/usr/bin/env python3


"""Combine the kmer-based coverage statistics from all samples

This script will take the kmer-based coverage information from each samples' coverage file and output the population-wide counts for each contig.

- `collect_kmer_coverage_stats` - Read and add the kmer counts from all samples
- `print_kmer_coverage` - Print the combined kmer coverage statistics for all contigs
"""


import logging
import gzip
import os
import subprocess


def collect_kmer_coverage_stats(input_file):
    """Combine the kmer coverage stats for all samples.

    Args:
        input_file (str): Text TSV file (Sample\tContig\tSum\tMean\tMedian\tHitrate\tVariance)

    Returns:
        allCoverage (dict):
            - key (str): contig ID
            - value (dict):
                - sum (int): sum of kmer hits
                - mean (float): mean kmer depth
                - median (float): median kmer depth
    """
    allCoverage = {}
    with gzip.open(input_file, "rt") as infh:
        infh.readline()
        for line in infh:
            l = line.strip().split()
            try:
                assert type(allCoverage[l[1]]) is dict
            except (AssertionError, KeyError):
                allCoverage[l[1]] = {"sum": 0, "mean": 0, "median": 0}
            allCoverage[l[1]]["sum"] += float(l[2])
            allCoverage[l[1]]["mean"] += float(l[3])
            allCoverage[l[1]]["median"] += float(l[4])
    return allCoverage


def print_kmer_coverage(allCoverage, output_file, lines_per_batch=1000):
    """Print the combined kmer coverage statistics from collect_kmer_coverage_stats().

    Args:
        output_file (str): Gzipped Text TSV filepath for writing
        allCoverage (dict):
            - key (str): contig ID
            - value (dict):
                - sum (int): sum of kmer hits
                - mean (float): mean kmer depth
                - median (float): median kmer depth
        lines_per_batch (int): Number of lines to compress and write at a time
    """
    with gzip.open(output_file, "wt", compresslevel=1) as file:
        batch = ["Contig\tSum\tMean\tMedian"]
        for contig in sorted(allCoverage.keys()):
            batch.append(
                "\t".join(
                    [
                        contig,
                        "{:.{}g}".format(allCoverage[contig]["sum"], 4),
                        "{:.{}g}".format(allCoverage[contig]["mean"], 4),
                        "{:.{}g}".format(allCoverage[contig]["median"], 4),
                    ]
                )
            )
            if len(batch) >= lines_per_batch:
                file.write("\n".join(batch) + "\n")
                batch = []
        if batch:
            file.write("\n".join(batch) + "\n")


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
    allCoverage = collect_kmer_coverage_stats(input_file)
    logging.debug("Printing all sample coverage")
    print_kmer_coverage(allCoverage, output_file)


if __name__ == "__main__":
    main(
        snakemake.input[0],
        snakemake.output.all_cov,
        snakemake.log[0],
        # pyspy=snakemake.params.pyspy,
        # pyspy_svg=snakemake.log.pyspy,
    )
