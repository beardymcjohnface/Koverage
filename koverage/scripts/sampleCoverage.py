#!/usr/bin/env python3


"""Calculate the coverage stats for a sample for each contig

This script will parse the raw count summary for a sample and calculate the output coverage stats for each contig.

- `slurp_variance` - Read in the variance and hitrate stats from the variance file (contigID \t hitrate \t variance)
- `calculate_coverage_stats_from_counts` - Read in the library size and the counts from minimapWrapper.py, calculate rpm, rpkm, and rpk
- `print_coverage_stats` - Finish cov stat calculations and print output
"""


import logging
import os
import subprocess


def calculate_coverage_stats_from_counts(lib_file, count_file):
    """Read in the library size and the counts from minimapWrapper.py, calculate rpm, rpkm, and rpk.
    Returns counts dictionary of stats, and rpkscale for finishing the cov stat calcs.

    Args:
        lib_file (str): filepath for library size file (one line, one column of the number of reads)
        count_file (str): filepath for counts file (contigID \t contig length \t read counts)

    Returns:
        counts (dict):
            - Key (str): contigID
            - value (dict):
                - count (int): number of mapped reads
                - rpm (float): reads per million
                - rpkm (float): reads per kilobase million
                - rpk (float): reads per kilobase
                - mean (str): mean read depth from count_file
                - median (str): median read depth from count_file
                - hitrate (str): hitrate read from count_file
                - variance (str): variance read from count_file
        rpkscale (float): sum of rpk / 1 million
    """
    with open(lib_file, "r") as f:
        # Count up the total reads in a sample and divide that number by 1,000,000 – this is our “per million” scaling factors
        rpmscale = int(f.readline().strip()) / 1000000
    allRpk = list()
    counts = dict()
    with open(count_file, "r") as t:
        for line in t:
            l = line.strip().split()
            lenkb = int(l[1]) / 1000
            counts[l[0]] = dict()
            counts[l[0]]["count"] = l[2]
            # Divide the read counts by the “per million” scaling factor. This normalizes for sequencing depth, giving you reads per million (RPM)
            counts[l[0]]["rpm"] = int(l[2]) / rpmscale
            # Divide the RPM values by the length of the gene, in kilobases. This gives you RPKM.
            counts[l[0]]["rpkm"] = counts[l[0]]["rpm"] / lenkb
            # Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
            rpk = int(l[2]) / lenkb
            counts[l[0]]["rpk"] = rpk
            allRpk.append(rpk)
            counts[l[0]]["mean"] = l[3]
            counts[l[0]]["median"] = l[4]
            counts[l[0]]["hitrate"] = l[5]
            counts[l[0]]["variance"] = l[6]
    # Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
    rpkscale = sum(allRpk) / 1000000
    return counts, rpkscale


def print_coverage_stats(**kwargs):
    """Take the counts, and rpkscale from calculate_coverage_stats_from_counts;
    Take the variance and hitrate from slurp_variance;
    print the sample output coverage stats
    output format = sample \t contig \t Count \t RPM \t RPKM \t RPK \t TPM \t Mean \t Median \t Hitrate \t Variance

    Args:
        **kwargs (dict):
            - counts (dict):
                - Key (str): contigID
                - value (dict):
                    - count (int): number of mapped reads
                    - rpm (float): reads per million
                    - rpkm (float): reads per kilobase million
                    - rpk (float): reads per kilobase
            - sample (str): sample name
            - mean (dict):
                - key (str): contigID
                - value (float): mean
            - median (dict):
                - key (str): contigID
                - value (float): median
            - variance (dict):
                - key (str): contigID
                - value (float): variance
            - hitrate (dict):
                - key (str): contigID
                - value (float): hitrate
            - rpkscale (float): sum of rpk / 1 million

    """
    with open(kwargs["output_file"], "w") as o:
        for contig in kwargs["counts"].keys():
            try:
                # Divide the RPK values by the “per million” scaling factor. This gives you TPM.
                tpm = kwargs["counts"][contig]["rpk"] / kwargs["rpkscale"]
            except ZeroDivisionError:
                tpm = float(0)
            o.write(
                "\t".join(
                    [
                        kwargs["sample"],
                        contig,
                        kwargs["counts"][contig]["count"],
                        "{:.{}g}".format(kwargs["counts"][contig]["rpm"], 4),
                        "{:.{}g}".format(kwargs["counts"][contig]["rpkm"], 4),
                        "{:.{}g}".format(kwargs["counts"][contig]["rpk"], 4),
                        "{:.{}g}".format(tpm, 4),
                        kwargs["counts"][contig]["mean"],
                        kwargs["counts"][contig]["median"],
                        kwargs["counts"][contig]["hitrate"],
                        kwargs["counts"][contig]["variance"] + "\n",
                    ]
                )
            )


def main(**kwargs):
    if kwargs["pyspy"]:
        subprocess.Popen(
            [
                "py-spy",
                "record",
                "-s",
                "-o",
                kwargs["pyspy_svg"],
                "--pid",
                str(os.getpid()),
            ]
        )
    logging.basicConfig(filename=kwargs["log_file"], filemode="w", level=logging.DEBUG)
    logging.debug("Reading in library size")
    counts, rpkscale = calculate_coverage_stats_from_counts(
        kwargs["lib_file"], kwargs["count_file"]
    )
    logging.debug("Calculating TPMs and printing")
    print_coverage_stats(
        output_file=kwargs["output_file"],
        counts=counts,
        rpkscale=rpkscale,
        sample=kwargs["sample"],
    )


if __name__ == "__main__":
    main(
        lib_file=snakemake.input.lib,
        count_file=snakemake.input.counts,
        log_file=snakemake.log[0],
        output_file=snakemake.output[0],
        sample=snakemake.wildcards.sample,
        pyspy=snakemake.params.pyspy,
        pyspy_svg=snakemake.log.pyspy,
    )
