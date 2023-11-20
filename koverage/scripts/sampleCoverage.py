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
import pickle
import numpy as np


def calculate_coverage_stats_from_counts(**kwargs):
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

    with open(kwargs["count_file"], "rb") as handle:
        # [[ctg,len],...]
        contig_lengths = pickle.load(handle)
        # row = ctg, col = bin
        contig_bin_counts = pickle.load(handle)

    contiglens = np.array([row[1] for row in contig_lengths], dtype=np.int32)
    contiglenkb = contiglens / 1000

    sums = np.sum(contig_bin_counts, axis=1)
    means = np.zeros(len(contiglens))
    medians = np.zeros(len(contiglens))
    hitrates = np.zeros(len(contiglens))
    variances = np.zeros(len(contiglens))

    for c in range(len(contiglens)):
        maxBin = int(contiglens[c] / kwargs["bin_width"])
        means[c] = np.mean(contig_bin_counts[c, 0:maxBin])
        medians[c] = np.median(contig_bin_counts[c, 0:maxBin])
        hitrates[c] = np.sum(contig_bin_counts[c, 0:maxBin] != 0) / maxBin
        variances[c] = np.var(contig_bin_counts[c, 0:maxBin], ddof=1)

    lib_size = np.sum(sums)
    rpmscale = lib_size / 1000000

    rpm = sums / rpmscale
    rpkm = rpm / contiglenkb
    rpk = sums / contiglenkb
    rpk_scale = np.sum(rpk) / 1000000
    if rpk_scale > 0:
        tpm = rpk / rpk_scale
    else:
        tpm = np.zeros(len(contiglens))

    variances = np.nan_to_num(variances, nan=0)
    rpm = np.nan_to_num(rpm, nan=0)
    rpkm = np.nan_to_num(rpkm, nan=0)

    with open(kwargs["output_file"], "w") as o:
        for c in range(len(contiglens)):
            o.write(
                "\t".join(
                    [
                        kwargs["sample"],
                        contig_lengths[c][0],
                        "{:d}".format(int(sums[c])),
                        "{:.{}g}".format(rpm[c], 4),
                        "{:.{}g}".format(rpkm[c], 4),
                        "{:.{}g}".format(rpk[c], 4),
                        "{:.{}g}".format(tpm[c], 4),
                        "{:.{}g}".format(means[c], 4),
                        "{:.{}g}".format(medians[c], 4),
                        "{:.{}g}".format(hitrates[c], 4),
                        "{:.{}g}".format(variances[c], 4) + "\n"
                    ]
                )
            )


# def print_coverage_stats(**kwargs):
#     """Take the counts, and rpkscale from calculate_coverage_stats_from_counts;
#     Take the variance and hitrate from slurp_variance;
#     print the sample output coverage stats
#     output format = sample \t contig \t Count \t RPM \t RPKM \t RPK \t TPM \t Mean \t Median \t Hitrate \t Variance
#
#     Args:
#         **kwargs (dict):
#             - counts (dict):
#                 - Key (str): contigID
#                 - value (dict):
#                     - count (int): number of mapped reads
#                     - rpm (float): reads per million
#                     - rpkm (float): reads per kilobase million
#                     - rpk (float): reads per kilobase
#             - sample (str): sample name
#             - mean (dict):
#                 - key (str): contigID
#                 - value (float): mean
#             - median (dict):
#                 - key (str): contigID
#                 - value (float): median
#             - variance (dict):
#                 - key (str): contigID
#                 - value (float): variance
#             - hitrate (dict):
#                 - key (str): contigID
#                 - value (float): hitrate
#             - rpkscale (float): sum of rpk / 1 million
#
#     """
#     with open(kwargs["output_file"], "w") as o:
#         for contig in kwargs["counts"].keys():
#             try:
#                 # Divide the RPK values by the “per million” scaling factor. This gives you TPM.
#                 tpm = kwargs["counts"][contig]["rpk"] / kwargs["rpkscale"]
#             except ZeroDivisionError:
#                 tpm = float(0)
#             o.write(
#                 "\t".join(
#                     [
#                         kwargs["sample"],
#                         contig,
#                         kwargs["counts"][contig]["count"],
#                         "{:.{}g}".format(kwargs["counts"][contig]["rpm"], 4),
#                         "{:.{}g}".format(kwargs["counts"][contig]["rpkm"], 4),
#                         "{:.{}g}".format(kwargs["counts"][contig]["rpk"], 4),
#                         "{:.{}g}".format(tpm, 4),
#                         kwargs["counts"][contig]["mean"],
#                         kwargs["counts"][contig]["median"],
#                         kwargs["counts"][contig]["hitrate"],
#                         kwargs["counts"][contig]["variance"] + "\n",
#                     ]
#                 )
#             )


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
    calculate_coverage_stats_from_counts(**kwargs)
    logging.debug("Calculating TPMs and printing")
    # print_coverage_stats(
    #     output_file=kwargs["output_file"],
    #     counts=counts,
    #     rpkscale=rpkscale,
    #     sample=kwargs["sample"],
    # )


if __name__ == "__main__":
    main(
        # lib_file=snakemake.input.lib,
        count_file=snakemake.input.counts,
        log_file=snakemake.log[0],
        output_file=snakemake.output[0],
        sample=snakemake.wildcards.sample,
        bin_width=snakemake.params.binwidth,
        pyspy=snakemake.params.pyspy,
        pyspy_svg=snakemake.log.pyspy,
    )
