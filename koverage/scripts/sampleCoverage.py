#!/usr/bin/env python3

"""Calculate the coverage stats for a sample for each contig

This script will parse the raw count summary for a sample and calculate the output coverage stats for each contig.

- `calculate_coverage_stats_from_counts` - Read in the library size and the counts from minimapWrapper.py, calculate rpm, rpkm, and rpk, write the counts
"""


import logging
import os
import subprocess
import pickle
import numpy as np


def calculate_coverage_stats_from_counts(**kwargs):
    """Read in the library size and the counts from minimapWrapper.py, calculate rpm, rpkm, and rpk, write counts for sample

    Kwargs:
        count_file (str): filepath to pickle file of contig lengths and np.array count objects
        bin_width (int): bin width size
        output_file (str): filepath of ouptut file for writing
        sample (str): sample name
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
                        "{:.{}g}".format(variances[c], 4) + "\n",
                    ]
                )
            )


def main(**kwargs):
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
    logging.basicConfig(filename=kwargs["log_file"], filemode="w", level=logging.DEBUG)
    logging.debug("Reading in library size")
    calculate_coverage_stats_from_counts(**kwargs)
    logging.debug("Calculating TPMs and printing")


if __name__ == "__main__":
    main(
        count_file=snakemake.input.counts,
        log_file=snakemake.log[0],
        output_file=snakemake.output[0],
        sample=snakemake.wildcards.sample,
        bin_width=snakemake.params.binwidth,
        # pyspy=snakemake.params.pyspy,
        # pyspy_svg=snakemake.log.pyspy,
    )
