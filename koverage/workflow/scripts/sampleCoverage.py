#!/usr/bin/env python3

import os
import atexit
import logging


def exitLogCleanup(*args):
    """Cleanup the logging file(s) prior to exiting"""
    for logFile in args:
        os.unlink(logFile)
    return None


atexit.register(exitLogCleanup, snakemake.log[0])
logging.basicConfig(filename=snakemake.log[0], filemode="w", level=logging.DEBUG)

logging.debug("Calculating sample RPKM etc")

logging.debug("Reading in library size")
with open(snakemake.input.r1, 'r') as f:
    rpmscale = int(f.readline().strip()) / 1000000                    # Count up the total reads in a sample and divide that number by 1,000,000 – this is our “per million” scaling factor.

logging.debug("Parsing contig counts and calculating rpm rpkm rpk")
allRpk = []
with open(input.tsv, 'r') as t:
    counts = dict()
    for line in t:
        l = line.strip().split()
        counts[l[0]] = dict()
        counts[l[0]]["rpm"] = int(l[2]) / rpmscale              # Divide the read counts by the “per million” scaling factor. This normalizes for sequencing depth, giving you reads per million (RPM)
        counts[l[0]]["rpkm"] = int(l[2]) / int(l[1])            # Divide the RPM values by the length of the gene, in kilobases. This gives you RPKM.
        rpk = int(l[1]) / (int(l[1]) / 1000)                    # Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
        counts[l[0]]["rpk"] = rpk
        allRpk.append(rpk)
rpkscale = sum(allRpk) / len(allRpk)                            # Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.

logging.debug("Calculating TPMs and printing")
with open(snakemake.output[0], 'w') as o:
    #o.write("sample\tcontig\tRPM\tRPKM\tRPK\tTPM\n")
    for contig in counts.keys():
        tpm = counts[contig]["rpk"] / rpkscale                  # Divide the RPK values by the “per million” scaling factor. This gives you TPM.
        o.write("\t".join([
            snakemake.wildcards.sample,
            contig,
            contig['rpm'],
            contig['rpkm'],
            contig['rpk'],
            tpm + "\n"
        ]))


