#!/usr/bin/env python3

import logging


logging.basicConfig(filename=snakemake.log[0], filemode="w", level=logging.DEBUG)


logging.debug("Slurping variance")


var = dict()
hitrate = dict()

with open(snakemake.input.var, 'r') as varfh:
    for line in varfh:
        l = line.strip().split()
        hitrate[l[0]] = l[1]
        var[l[0]] = l[2]


logging.debug("Reading in library size")


with open(snakemake.input.lib, 'r') as f:
    rpmscale = int(f.readline().strip()) / 1000000                    # Count up the total reads in a sample and divide that number by 1,000,000 – this is our “per million” scaling factor.


logging.debug("Parsing contig counts and calculating rpm rpkm rpk")


allRpk = list()

with open(snakemake.input.counts, 'r') as t:
    counts = dict()
    for line in t:
        l = line.strip().split()
        counts[l[0]] = dict()
        counts[l[0]]["count"] = l[2]
        counts[l[0]]["rpm"] = int(l[2]) / rpmscale              # Divide the read counts by the “per million” scaling factor. This normalizes for sequencing depth, giving you reads per million (RPM)
        counts[l[0]]["rpkm"] = int(l[2]) / int(l[1])            # Divide the RPM values by the length of the gene, in kilobases. This gives you RPKM.
        rpk = int(l[2]) / (int(l[1]) / 1000)                    # Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
        counts[l[0]]["rpk"] = rpk
        allRpk.append(rpk)
rpkscale = sum(allRpk) / len(allRpk)                            # Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.


logging.debug("Calculating TPMs and printing")


with open(snakemake.output[0], 'w') as o:
    #o.write("sample\tcontig\tCount\tRPM\tRPKM\tRPK\tTPM\tHitrate\tVariance\n")
    for contig in counts.keys():
        try:
            tpm = counts[contig]["rpk"] / rpkscale              # Divide the RPK values by the “per million” scaling factor. This gives you TPM.
        except ZeroDivisionError:
            tpm = float(0)
        o.write("\t".join([
            snakemake.wildcards.sample,
            contig,
            counts[contig]["count"],
            "{:.{}g}".format(counts[contig]["rpm"], 4),
            "{:.{}g}".format(counts[contig]["rpkm"], 4),
            "{:.{}g}".format(counts[contig]["rpk"], 4),
            "{:.{}g}".format(tpm, 4),
            hitrate[contig],
            var[contig] + "\n"
        ]))
