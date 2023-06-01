#!/usr/bin/env python3

import logging


logging.basicConfig(filename=snakemake.log[0], filemode="w", level=logging.DEBUG)


logging.debug("Collecting combined coverage stats")


allCoverage = dict()
with open(snakemake.input[0], "r") as infh:
    infh.readline()
    for line in infh:
        l = line.strip().split()
        try:
            assert(type(allCoverage[l[1]]) is dict)
        except (AssertionError, KeyError) as err:
            allCoverage[l[1]] = {"mean":0,"median":0}
        allCoverage[l[1]]["mean"] += float(l[2])
        allCoverage[l[1]]["median"] += float(l[3])


logging.debug("Printing all sample coverage")


with open(snakemake.output.all_cov, "w") as outCov:
    outCov.write("Contig\tMean\tMedian\n")
    for contig in sorted(allCoverage.keys()):
        outCov.write("\t".join([
            contig,
            "{:.{}g}".format(allCoverage[contig]["mean"], 4),
            "{:.{}g}".format(allCoverage[contig]["median"], 4) + "\n"
        ]))
