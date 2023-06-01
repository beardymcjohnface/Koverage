#!/usr/bin/env python3

import logging
import gzip


logging.basicConfig(filename=snakemake.log[0], filemode="w", level=logging.DEBUG)


logging.debug("Collecting combined coverage stats")


allCoverage = dict()
with gzip.open(snakemake.input[0], "rt") as infh:
    infh.readline()
    for line in infh:
        l = line.strip().split()
        try:
            assert(type(allCoverage[l[1]]) is dict)
        except (AssertionError, KeyError) as err:
            allCoverage[l[1]] = {"sum":0,"mean":0,"median":0}
        allCoverage[l[1]]["sum"] += float(l[2])
        allCoverage[l[1]]["mean"] += float(l[3])
        allCoverage[l[1]]["median"] += float(l[4])


logging.debug("Printing all sample coverage")


with gzip.open(snakemake.output.all_cov, "wt", compresslevel=1) as file:
    lines_per_batch = 1000
    batch = ["Contig\tSum\tMean\tMedian"]
    for contig in sorted(allCoverage.keys()):
        batch.append("\t".join([
            contig,
            "{:.{}g}".format(allCoverage[contig]["sum"], 4),
            "{:.{}g}".format(allCoverage[contig]["mean"], 4),
            "{:.{}g}".format(allCoverage[contig]["median"], 4)
        ]))
        if len(batch) >= lines_per_batch:
            file.write('\n'.join(batch) + '\n')
            batch = []
    if batch:
        file.write('\n'.join(batch) + '\n')
