#!/usr/bin/env python3

import logging


logging.basicConfig(filename=snakemake.log[0], filemode="w", level=logging.DEBUG)


logging.debug("Collecting combined coverage stats")

allCoverage = dict()
with open(snakemake.input[0], 'r') as infh:
    infh.readline()
    for line in infh:
        l = line.strip().split()
        try:
            assert(type(allCoverage[l[1]]) is dict)
        except (AssertionError, KeyError) as err:
            allCoverage[l[1]] = {'count':0,'rpm':0,'rpkm':0,'rpk':0,'tpm':0}
        allCoverage[l[1]]['count'] += int(l[2])
        allCoverage[l[1]]['rpm'] += float(l[3])
        allCoverage[l[1]]['rpkm'] += float(l[4])
        allCoverage[l[1]]['rpk'] += float(l[5])
        allCoverage[l[1]]['tpm'] += float(l[6])

logging.debug("Printing all sample coverage")

with open(snakemake.output.all_cov, 'w') as outCov:
    outCov.write("Contig\tCount\tRPM\tRPKM\tRPK\tTPM\n")
    for contig in sorted(allCoverage.keys()):
        outCov.write("\t".join([
            contig,
            str(allCoverage[contig]['count']),
            "{:.{}g}".format(allCoverage[contig]["rpm"], 4),
            "{:.{}g}".format(allCoverage[contig]["rpkm"], 4),
            "{:.{}g}".format(allCoverage[contig]["rpk"], 4),
            "{:.{}g}".format(allCoverage[contig]["tpm"], 4) + "\n"
        ]))

