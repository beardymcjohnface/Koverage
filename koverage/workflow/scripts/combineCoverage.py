#!/usr/bin/env python3

import logging


logging.basicConfig(filename=snakemake.log[0], filemode="w", level=logging.DEBUG)


logging.debug("Collecting combined coverage stats")

allCoverage = dict()
with open(snakemake.input[0], 'r') as infh:
    for line in infh:
        l = line.strip().split()
        try:
            assert(type(allCoverage[l[1]]) is dict)
        except (AssertionError, KeyError) as err:
            allCoverage[l[1]] = {'rpm':0,'rpkm':0,'rpk':0,'tpm':0}
        allCoverage[l[1]]['rpm'] += int(l[2])
        allCoverage[l[1]]['rpkm'] += int(l[3])
        allCoverage[l[1]]['rpk'] += int(l[4])
        allCoverage[l[1]]['tpm'] += int(l[5])


logging.debug("Printing all sample coverage")

with open(snakemake.output.all_cov, 'w') as outCov:
    for contig in sorted(allCoverage.keys()):
        outCov.write("\t".join([
            str(allCoverage[contig]['rpm']),
            str(allCoverage[contig]['rpkm']),
            str(allCoverage[contig]['rpk']),
            str(allCoverage[contig]['tpm']) + "\n"
        ]))

