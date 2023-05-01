#!/usr/bin/env python3

import logging
from scipy.stats import kurtosis


logging.basicConfig(filename=snakemake.log[0], filemode="w", level=logging.DEBUG)


logging.debug("Collecting combined coverage stats")


def dumpContig(ctg, dict, fh):
    for depth in sorted(dict.keys()):
        fh.write(f"{ctg}\t{depth}\t{str(dict[depth])}\n")


def dumpKurtosis(ctg, list, fh):
    kurt = str(kurtosis(list))
    fh.write(f"{ctg}\t{kurt}\n")


def initDepth():
    d = dict()
    for i in range(snakemake.params.max_depth):
        d[i] = 0
    return d


logging.debug("Parsing mpileup and printing depth on the fly")


currentdepth = initDepth()
currentcontig = str()
currentKurt = list()

outhist = open(snakemake.output.hist, 'w')
outkurt = open(snakemake.output.kurt, 'w')

with open(snakemake.input[0], 'r') as infh:
    for line in infh:
        l = line.strip().split()
        l[1] = int(l[1])
        l[2] = int(l[2])
        if not currentcontig == l[0]:
            if currentcontig:
                dumpContig(currentcontig, currentdepth, outhist)
                dumpKurtosis(currentcontig, currentKurt, outkurt)
                currentdepth = initDepth()
                currentKurt = list()
            currentcontig = l[0]
        currentKurt.append(l[1])
        if l[2] > snakemake.params.max_depth:
            l[2] = snakemake.params.max_depth
        currentdepth[l[2]] += 1


dumpContig(currentcontig, currentdepth, outhist)
dumpKurtosis(currentcontig, currentKurt, outkurt)

outhist.close()
outkurt.close()
