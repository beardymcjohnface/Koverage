#!/usr/bin/env python3

import logging
from scipy.stats import kurtosis


logging.basicConfig(filename=snakemake.log[0], filemode="w", level=logging.DEBUG)


logging.debug("Collecting combined coverage stats")


def dumpContig(ctg, dict, fh):
    for depth in sorted(dict.keys()):
        fh.write(f"{ctg}\t{depth}\t{dict[depth]}\n")


def dumpKurtosis(ctg, list, fh):
    kurt = kurtosis(list)
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
        if currentcontig == l[0]:
            currentKurt.append(int(l[1]))
            if l[2] > snakemake.params.max_depth:
                l[2] = snakemake.params.max_depth
            currentdepth[l[2]] += 1
        else:
            if currentcontig:
                dumpContig(currentcontig, currentdepth, outhist)
                dumpKurtosis(currentcontig, currentKurt, outkurt)
                currentdepth = initDepth()
                currentcontig = l[0]
                currentKurt = list()

dumpContig(currentcontig, currentdepth, outhist)
dumpKurtosis(currentcontig, currentKurt, outkurt)

outhist.close()
outkurt.close()
