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

logging.debug("Collecting combined coverage stats")


def dumpContig(ctg, dict, fh):
    for depth in dict.keys().sorted():
        fh.write(f"{ctg}\t{depth}\t{dict[depth]}\n")


def initDepth():
    d = dict()
    for i in range(snakemake.params.maxDepth):
        d[i] = 0
    return d


logging.debug("Parsing mpileup and printing depth on the fly")

currentdepth = initDepth()
currentcontig = str()

with open(snakemake.input[0], 'r') as infh:
    for line in infh:
        l = line.strip().split()
        if currentcontig == l[0]:
            if l[1] > snakemake.params.maxDepth:
                l[1] = snakemake.params.maxDepth
            currentdepth[l[1]] += 1
        else:
            if currentcontig:
                dumpContig(currentcontig, currentdepth, infh)
                currentdepth = initDepth()
                currentcontig = l[0]

dumpContig(currentcontig, currentdepth, infh)
