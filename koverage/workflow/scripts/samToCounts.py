#!/usr/bin/env python3

import os
import atexit
import logging
import re


def exitLogCleanup(*args):
    """Cleanup the logging file(s) prior to exiting"""
    for logFile in args:
        os.unlink(logFile)
    return None


atexit.register(exitLogCleanup, snakemake.log[0])
logging.basicConfig(filename=snakemake.log[0], filemode="w", level=logging.DEBUG)

logging.debug("")

sam = open(snakemake.input[0], 'r')
if snakemake.params:
    outSam = open(snakemake.output.sam, 'a')

contigLens = dict()
contigCounts = dict()

# parse header
for line in sam:
    if snakemake.params:
        outSam.write(line) # shouldn't need to catch nameerror exceptions
    if line.startswith("@"):
        if line.startswith("@SQ"):
            l = line.strip().split()
            c = l[1].split(":")
            n = l[2].split(":")
            contigLens[c[1]] = n[1]
            contigCounts[c[1]] = 0
    else:
        l = line.strip().split()
        contigCounts[l[2]] = contigCounts[l[2]] + 1
        break

# parse body and echo to output pipe
if snakemake.params:
    for line in sam:
        outSam.write(line) # shouldn't need to catch nameerror exceptions
        l = line.strip().split()
        contigCounts[l[2]] = contigCounts[l[2]] + 1

# parse body, don't echo to output pipe - yes it's redundant coding but it should be faster
else:
    for line in sam:
        l = line.strip().split()
        contigCounts[l[2]] = contigCounts[l[2]] + 1

# print output
with open(snakemake.output.tsv, 'w') as outTsv:
    for contig in contigCounts.keys():
        outTsv.write(f"{contig}\t{contigLens[contig]}\t{contigCounts[contig]}\n")
