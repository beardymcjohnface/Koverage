#!/usr/bin/env python3

import subprocess
import threading
import queue
from statistics import variance
import os
import logging
import sys


### TODO: CAPTURE AND CHECK ERROR CODES FOR SYSTEM COMMANDS
### TODO: CAPTURE STDERR AND SAVE TO LOG FOR SYSTEM COMMANDS
### TODO: Other teadious crap that Snakemake usually handles?


# """For testing as a standalone script"""
# import attrmap as ap
# snakemake = ap.AttrMap()
# # test inputs
# snakemake.input.assembly = "ref.fa"
# snakemake.input.r1 = "test.r1.fastq.gz"
# snakemake.input.r2 = "test.r2.fastq.gz"
# snakemake.threads = 8
# snakemake.resources.mem_mb = 16000
# snakemake.params.bams = True
# # test outputs
# snakemake.params.max_depth = 300
# snakemake.params.bin_width = 50
# snakemake.output.var = "test.variance"
# snakemake.params.bamfile = "test.bam"
# snakemake.output.counts = "test.count"
# snakemake.output.lib = "test.lib"
# snakemake.log = ["test.log"]


logging.basicConfig(filename=snakemake.log[0], filemode="w", level=logging.DEBUG)


def mm_to_counts_bam(pipe, count_queue, bam_queue):
    """Capture the minimap SAM output and add it to queues for read counts, and for sorting and saving the bam"""
    for line in iter(pipe.stdout.readline, b''):
        line = line.decode()
        count_queue.put(line)
        bam_queue.put(line)
    for q in [count_queue, bam_queue]:
        q.put(None)


def mm_to_counts(pipe, count_queue):
    """Capture the minimap SAM output and add it to a queue for read counts"""
    for line in iter(pipe.stdout.readline, b''):
        line = line.decode()
        count_queue.put(line)
    count_queue.put(None)


def sort_save_bam(bam_queue):
    """Sort and save the bam file"""
    scmd = [
        "samtools",
        "sort",
        "-@", str(snakemake.threads),
        "-O", "SAM",
        "-m", str(int(snakemake.resources.mem_mb / snakemake.threads / 2)) + "M",
    ]
    vcmd = [
        "samtools",
        "view",
        "-b",
        "-h",
        "-o",
        snakemake.params.bamfile
    ]
    logging.debug(f"Starting samtools sort: {' '.join(scmd)}\n")
    pipe_sort = subprocess.Popen(scmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    logging.debug(f"Starting samtools view: {' '.join(vcmd)}\n")
    pipe_view = subprocess.Popen(vcmd, stdin=pipe_sort.stdout, stderr=subprocess.PIPE)
    while True:
        line = bam_queue.get()
        if line is None:
            break
        pipe_sort.stdin.write(line.encode())
    pipe_sort.stdin.close()
    pipe_sort.wait()
    pipe_view.wait()
    if not pipe_sort.returncode != 0:
        logging.debug(f"\nERROR: Pipe failure for:\n{' '.join(scmd)}\n")
        logging.debug(f"STDERR: {pipe_sort.stderr.read().decode()}")
        sys.exit(1)
    if not pipe_view.returncode != 0:
        logging.debug(f"\nERROR: Pipe failure for:\n{' '.join(vcmd)}\n")
        logging.debug(f"STDERR: {pipe_view.stderr.read().decode()}")
        sys.exit(1)


def count_reads(count_queue):
    """Parse the minimap SAM output -> read counts, contig lens, estimated coverage variance"""
    ctglen = dict()                         # contig lens
    ctgcnt = dict()                         # contig counts
    ctgvar = dict()                         # contig binned start coord histogram
    rcnt = 0                                # total read counts
    while True:
        line = count_queue.get()
        if line is None:
            break
        l = line.strip().split()
        if line.startswith("@"):            # SAM header
            if line.startswith("@SQ"):      # Contig part of header
                c = l[1].split(":")
                n = l[2].split(":")
                ctglen[c[1]] = n[1]
                ctgcnt[c[1]] = 0
                ctgvar[c[1]] = [0] * (int(int(n[1]) / snakemake.params.bin_width) + 1)
        else:                               # SAM body
            rcnt += 1
            if not l[2] == "*":             # * = unmapped
                ctgcnt[l[2]] += 1
                ctgvar[l[2]][int(int(l[3]) / snakemake.params.bin_width)] += 1
    with open(snakemake.output.counts, 'w') as outfh:
        for c in ctgcnt.keys():
            outfh.write(f"{c}\t{ctglen[c]}\t{ctgcnt[c]}\n")
    with open(snakemake.output.lib, 'w') as outfh:
        outfh.write(f"{str(rcnt)}\n")
    with open(snakemake.output.var, 'w') as outfh:
        for c in ctgvar.keys():
            var = variance(ctgvar[c])
            outfh.write(f"{c}\t{var}\n")


# Start minimap
mm2cmd = [
    "minimap2",
    "yeet",
    "-t",
    str(snakemake.threads),
    "-ax",
    "sr",
    "--secondary=no",
    snakemake.input.assembly,
    snakemake.input.r1,
    snakemake.input.r2
]
logging.debug(f"Starting minimap2: {' '.join(mm2cmd)}\n")
pipe_minimap = subprocess.Popen(mm2cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


# Create queue for counts
queue_counts = queue.Queue()


if snakemake.params.bams:
    queue_bam = queue.Queue()
    thread_reader = threading.Thread(target=mm_to_counts_bam, args=(pipe_minimap, queue_counts, queue_bam))
    thread_reader.start()
    thread_parser_bam = threading.Thread(target=sort_save_bam, args=(queue_bam,))
    thread_parser_bam.start()
else:
    thread_reader = threading.Thread(target=mm_to_counts, args=(pipe_minimap, queue_counts))
    thread_reader.start()
    with open(snakemake.output.bam, 'a') as b:
        os.utime(snakemake.output.bam, None)


# Read from q2 and get read counts
thread_parser_counts = threading.Thread(target=count_reads, args=(queue_counts,))
thread_parser_counts.start()


# wait for workers to finish
thread_parser_counts.join()
if snakemake.params.bams:
    thread_parser_bam.join()


# check minimap2 finished ok
pipe_minimap.stdout.close()
pipe_minimap.wait()
if pipe_minimap.returncode != 0:
    logging.debug(f"\nERROR: Pipe failure for:\n{' '.join(mm2cmd)}\n")
    logging.debug(f"STDERR: {pipe_minimap.stderr.read().decode()}")
    sys.exit(1)
