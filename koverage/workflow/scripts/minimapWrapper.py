#!/usr/bin/env python3

import subprocess
import threading
import queue
from statistics import variance
import os
import logging
import sys
import zstandard as zstd


### TODO: Any missing teadious crap that Snakemake usually handles?


logging.basicConfig(filename=snakemake.log[0], filemode="w", level=logging.DEBUG)


def worker_mm_to_count_paf_queues(pipe, count_queue, paf_queue):
    for line in iter(pipe.stdout.readline, b""):
        line = line.decode()
        count_queue.put(line)
        paf_queue.put(line)
    for q in [count_queue, paf_queue]:
        q.put(None)


def worker_mm_to_count_queues(pipe, count_queue):
    """Capture the minimap paf output and add it to a queue for read counts"""
    for line in iter(pipe.stdout.readline, b""):
        line = line.decode()
        count_queue.put(line)
    count_queue.put(None)


def worker_paf_writer(paf_queue):
    cctx = zstd.ZstdCompressor()
    with open(snakemake.output.paf, "wb") as output_f:
        chunk_size = 100
        lines = []
        while True:
            line = paf_queue.get()
            if line is None:
                break
            lines.append(line)
            if len(lines) >= chunk_size:
                compressed_chunk = cctx.compress("".join(lines).encode())
                output_f.write(compressed_chunk)
                lines = []
        if lines:
            compressed_chunk = cctx.compress("".join(lines).encode())
            output_f.write(compressed_chunk)


def worker_count_and_print(count_queue):
    """Parse the minimap SAM output -> read counts, contig lens, estimated coverage variance"""
    ctglen = dict()                         # contig lens
    ctgcnt = dict()                         # contig counts
    ctgvar = dict()                         # contig binned start coord histogram
    rcnt = 0                                # total mapped read counts
    while True:
        line = count_queue.get()
        if line is None:
            break
        l = line.strip().split()
        try:
            ctgcnt[l[5]] += 1
        except KeyError:
            ctgcnt[l[5]] = 1
            ctgvar[l[5]] = [0] * (int(int(l[6]) / snakemake.params.bin_width) + 1)
            ctglen[l[5]] = l[6]
        ctgvar[l[5]][int(int(l[7]) / snakemake.params.bin_width)] += 1
        rcnt += 1
    with open(snakemake.output.counts, "w") as outfh:
        for c in ctgcnt.keys():
            outfh.write(f"{c}\t{ctglen[c]}\t{ctgcnt[c]}\n")
    with open(snakemake.output.lib, "w") as outfh:
        outfh.write(f"{str(rcnt)}\n")
    with open(snakemake.output.var, "w") as outfh:
        for c in ctgvar.keys():
            hitrate = "{:.{}g}".format((len(ctgvar[c]) - ctgvar[c].count(0)) / len(ctgvar[c]), 4)
            var = "{:.{}g}".format(variance(ctgvar[c]), 4)
            outfh.write(f"{c}\t{hitrate}\t{var}\n")


# Start minimap
mm2cmd = [
    "minimap2",
    "-t",
    str(snakemake.threads),
    "-x",
    snakemake.params.minimap,
    "--secondary=no",
    snakemake.input.ref,
    snakemake.input.r1,
]
if snakemake.params.r2 != str():
    mm2cmd.append(snakemake.params.r2)

logging.debug(f"Starting minimap2: {' '.join(mm2cmd)}\n")
pipe_minimap = subprocess.Popen(mm2cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


# Create queue for counts
queue_counts = queue.Queue()


if snakemake.params.pafs:
    queue_paf = queue.Queue()
    thread_reader = threading.Thread(target=worker_mm_to_count_paf_queues, args=(pipe_minimap, queue_counts, queue_paf))
    thread_reader.start()
    thread_parser_paf = threading.Thread(target=worker_paf_writer, args=(queue_paf,))
    thread_parser_paf.start()
else:
    thread_reader = threading.Thread(target=worker_mm_to_count_queues, args=(pipe_minimap, queue_counts))
    thread_reader.start()
    with open(snakemake.output.paf, "a") as b:
        os.utime(snakemake.output.paf, None)


# Read from q2 and get read counts
thread_parser_counts = threading.Thread(target=worker_count_and_print, args=(queue_counts,))
thread_parser_counts.start()


# wait for workers to finish
thread_parser_counts.join()
if snakemake.params.pafs:
    thread_parser_paf.join()


# check minimap2 finished ok
pipe_minimap.stdout.close()
pipe_minimap.wait()
if pipe_minimap.returncode != 0:
    logging.debug(f"\nERROR: Pipe failure for:\n{' '.join(mm2cmd)}\n")
    logging.debug(f"STDERR: {pipe_minimap.stderr.read().decode()}")
    sys.exit(1)
