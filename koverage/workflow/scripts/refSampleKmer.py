#!/usr/bin/env python3

import threading
import queue
import gzip
import logging
import random


logging.basicConfig(filename=snakemake.log[0], filemode="w", level=logging.DEBUG)


ksize = snakemake.params.ksize
kspace = snakemake.params.kspace
kmin = snakemake.params.kmin
kmax = snakemake.params.kmax


def parse_fasta(file):
    if file.endswith(".gz"):
        with gzip.open(file, 'rt') as f:
            for line in f:
                yield line.strip()
    else:
        with open(file, 'r') as f:
            for line in f:
                yield line.strip()


def contigs_to_queue(file, queue):
    id = str()
    seq = str()
    for line in parse_fasta(file):
        if line.startswith(">"):
            if seq:
                queue.put([id,seq])
            seq = line
            id = str()
        else:
            seq += line
    if seq:
        queue.put([id,seq])
    for _ in range(snakemake.threads):
        queue.put(None)


def string_to_kmers(seq):
    nkmer = len(seq) / kspace
    if nkmer < kmin:
        nkmer = kmin
    elif nkmer > kmax:
        nkmer = kmax
    kmers = set()
    for i in range(len(seq) - ksize + 1):
        kmers.add(seq[i:i + ksize])
    kmers = list(kmers)
    random.shuffle(kmers)
    return kmers[0:nkmer]


def process_contigs(in_queue, out_queue):
    while True:
        item = in_queue.get()      # item[0] = id, item[1] = seq
        if item is None:
            break
        out_queue.put(' '.join([item[0]] + string_to_kmers(item[1])))
    out_queue.put(None)


def output_printer(queue):
    with open(snakemake.output[0], 'w') as outfh:
        while True:
            item = queue.get()
            if item is None:
                break
            outfh.write(item + "\n")


# create queue
contig_queue = queue.Queue()
output_queue = queue.Queue()


# create fasta parser
thread_parser_counts = threading.Thread(target=contigs_to_queue, args=(snakemake.input[0], contig_queue,))
thread_parser_counts.start()


# spin up some workers
threads = list()
for _ in range(snakemake.threads):
    thread_worker = threading.Thread(target=process_contigs, args=(contig_queue,output_queue,))
    thread_worker.start()
    threads.append(thread_worker)


# create the output writer worker
thread_printer_output = threading.Thread(target=output_printer, args=(output_queue,))
thread_printer_output.start()


# wait for everything to finish
for t in [thread_parser_counts, thread_printer_output] + threads:
    t.join()

