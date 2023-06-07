#!/usr/bin/env python3

import threading
import queue
import gzip
import logging
import zstandard as zstd
import time


def parse_fasta(file):
    if file.endswith(".gz"):
        with gzip.open(file, 'rt') as f:
            for line in f:
                yield line.strip()
    else:
        with open(file, 'r') as f:
            for line in f:
                yield line.strip()


def contigs_to_queue(file, queue_put, available_threads, queue_hold=1000):
    id = str()
    seq = str()
    for line in parse_fasta(file):
        if line.startswith(">"):
            if seq:
                if queue_put.qsize() > queue_hold:
                    time.sleep(1)
                queue_put.put({"id":id,"seq":seq})
            l = line.strip().split()
            id = l[0].replace('>','')
            seq = str()
        else:
            seq += line.strip()
    if seq:
        queue_put.put({"id":id,"seq":seq})
    for _ in range(available_threads):
        queue_put.put(None)


def string_to_kmers(seq, **kwargs):
    nkmer = int(len(seq) / kwargs["kspace"])
    imax = len(seq) - (kwargs["ksize"] + 1)
    if nkmer < kwargs["kmin"]:
        nkmer = kwargs["kmin"]
    elif nkmer > kwargs["kmax"]:
        nkmer = kwargs["kmax"]
    kpad = int(imax / nkmer)
    if kpad < 1:
        kpad = 1
    kmers = set()
    for i in range(nkmer):
        start = i * kpad
        if start < imax:
            kmers.add(seq[start:start + kwargs["ksize"]])
    kmers = list(kmers)
    return kmers


def process_contigs(in_queue, out_queue, **kwargs):
    while True:
        item = in_queue.get()
        if item is None:
            break
        outKmer = ' '.join(string_to_kmers(item["seq"], **kwargs))
        out_queue.put(f"{item['id']} {outKmer}\n")
    out_queue.put(None)


def output_printer(queue, outfile, chunk_size=100):
    cctx = zstd.ZstdCompressor()
    with open(outfile, 'wb') as outfh:
        lines = []
        while True:
            item = queue.get()
            if item is None:
                break
            lines.append(item)
            if len(lines) >= chunk_size:
                compressed_chunk = cctx.compress(''.join(lines).encode())
                outfh.write(compressed_chunk)
                lines = []
        if lines:
            compressed_chunk = cctx.compress(''.join(lines).encode())
            outfh.write(compressed_chunk)


def main(**kwargs):
    logging.basicConfig(filename=kwargs["log_file"], filemode="w", level=logging.DEBUG)
    # create queue
    contig_queue = queue.Queue()
    output_queue = queue.Queue()
    # create fasta parser
    thread_parser_counts = threading.Thread(target=contigs_to_queue, args=(kwargs["input_file"], contig_queue, kwargs["threads"]))
    thread_parser_counts.daemon = True
    thread_parser_counts.start()
    # spin up some workers
    threads = list()
    for _ in range(kwargs["threads"]):
        thread_worker = threading.Thread(target=process_contigs, args=(contig_queue,output_queue,), kwargs=kwargs)
        thread_worker.daemon = True
        thread_worker.start()
        threads.append(thread_worker)
    # create the output writer worker
    thread_printer_output = threading.Thread(target=output_printer, args=(output_queue,kwargs["output_file"]))
    thread_printer_output.daemon = True
    thread_printer_output.start()
    # wait for everything to finish
    for t in [thread_parser_counts, thread_printer_output] + threads:
        t.join()


if __name__ == "__main__":
    main(log_file = snakemake.log[0],
         input_file=snakemake.input[0],
         output_file=snakemake.output[0],
         threads=snakemake.threads,
         ksize = snakemake.params.ksize,
         kspace = snakemake.params.kspace,
         kmin = snakemake.params.kmin,
         kmax = snakemake.params.kmax)
