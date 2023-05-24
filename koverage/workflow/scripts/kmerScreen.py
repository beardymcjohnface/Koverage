#!/usr/bin/env python3

import subprocess
import io
import threading
import queue
import logging
import zstandard as zstd
from statistics import variance
import numpy as np
import sys


logging.basicConfig(filename=snakemake.log[0], filemode="w", level=logging.DEBUG)


def jellyfish_worker(in_queue, out_queue):
    cmd = [
        "jellyfish",
        "query",
        "-i",
        snakemake.input.db
    ]
    logging.debug(f"Starting interactive jellyfish session: {' '.join(cmd)}\n")
    pipe_jellyfish = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    while True:
        item = in_queue.get()
        if item is None:
            break
        counts = list()
        for k in item[1:]:
            pipe_jellyfish.stdin.write(f'{k}\n'.encode())
        pipe_jellyfish.stdin.flush()
        for _ in item[1:]:
            counts.append(int(pipe_jellyfish.stdout.readline().decode()))
        v = str(variance(counts))
        m = str(np.mean(counts))
        d = str(np.median(counts))
        out_line = ' '.join([item[0],m,d,v]) + "\n"
        out_queue.put(out_line)
    pipe_jellyfish.stdin.close()
    pipe_jellyfish.stdout.close()
    pipe_jellyfish.wait()
    if pipe_jellyfish.returncode != 0:
        logging.debug(f"\nERROR: Jellyfish failure for:\n{' '.join(cmd)}\n")
        logging.debug(f"STDERR: {pipe_jellyfish.stderr.read().decode()}")
        sys.exit(1)
    out_queue.put(None)


def output_print_worker(out_queue):
    cctx = zstd.ZstdCompressor()
    with open(snakemake.output[0], 'w') as out_fh:
        chunk_size = 100
        lines = []
        while True:
            item = out_queue.get()
            if item is None:
                break
            lines.append(item)
            if len(lines) >= chunk_size:
                compressed_chunk = cctx.compress(''.join(lines).encode())
                out_fh.write(compressed_chunk)
                lines = []
        if lines:
            compressed_chunk = cctx.compress(''.join(lines).encode())
            out_fh.write(compressed_chunk)


def ref_parser_worker(in_queue):
    with open(snakemake.input.ref, 'rb') as in_fh:
        dctx = zstd.ZstdDecompressor()
        with dctx.stream_reader(in_fh) as reader:
            wrap = io.TextIOWrapper(io.BufferedReader(reader), encoding='utf8')
            for line in wrap:
                line = line.strip()
                l = line.strip().split()
                in_queue.put(l)
    for _ in range(snakemake.threads):
        in_queue.put(None)


# open queues
queue_in = queue.Queue()
queue_out = queue.Queue()


# start jellyfish workers
jf_workers = list()
for _ in range(snakemake.threads):
    jf_thread = threading.Thread(target=jellyfish_worker, args=(queue_in, queue_out,))
    jf_thread.daemon = True
    jf_thread.start()
    jf_workers.append(jf_thread)


# start print worker
print_worker = threading.Thread(target=output_print_worker, args=(queue_out,))
print_worker.daemon = True
print_worker.start()


# start parser worker
parse_worker = threading.Thread(target=ref_parser_worker, args=(queue_in,))
parse_worker.daemon = True
parse_worker.start()


# join jellyfish workers
for t in jf_workers + [print_worker, parse_worker]:
    t.join()
