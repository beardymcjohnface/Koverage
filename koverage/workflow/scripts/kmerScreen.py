#!/usr/bin/env python3

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
            counts.append(int(pipe_jellyfish.stdout.readline().decode()))
        v = variance(counts)
        m = np.mean(counts)
        d = np.median(counts)
        out_queue.put(' '.join([item[0],m,d,v]) + "\n")
    pipe_jellyfish.stdin.close()
    pipe_jellyfish.stdout.close()
    pipe_jellyfish.wait()
    if pipe_jellyfish.returncode != 0:
        logging.debug(f"\nERROR: Jellyfish failure for:\n{' '.join(cmd)}\n")
        logging.debug(f"STDERR: {pipe_jellyfish.stderr.read().decode()}")
        sys.exit(1)


def output_print_worker(out_queue):
    with open(snakemake.output[0], 'r') as out_fh:
        while True:
            line = out_queue.get()
            if line is None:
                break
            out_fh.write(line)


def ref_parser_worker(in_queue):
    with open(snakemake.input.ref, 'rb') as in_fh:
        dctx = zstd.ZstdDecompressor()
        with dctx.stream_reader(in_fh) as reader:
            for line in reader:
                line = line.decode().strip()
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


# start parser worker
parse_worker = threading.Thread(target=ref_parser_worker, args=(queue_in,))


# join jellyfish workers
for t in jf_workers + [print_worker, parse_worker]:
    t.join()
