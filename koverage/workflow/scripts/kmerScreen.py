#!/usr/bin/env python3

import subprocess
import io
import threading
import queue
import logging
import zstandard as zstd
import numpy as np
import sys


def trimmed_variance(data):
    trim_size = int(len(data) * 0.05)
    sorted_data = np.sort(data)[::-1]
    trimmed_data = sorted_data[trim_size:]
    variance = np.var(trimmed_data, dtype=np.float64, ddof=1)
    return variance


def output_print_worker(out_queue=None, out_file=None):
    cctx = zstd.ZstdCompressor()
    with open(out_file, 'wb') as out_fh:
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


def process_counts(kmer_counts, sample_name, contig_name):
    sum_kmer = "{:.{}g}".format(np.sum(kmer_counts), 4)
    if sum_kmer != "0":
        mean_kmer = "{:.{}g}".format(np.mean(kmer_counts), 4)
        median_kmer = "{:.{}g}".format(np.median(kmer_counts), 4)
        hitrate_kmer = "{:.{}g}".format((len(kmer_counts) - kmer_counts.count(0)) / len(kmer_counts), 4)
        variance_kmer = "{:.{}g}".format(trimmed_variance(kmer_counts), 4)
        out_line = '\t'.join([
            sample_name,
            contig_name,
            sum_kmer,
            mean_kmer,
            median_kmer,
            hitrate_kmer,
            variance_kmer + "\n"
        ])
        return out_line
    else:
        return None


def ref_parser_worker(
        ref_kmers=None,
        jellyfish_db=None,
        out_queue=None,
        sample_name=None,
        cmd=["jellyfish", "query", "-i"]):
    if jellyfish_db:
        cmd.append(jellyfish_db)
    logging.debug(f"Starting interactive jellyfish session: {' '.join(cmd)}\n")
    pipe_jellyfish = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    with open(ref_kmers, 'rb') as in_fh:
        dctx = zstd.ZstdDecompressor()
        with dctx.stream_reader(in_fh) as reader:
            wrap = io.TextIOWrapper(io.BufferedReader(reader), encoding='utf8')
            for line in wrap:
                line = line.strip()
                l = line.strip().split()
                kmer_counts = list()
                for k in l[1:]:
                    pipe_jellyfish.stdin.write(f'{k}\n'.encode())
                pipe_jellyfish.stdin.flush()
                for _ in l[1:]:
                    kmer_counts.append(int(pipe_jellyfish.stdout.readline().decode()))
                    out_line = process_counts(kmer_counts, sample_name, l[0])
                    if out_line:
                        out_queue.put(out_line)
    pipe_jellyfish.stdin.close()
    pipe_jellyfish.stdout.close()
    pipe_jellyfish.wait()
    if pipe_jellyfish.returncode != 0:
        logging.debug(f"\nERROR: Jellyfish failure for:\n{' '.join(cmd)}\n")
        logging.debug(f"STDERR: {pipe_jellyfish.stderr.read().decode()}")
        sys.exit(1)
    out_queue.put(None)


def main(**kwargs):
    logging.basicConfig(filename=log_file, filemode="w", level=logging.DEBUG)
    # open printing queue
    queue_out = queue.Queue()
    # start print worker
    print_worker = threading.Thread(
        target=output_print_worker,
        kwargs={"out_queue":queue_out, "out_file": kwargs["out_file"]})
    print_worker.daemon = True
    print_worker.start()
    # start parser worker
    parse_worker = threading.Thread(
        target=ref_parser_worker,
        kwargs={"ref_kmers":kwargs["ref_kmers"], "jellyfish_db":kwargs["jellyfish_db"], "out_queue":queue_out})
    parse_worker.daemon = True
    parse_worker.start()
    # join jellyfish workers
    for t in [print_worker, parse_worker]:
        t.join()

if __name__ == "__main__":
    main(jellyfish_db=snakemake.input.db,
         log_file=snakemake.log[0],
         ref_kmers=snakemake.input.ref,
         sample_name=snakemake.wildcards.sample,
         out_file=snakemake.output[0])