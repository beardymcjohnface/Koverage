#!/usr/bin/env python3

"""Screen the sample for reference sampled kmers

This script will parse the reference sampled kmers and query them from the sample jellyfish db

- `trimmed_variance` - Calculate variance from a list of integers
- `output_print_worker` - Take output lines from a queue and print to zstandard-compressed TSV
- `process_counts` - Process the kmer depths of the ref sampled kmers and the sample jellyfish database
"""


import os
import subprocess
import io
import threading
import queue
import logging
import zstandard as zstd
import numpy as np
import sys


def trimmed_variance(data, trim_frac=0.05):
    """Calculate the variance, minus the top x percent of outliers

    Args:
        data (list): list of int kmer depths for calculating variance
        trim_frac (float): fraction of top depths to trim (outlier handling)

    Returns:
        variance (float): variance of the trimmed kmer depths
    """
    trim_size = int(len(data) * trim_frac)
    sorted_data = np.sort(data)[::-1]
    trimmed_data = sorted_data[trim_size:]
    variance = np.var(trimmed_data, dtype=np.float64, ddof=1)
    return variance


def output_print_worker(out_queue=None, out_file=None):
    """Worker to take the output lines for printing, compress with gzip, and print to the output file.

    Args:
        out_queue (Queue): Queue with lines for printing to output file
        out_file (str): Output file for writing gzipped output
    """
    cctx = zstd.ZstdCompressor()
    with open(out_file, "wb") as out_fh:
        chunk_size = 100
        lines = []
        while True:
            item = out_queue.get()
            if item is None:
                break
            lines.append(item)
            if len(lines) >= chunk_size:
                compressed_chunk = cctx.compress("".join(lines).encode())
                out_fh.write(compressed_chunk)
                lines = []
        if lines:
            compressed_chunk = cctx.compress("".join(lines).encode())
            out_fh.write(compressed_chunk)


def process_counts(kmer_counts, sample_name, contig_name):
    """Process the kmer depths of the ref sampled kmers and the sample jellyfish database.

    Args:
        kmer_counts (list): list of kmer depths
        sample_name (str): name of the sample
        contig_name (str): contig ID

    Returns:
        out_line (str): output line for printing to output file, or None
    """
    sum_kmer = "{:.{}g}".format(np.sum(kmer_counts), 4)
    if sum_kmer != "0":
        mean_kmer = "{:.{}g}".format(np.mean(kmer_counts), 4)
        median_kmer = "{:.{}g}".format(np.median(kmer_counts), 4)
        hitrate_kmer = "{:.{}g}".format(
            (len(kmer_counts) - kmer_counts.count(0)) / len(kmer_counts), 4
        )
        variance_kmer = "{:.{}g}".format(trimmed_variance(kmer_counts), 4)
        out_line = "\t".join(
            [
                sample_name,
                contig_name,
                sum_kmer,
                mean_kmer,
                median_kmer,
                hitrate_kmer,
                variance_kmer + "\n",
            ]
        )
        return out_line
    else:
        return None


def ref_kmer_parser_worker(
    ref_kmers=None, jellyfish_db=None, out_queue=None, sample_name=None, cmd=None
):
    """Parse the processed reference kmer file (zstd-compressed) and query kmers from the Jellyfish database.

    Args:
        ref_kmers (str): Filepath of the sampled kmers from ref fasta (zstd-compressed; "contigID\tkmer\tkmer\tkmer...")
        jellyfish_db (str): Filepath of the jellyfish database for the sample
        out_queue (Queue): Queue of the lines of output for compression and writing to the output file.
        sample_name (str): Name of the sample
        cmd (list): jellyfish command. The command is passed here to allow for unit testing without invoking jellyfish.
    """
    if jellyfish_db:
        cmd.append(jellyfish_db)
    logging.debug("Starting interactive jellyfish session: " + ' '.join(cmd) + "\n")
    pipe_jellyfish = subprocess.Popen(
        cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    with open(ref_kmers, "rb") as in_fh:
        dctx = zstd.ZstdDecompressor()
        with dctx.stream_reader(in_fh) as reader:
            wrap = io.TextIOWrapper(io.BufferedReader(reader), encoding="utf8")
            for line in wrap:
                line = line.strip()
                l = line.strip().split()
                kmer_counts = list()
                for k in l[1:]:
                    k += "\n"
                    pipe_jellyfish.stdin.write(k.encode())
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
        logging.debug("\nERROR: Jellyfish failure for:\n" + ' '.join(cmd) + "\n")
        logging.debug("STDERR: " + pipe_jellyfish.stderr.read().decode())
        sys.exit(1)
    out_queue.put(None)


def main(**kwargs):
    # if kwargs["pyspy"]:
    #     subprocess.Popen(
    #         [
    #             "py-spy",
    #             "record",
    #             "-s",
    #             "-o",
    #             kwargs["pyspy_svg"],
    #             "--pid",
    #             str(os.getpid()),
    #         ]
    #     )
    logging.basicConfig(filename=kwargs["log_file"], filemode="w", level=logging.DEBUG)
    # open printing queue
    queue_out = queue.Queue()
    # start print worker
    print_worker = threading.Thread(
        target=output_print_worker,
        kwargs={"out_queue": queue_out, "out_file": kwargs["out_file"]},
    )
    print_worker.daemon = True
    print_worker.start()
    # start parser worker
    parse_worker = threading.Thread(
        target=ref_kmer_parser_worker,
        kwargs={
            "ref_kmers": kwargs["ref_kmers"],
            "jellyfish_db": kwargs["jellyfish_db"],
            "out_queue": queue_out,
            "sample_name": kwargs["sample_name"],
            "cmd": ["jellyfish", "query", "-i"],
        },
    )
    parse_worker.daemon = True
    parse_worker.start()
    # join jellyfish workers
    for t in [print_worker, parse_worker]:
        t.join()


if __name__ == "__main__":
    main(
        jellyfish_db=snakemake.input.db,
        log_file=snakemake.log[0],
        ref_kmers=snakemake.input.ref,
        sample_name=snakemake.wildcards.sample,
        out_file=snakemake.output[0],
        # pyspy=snakemake.params.pyspy,
        # pyspy_svg=snakemake.log.pyspy,
    )
