#!/usr/bin/env python3


"""Run minimap2, parse its output, calculate counts on the fly

This script will run minimap2 of a sample's reads against the reference FASTA.
We use a wrapper instead of a snakemake rule to avoid additional read/writes for every sample.
PAF files of alignments can optionally be saved.

- `worker_mm_to_count_paf_queues` - read minimap2 output and pass to queues for processing and saving PAF
- `worker_mm_to_count_queues` - read minimap2 output and pass to queue for processing only
- `worker_paf_writer` - read minimap2 output from queue and write to zstandard-zipped file
- `worker_count_and_print` - read minimap2 output from queue, calculate counts, print to output files
- `build_mm2cmd` - return the minimap2 command based on presence of R2 file
- `start_workers` - start queues and worker threads
"""


import subprocess
import threading
import queue
from statistics import variance
import numpy as np
import os
import logging
import sys
import zstandard as zstd


def worker_mm_to_count_paf_queues(pipe, count_queue, paf_queue):
    """Read minimap2 output and slot into queues for collecting coverage counts, and saving the paf file.

    Args:
        pipe (pipe): minimap2 pipe for reading
        count_queue (Queue): queue for putting for counts
        paf_queue (Queue): queue for putting for saving paf
    """

    for line in iter(pipe.stdout.readline, b""):
        line = line.decode()
        count_queue.put(line)
        paf_queue.put(line)

    for q in [count_queue, paf_queue]:
        q.put(None)


def worker_mm_to_count_queues(pipe, count_queue):
    """Read minimap2 output and slot into queues for collecting coverage counts

    Args:
    pipe (pipe): minimap2 pipe for reading
    count_queue (Queue): queue for putting for counts
    """

    for line in iter(pipe.stdout.readline, b""):
        line = line.decode()
        count_queue.put(line)

    count_queue.put(None)


def worker_paf_writer(paf_queue, paf_dir, sample, chunk_size=100):
    """Read minimap2 output from queue and write to zstd-zipped file

    Args:
        paf_queue (Queue): queue of minimap2 output for reading
        paf_dir (str): dir for saving paf files
    """

    cctx = zstd.ZstdCompressor()
    os.makedirs(paf_dir, exist_ok=True)
    output_f = open(os.path.join(paf_dir, sample + ".paf.zst"), "wb")
    lines = []

    while True:
        line = paf_queue.get()
        if line is None:
            break
        lines.append(line.encode())
        if len(lines) >= chunk_size:
            compressed_chunk = cctx.compress(b"".join(lines))
            output_f.write(compressed_chunk)
            lines = []

    if lines:
        compressed_chunk = cctx.compress(b"".join(lines))
        output_f.write(compressed_chunk)
        output_f.flush()

    output_f.close()


def worker_count_and_print(count_queue, **kwargs):
    """Collect the counts from minimap2 queue and calc counts on the fly

    Args:
        count_queue (Queue): queue of minimap2 output for reading
        **kwargs (dict):
            - bin_width (int): Width of bins for hitrate and variance estimation
            - output_counts (str): filepath for writing output count stats
            - output_lib (str): filepath for writing library size
    """

    contig_lengths = dict()
    contig_counts = dict()
    contig_bin_counts = dict()
    total_count = 0

    while True:
        line = count_queue.get()
        if line is None:
            break
        l = line.strip().split()

        try:
            contig_counts[l[5]] += 1
        except KeyError:
            contig_counts[l[5]] = 1
            contig_bin_counts[l[5]] = [0] * (int(int(l[6]) / kwargs["bin_width"]) + 1)
            contig_lengths[l[5]] = l[6]

        for i in range(
            int(int(l[7]) / kwargs["bin_width"]), int(int(l[6]) / kwargs["bin_width"])
        ):
            contig_bin_counts[l[5]][i] += 1
        total_count += 1

    with open(kwargs["output_counts"], "w") as out_counts:
        for c in contig_counts.keys():
            ctg_mean = "{:.{}g}".format(np.mean(contig_bin_counts[c]), 4)
            ctg_median = "{:.{}g}".format(np.median(contig_bin_counts[c]), 4)
            ctg_hitrate = "{:.{}g}".format(
                (len(contig_bin_counts[c]) - contig_bin_counts[c].count(0))
                / len(contig_bin_counts[c]),
                4,
            )
            contig_bin_counts[c] = [x / kwargs["bin_width"] for x in contig_bin_counts[c]]
            ctg_variance = "{:.{}g}".format(variance(contig_bin_counts[c]), 4)
            out_counts.write(
                "\t".join(
                    [
                        c,
                        str(contig_lengths[c]),
                        str(contig_counts[c]),
                        ctg_mean,
                        ctg_median,
                        ctg_hitrate,
                        ctg_variance + "\n",
                    ]
                )
            )

    with open(kwargs["output_lib"], "w") as out_lib:
        out_lib.write(f"{str(total_count)}\n")


def build_mm2cmd(**kwargs):
    """Return the minimap2 command

    Args:
        **kwargs (dict):
            - threads (int): Number of worker threads to use
            - minimap_mode (str): Mapping preset for minimap2
            - ref_idx (str): Reference indexed file
            - r1_file (str): Forward reads file
            - r2_file (str): Reverse reads file (or "" for SE reads/longreads)

    Returns:
        mm2cmd (list): minimap2 command for opening with subprocess
    """

    mm2cmd = [
        "minimap2",
        "-t",
        str(kwargs["threads"]),
        "-x",
        kwargs["minimap_mode"],
        "--secondary=no",
        kwargs["ref_idx"],
        kwargs["r1_file"],
    ]

    if kwargs["r2_file"] != str():
        mm2cmd.append(kwargs["r2_file"])

    return mm2cmd


def start_workers(queue_counts, queue_paf, pipe_minimap, **kwargs):
    """Start workers for reading the minimap output and parsing to queue(s) for processing

    Args:
        queue_counts (Queue): queue to use for putting minimap2 output for collecting counts
        pipe_minimap (pipe): subprocess pipe for minimap2 for reading
        **kwargs (dict):
            - paf_file (str): PAF file for writing
            - save_pafs (bool): flag for if PAF files should be saved
    """

    thread_parser_paf = None

    if kwargs["save_pafs"]:
        thread_reader = threading.Thread(
            target=worker_mm_to_count_paf_queues,
            args=(pipe_minimap, queue_counts, queue_paf),
        )
        thread_reader.start()
        thread_parser_paf = threading.Thread(
            target=worker_paf_writer,
            args=(queue_paf, kwargs["paf_dir"], kwargs["sample"]),
        )
        thread_parser_paf.start()
    else:
        thread_reader = threading.Thread(
            target=worker_mm_to_count_queues, args=(pipe_minimap, queue_counts)
        )
        thread_reader.start()

    return thread_reader, thread_parser_paf


def main(**kwargs):
    if kwargs["pyspy"]:
        subprocess.Popen(
            [
                "py-spy",
                "record",
                "-s",
                "-o",
                kwargs["pyspy_svg"],
                "--pid",
                str(os.getpid()),
            ]
        )

    logging.basicConfig(filename=kwargs["log_file"], filemode="w", level=logging.DEBUG)
    mm2cmd = build_mm2cmd(**kwargs)
    logging.debug(f"Starting minimap2: {' '.join(mm2cmd)}\n")
    pipe_minimap = subprocess.Popen(
        mm2cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )

    # Create queue for counts
    queue_counts = queue.Queue()
    queue_paf = queue.Queue()
    thread_reader, thread_parser_paf = start_workers(
        queue_counts, queue_paf, pipe_minimap, **kwargs
    )

    # Read from q2 and get read counts
    thread_parser_counts = threading.Thread(
        target=worker_count_and_print, args=(queue_counts,), kwargs=kwargs
    )
    thread_parser_counts.start()

    # wait for workers to finish
    thread_parser_counts.join()
    if thread_parser_paf:
        thread_parser_paf.join()

    # check minimap2 finished ok
    pipe_minimap.stdout.close()
    pipe_minimap.wait()
    if pipe_minimap.returncode != 0:
        logging.debug(f"\nERROR: Pipe failure for:\n{' '.join(mm2cmd)}\n")
        logging.debug(f"STDERR: {pipe_minimap.stderr.read().decode()}")
        sys.exit(1)

    # Join reader
    thread_reader.join()


if __name__ == "__main__":
    main(
        threads=snakemake.threads,
        log_file=snakemake.log.err,
        minimap_mode=snakemake.params.minimap,
        ref_idx=snakemake.input.ref,
        r1_file=snakemake.input.r1,
        r2_file=snakemake.params.r2,
        save_pafs=snakemake.params.pafs,
        paf_dir=snakemake.params.paf_dir,
        sample=snakemake.wildcards.sample,
        bin_width=snakemake.params.bin_width,
        output_counts=snakemake.output.counts,
        output_lib=snakemake.output.lib,
        pyspy=snakemake.params.pyspy,
        pyspy_svg=snakemake.log.pyspy,
    )
