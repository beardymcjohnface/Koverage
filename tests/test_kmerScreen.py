import os
import tempfile
import zstandard as zstd
import numpy as np
import queue
from koverage.scripts.kmerScreen import (
    trimmed_variance,
    output_print_worker,
    process_counts,
    ref_kmer_parser_worker,
)


def test_trimmed_variance():
    data = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    expected_variance = 9.166666666666666
    actual_variance = trimmed_variance(data)
    assert np.isclose(actual_variance, expected_variance)


def test_output_print_worker():
    with tempfile.NamedTemporaryFile(delete=False) as temp_file:
        temp_file_path = temp_file.name
    out_queue = queue.Queue()
    out_queue.put("Line 1\n")
    out_queue.put("Line 2\n")
    out_queue.put(None)
    output_print_worker(out_queue=out_queue, out_file=temp_file_path)
    with open(temp_file_path, "rb") as result_file:
        compressed_data = result_file.read()
    dctx = zstd.ZstdDecompressor()
    decompressed_data = dctx.decompress(compressed_data).decode()
    lines = decompressed_data.split("\n")
    assert len(lines) == 3
    assert lines[0] == "Line 1"
    assert lines[1] == "Line 2"
    os.remove(temp_file_path)


def test_process_counts():
    kmer_counts = [1, 2, 3, 4, 5]
    expected_output = "sample\tcontig\t15\t3\t3\t1\t2.5\n"
    sample_name = "sample"
    contig_name = "contig"
    output = process_counts(kmer_counts, sample_name, contig_name)
    assert output == expected_output


def test_process_counts_with_zero_sum():
    kmer_counts = [0, 0, 0]
    sample_name = "sample"
    contig_name = "contig"
    output = process_counts(kmer_counts, sample_name, contig_name)
    assert output is None


def test_ref_kmer_parser_worker():
    with tempfile.NamedTemporaryFile(suffix=".zst", delete=False) as temp_file:
        compressor = zstd.ZstdCompressor()
        with compressor.stream_writer(temp_file) as compressed_file:
            for line in ["contig1 1 2 3 4 5", "contig2 1 2 3 4 5"]:
                compressed_file.write(line.encode() + b"\n")
    queue_out = queue.Queue()
    ref_kmer_parser_worker(
        ref_kmers=temp_file.name,
        jellyfish_db=None,
        out_queue=queue_out,
        sample_name="sample",
        cmd=["cat"],
    )
    expected_line1 = "sample\tcontig1\t15\t3\t3\t1\t2.5\n"
    expected_line2 = "sample\tcontig2\t15\t3\t3\t1\t2.5\n"
    actual_line1 = queue_out.queue[0]
    actual_line2 = queue_out.queue[1]
    assert expected_line1 == actual_line1
    assert expected_line2 == actual_line2
    assert len(queue_out.queue) == 3
    os.remove(temp_file.name)
