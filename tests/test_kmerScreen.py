import os
import pytest
import tempfile
import zstandard as zstd
import numpy as np
from queue import Queue
from koverage.workflow.scripts.kmerScreen import (
    trimmed_variance, output_print_worker, process_counts, ref_parser_worker)


def test_trimmed_variance():
    data = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    expected_variance = 6.666666666666667
    actual_variance = trimmed_variance(data)
    assert np.isclose(actual_variance, expected_variance)


def test_output_print_worker():
    with tempfile.NamedTemporaryFile(delete=False) as temp_file:
        temp_file_path = temp_file.name
    queue = Queue()
    queue.put("Line 1")
    queue.put("Line 2")
    queue.put(None)
    output_print_worker(out_queue=queue, out_file=temp_file_path)
    with open(temp_file_path, 'rb') as result_file:
        compressed_data = result_file.read()
    dctx = zstd.ZstdDecompressor()
    decompressed_data = dctx.decompress(compressed_data).decode()
    lines = decompressed_data.split('\n')
    assert len(lines) == 3
    assert lines[0] == "Line 1"
    assert lines[1] == "Line 2"
    os.remove(temp_file_path)


def test_process_counts():
    kmer_counts = [1, 2, 3, 4, 5]
    expected_output = 'sample\tcontig\t15\t3\t3\t1\t2.5\n'
    sample_name = 'sample'
    contig_name = 'contig'
    output = process_counts(kmer_counts, sample_name, contig_name)
    assert output == expected_output


def test_process_counts_with_zero_sum():
    kmer_counts = [0, 0, 0]
    output = process_counts(kmer_counts)
    assert output is None


# still need a test for ref_parser_worker