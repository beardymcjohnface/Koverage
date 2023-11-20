import pytest
import tempfile
import os
import pickle
import numpy as np
from numpy.testing import assert_array_equal
from queue import Queue
import zstandard as zstd
from unittest.mock import MagicMock, patch, mock_open, call

from koverage.scripts.minimapWrapper import (
    worker_mm_to_count_paf_queues,
    worker_mm_to_count_queues,
    worker_paf_writer,
    worker_count_and_print,
    build_mm2cmd,
    start_workers,
    contig_lens_from_fai,
)


@pytest.fixture
def fasta_content():
    return "seq1 50 0 11 11\n" "seq2 125 50 22 22\n" "seq3 100 150 33 33\n"


@pytest.fixture
def fasta_lens():
    return [("seq1", 50), ("seq2", 125), ("seq3", 100)]


@pytest.fixture
def minimap_pipe():
    out = [
        "col1\tcol2\tcol3\tcol4\tcol5\t0\t50\t25\tcol9\tcol10\tcol11\tcol12\n",
        "col1\tcol2\tcol3\tcol4\tcol5\t0\t50\t25\tcol9\tcol10\tcol11\tcol12\n",
        "col1\tcol2\tcol3\tcol4\tcol5\t0\t50\t25\tcol9\tcol10\tcol11\tcol12\n",
        "col1\tcol2\tcol3\tcol4\tcol5\t1\t125\t25\tcol9\tcol10\tcol11\tcol12\n",
        "col1\tcol2\tcol3\tcol4\tcol5\t1\t125\t125\tcol9\tcol10\tcol11\tcol12\n",
        "col1\tcol2\tcol3\tcol4\tcol5\t2\t100\t25\tcol9\tcol10\tcol11\tcol12\n",
        "col1\tcol2\tcol3\tcol4\tcol5\t2\t100\t25\tcol9\tcol10\tcol11\tcol12\n",
        "col1\tcol2\tcol3\tcol4\tcol5\t2\t100\t75\tcol9\tcol10\tcol11\tcol12\n",
    ]
    return out


@pytest.fixture
def numpy_count_arr():
    out = np.array([[3, 0, 0], [1, 0, 1], [2, 1, 0]], dtype=np.int32)
    return out


def test_worker_mm_to_count_paf_queues():
    count_queue = Queue()
    paf_queue = Queue()
    pipe = MagicMock()
    pipe.stdout.readline.side_effect = [
        b"line1\n",
        b"line2\n",
        b"line3\n",
        b"",
    ]
    worker_mm_to_count_paf_queues(pipe, count_queue, paf_queue)
    assert count_queue.get() == "line1\n"
    assert paf_queue.get() == "line1\n"
    assert count_queue.get() == "line2\n"
    assert paf_queue.get() == "line2\n"
    assert count_queue.get() == "line3\n"
    assert paf_queue.get() == "line3\n"
    assert count_queue.get() is None
    assert paf_queue.get() is None


def test_worker_mm_to_count_queues():
    count_queue = Queue()
    pipe = MagicMock()
    pipe.stdout.readline.side_effect = [
        b"line1\n",
        b"line2\n",
        b"line3\n",
        b"",
    ]
    worker_mm_to_count_queues(pipe, count_queue)
    assert count_queue.get() == "line1\n"
    assert count_queue.get() == "line2\n"
    assert count_queue.get() == "line3\n"
    assert count_queue.get() is None


def test_worker_paf_writer(tmp_path):
    paf_queue = Queue()
    paf_dir = tmp_path / "output"
    sample = "sample"
    paf_file = os.path.join(paf_dir, sample + ".paf.zst")
    paf_queue.put("line1\n")
    paf_queue.put("line2\n")
    paf_queue.put("line3\n")
    paf_queue.put(None)
    worker_paf_writer(paf_queue, paf_dir, sample)
    with open(paf_file, "rb") as f:
        compressed_content = f.read()
    assert compressed_content.startswith(b"\x28\xb5\x2f\xfd")
    assert b"line1\nline2\nline3\n" in compressed_content


def test_worker_paf_writer_empty_queue(tmp_path):
    paf_queue = Queue()
    paf_dir = tmp_path / "output"
    sample = "sample"
    paf_file = os.path.join(paf_dir, sample + ".paf.zst")
    paf_queue.put(None)
    worker_paf_writer(paf_queue, paf_dir, sample)
    assert os.stat(paf_file).st_size == 0


def test_worker_paf_writer_chunksize(tmp_path):
    paf_queue = Queue()
    paf_dir = tmp_path / "output"
    sample = "sample"
    paf_file = os.path.join(paf_dir, sample + ".paf.zst")
    paf_queue.put("line1\n")
    paf_queue.put("line2\n")
    paf_queue.put("line3\n")
    paf_queue.put(None)
    worker_paf_writer(paf_queue, paf_dir, sample, chunk_size=2)
    dctx = zstd.ZstdDecompressor()
    with open(paf_file, "rb") as in_fh:
        with dctx.stream_reader(in_fh) as f:
            decompressed_data = f.read()
    decoded_output = decompressed_data.decode()
    expected_output = "line1\nline2\nline3\n"
    assert decoded_output == expected_output


def test_worker_count_and_print(tmp_path, fasta_lens, minimap_pipe, numpy_count_arr):
    count_queue = Queue()
    output_counts = tmp_path / "counts.pkl"

    for line in minimap_pipe:
        count_queue.put(line)
    count_queue.put(None)

    worker_count_and_print(
        count_queue, fasta_lens, output_counts=output_counts, bin_width=50
    )

    with open(output_counts, "rb") as handle:
        contig_lengths = pickle.load(handle)
        contig_bin_counts = pickle.load(handle)

    assert contig_lengths == fasta_lens
    assert_array_equal(contig_bin_counts, numpy_count_arr)


def test_worker_count_and_print_empty_queue(tmp_path, fasta_lens):
    count_queue = Queue()
    output_counts = tmp_path / "counts.pkl"
    count_queue.put(None)
    numpy_empty_counts = np.zeros([3, 3], dtype=np.int32)

    worker_count_and_print(
        count_queue, fasta_lens, output_counts=output_counts, bin_width=50
    )

    with open(output_counts, "rb") as handle:
        contig_lengths = pickle.load(handle)
        contig_bin_counts = pickle.load(handle)

    assert contig_lengths == fasta_lens
    assert_array_equal(contig_bin_counts, numpy_empty_counts)


def test_worker_count_and_print_mock_open(tmp_path, minimap_pipe, fasta_lens):
    count_queue = Queue()
    output_counts = tmp_path / "counts.pkl"

    for line in minimap_pipe:
        count_queue.put(line)
    count_queue.put(None)

    with patch("builtins.open", mock_open()) as mock_file:
        worker_count_and_print(
            count_queue,
            fasta_lens,
            output_counts=output_counts,
            bin_width=50,
        )

    mock_file.assert_any_call(output_counts, "wb")


def test_build_mm2cmd():
    kwargs = {
        "threads": 4,
        "minimap_mode": "sr",
        "ref_idx": "ref.idx",
        "r1_file": "reads1.fastq",
        "r2_file": "reads2.fastq",
    }
    expected_cmd = [
        "minimap2",
        "-t",
        "4",
        "-x",
        "sr",
        "--secondary=no",
        "ref.idx",
        "reads1.fastq",
        "reads2.fastq",
    ]
    assert build_mm2cmd(**kwargs) == expected_cmd
    kwargs["r2_file"] = ""
    expected_cmd = [
        "minimap2",
        "-t",
        "4",
        "-x",
        "sr",
        "--secondary=no",
        "ref.idx",
        "reads1.fastq",
    ]
    assert build_mm2cmd(**kwargs) == expected_cmd
    kwargs = {
        "threads": 8,
        "minimap_mode": "map-ont",
        "ref_idx": "ref.fa",
        "r1_file": "reads.fq",
        "r2_file": "",
    }
    expected_cmd = [
        "minimap2",
        "-t",
        "8",
        "-x",
        "map-ont",
        "--secondary=no",
        "ref.fa",
        "reads.fq",
    ]
    assert build_mm2cmd(**kwargs) == expected_cmd


def test_start_workers_mock_thread():
    kwargs = {"save_pafs": True, "paf_dir": "output", "sample": "sample"}
    queue_counts = Queue()
    paf_queue = Queue()
    pipe_minimap = MagicMock()
    mock_thread_reader = MagicMock()
    mock_thread_parser_paf = MagicMock()
    with patch("threading.Thread") as mock_thread:
        mock_thread.side_effect = [mock_thread_reader, mock_thread_parser_paf]
        start_workers(queue_counts, paf_queue, pipe_minimap, **kwargs)
    mock_thread.assert_has_calls(
        [
            call(
                target=worker_mm_to_count_paf_queues,
                args=(pipe_minimap, queue_counts, paf_queue),
            ),
            call(
                target=worker_paf_writer,
                args=(paf_queue, kwargs["paf_dir"], kwargs["sample"]),
            ),
        ]
    )
    kwargs["save_pafs"] = False
    with patch("threading.Thread") as mock_thread:
        mock_thread.side_effect = [mock_thread_reader, mock_thread_parser_paf]
        start_workers(queue_counts, paf_queue, pipe_minimap, **kwargs)
    mock_thread.assert_has_calls(
        [call(target=worker_mm_to_count_queues, args=(pipe_minimap, queue_counts))]
    )


def test_contig_lens_from_fai(fasta_content, fasta_lens):
    with tempfile.NamedTemporaryFile(delete=False) as temp_file:
        temp_file.write(fasta_content.encode())
        temp_file.seek(0)
        result = contig_lens_from_fai(temp_file.name)
        assert result == fasta_lens


def test_contig_lens_from_fai_empty():
    with tempfile.NamedTemporaryFile(delete=False) as temp_file:
        temp_file.write("".encode())
        temp_file.seek(0)
        result = contig_lens_from_fai(temp_file.name)
        assert result == []


def test_contig_lens_from_fai_invalid():
    with tempfile.NamedTemporaryFile(delete=False) as temp_file:
        temp_file.write("ACGTACGT\n".encode())
        temp_file.seek(0)
        result = contig_lens_from_fai(temp_file.name)
        assert result == []
