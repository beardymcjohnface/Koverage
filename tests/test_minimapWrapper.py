import pytest
import os
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
)


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


def test_worker_count_and_print(tmp_path):
    count_queue = Queue()
    output_counts = tmp_path / "counts.txt"
    output_lib = tmp_path / "lib.txt"
    input_lines = [
        "col1\tcol2\tcol3\tcol4\tcol5\tcol6\t50\t5\tcol9\tcol10\tcol11\tcol12\n",
        "col1\tcol2\tcol3\tcol4\tcol5\tcol6\t50\t5\tcol9\tcol10\tcol11\tcol12\n",
        "col1\tcol2\tcol3\tcol4\tcol5\tcol6\t50\t20\tcol9\tcol10\tcol11\tcol12\n",
    ]
    expected_counts_content = "col6\t50\t3\t2.167\t2.5\t0.8333\t0.01367\n"
    expected_lib_content = "3\n"
    for line in input_lines:
        count_queue.put(line)
    count_queue.put(None)
    worker_count_and_print(
        count_queue, output_counts=output_counts, output_lib=output_lib, bin_width=10
    )
    with open(output_counts, "r") as f:
        counts_content = f.read()
    with open(output_lib, "r") as f:
        lib_content = f.read()
    assert counts_content == expected_counts_content
    assert lib_content == expected_lib_content


def test_worker_count_and_print_empty_queue(tmp_path):
    count_queue = Queue()
    output_counts = tmp_path / "counts.txt"
    output_lib = tmp_path / "lib.txt"
    count_queue.put(None)
    worker_count_and_print(
        count_queue, output_counts=output_counts, output_lib=output_lib, bin_width=10
    )
    assert os.stat(output_counts).st_size == 0
    assert os.stat(output_lib).st_size == 2


def test_worker_count_and_print_mock_open(tmp_path):
    count_queue = Queue()
    output_counts = tmp_path / "counts.txt"
    output_lib = tmp_path / "lib.txt"
    input_lines = [
        "col1\tcol2\tcol3\tcol4\tcol5\tcol6\t50\t5\tcol9\tcol10\tcol11\tcol12\n",
        "col1\tcol2\tcol3\tcol4\tcol5\tcol6\t50\t5\tcol9\tcol10\tcol11\tcol12\n",
    ]
    for line in input_lines:
        count_queue.put(line)
    count_queue.put(None)
    with patch("builtins.open", mock_open()) as mock_file:
        worker_count_and_print(
            count_queue,
            output_counts=output_counts,
            output_lib=output_lib,
            bin_width=10,
        )
    mock_file.assert_any_call(output_counts, "w")
    mock_file.assert_any_call(output_lib, "w")


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
            call(target=worker_paf_writer, args=(paf_queue, kwargs["paf_dir"], kwargs["sample"])),
        ]
    )
    kwargs["save_pafs"] = False
    with patch("threading.Thread") as mock_thread:
        mock_thread.side_effect = [mock_thread_reader, mock_thread_parser_paf]
        start_workers(queue_counts, paf_queue, pipe_minimap, **kwargs)
    mock_thread.assert_has_calls(
        [call(target=worker_mm_to_count_queues, args=(pipe_minimap, queue_counts))]
    )
