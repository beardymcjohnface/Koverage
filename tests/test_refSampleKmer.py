import pytest
import gzip
from queue import Queue
import zstandard as zstd
from koverage.scripts.refSampleKmer import (
    parse_fasta,
    contigs_to_queue,
    string_to_kmers,
    process_contigs,
    output_printer,
)


@pytest.fixture(scope="function")
def temp_files(tmp_path):
    file_path = tmp_path / "test.fasta"
    expected_output = [">sequence1", "ACGT", ">sequence2", "TGCA"]
    file_path.write_text("\n".join(expected_output) + "\n")
    gzipped_file_path = tmp_path / "test.fasta.gz"
    with gzip.open(str(gzipped_file_path), "wt") as f:
        f.write(">sequence1\nACGT\n>sequence2\nTGCA\n")
    yield str(file_path), str(gzipped_file_path), expected_output


def test_regular_file(temp_files):
    file_path, _, expected_output = temp_files
    result = list(parse_fasta(str(file_path)))
    assert result == expected_output


def test_gzipped_file(temp_files):
    _, file_path, expected_output = temp_files
    result = list(parse_fasta(str(file_path)))
    assert result == expected_output


def contigs_to_queue_process(file_path, expected_output):
    queue = Queue()
    available_threads = 2
    expected_output = [
        {"id": expected_output[0][1:], "seq": expected_output[1]},
        {"id": expected_output[2][1:], "seq": expected_output[3]},
    ]
    contigs_to_queue(file_path, queue, available_threads, queue_hold=1)
    output = []
    while True:
        item = queue.get()
        if item is None:
            break
        output.append(item)
    assert output == expected_output


def test_contigs_to_queue_regular_file(temp_files):
    file_path, _, expected_output = temp_files
    contigs_to_queue_process(file_path, expected_output)


def test_contigs_to_queue_gzipped_file(temp_files):
    _, gzipped_file_path, expected_output = temp_files
    contigs_to_queue_process(gzipped_file_path, expected_output)


def test_string_to_kmers():
    seq = "ATCGATCGATCG"
    kwargs = {"kspace": 3, "ksize": 2, "kmin": 2, "kmax": 5}
    expected_output = {"AT", "CG"}
    result = set(string_to_kmers(seq, **kwargs))
    assert result == expected_output


def sort_output_process_contigs(string):
    elements = string.split()
    elements[1:] = sorted(elements[1:])
    sorted_string = " ".join(elements) + "\n"
    return sorted_string


def test_process_contigs():
    in_queue = Queue()
    out_queue = Queue()
    kwargs = {"kspace": 3, "ksize": 2, "kmin": 2, "kmax": 5}
    contigs = [
        {"id": "sequence1", "seq": "ATCGATCGATCG"},
        {"id": "sequence2", "seq": "ACGTACGTACGT"},
    ]
    expected_output = ["sequence1 AT CG\n", "sequence2 AC GT\n"]
    for contig in contigs:
        in_queue.put(contig)
    in_queue.put(None)
    process_contigs(in_queue, out_queue, **kwargs)
    output = []
    while True:
        item = out_queue.get()
        if item is None:
            break
        output.append(item)
    output[0] = sort_output_process_contigs(output[0])
    output[1] = sort_output_process_contigs(output[1])
    assert output == expected_output


def test_output_printer(tmp_path):
    file_path = tmp_path / "output.zst"
    expected_output = "Line 1\nLine 2\n"
    dctx = zstd.ZstdDecompressor()
    for chunk_size in [1, 10]:
        mock_queue = Queue()
        mock_queue.put("Line 1\n")
        mock_queue.put("Line 2\n")
        mock_queue.put(None)
        output_printer(mock_queue, str(file_path), chunk_size=chunk_size)
        with open(str(file_path), "rb") as in_fh:
            with dctx.stream_reader(in_fh) as f:
                decompressed_data = f.read().decode()
        assert decompressed_data == expected_output
