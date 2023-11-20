import pytest
import pickle
import numpy as np
import tempfile
import os
from unittest.mock import mock_open, patch, call

from koverage.scripts.sampleCoverage import (
    calculate_coverage_stats_from_counts,
)


def dump_pickle(contig_lens, contig_bin_counts, **kwargs):
    with open(kwargs["count_file"], "wb") as handle:
        pickle.dump(contig_lens, handle, protocol=pickle.HIGHEST_PROTOCOL)
        pickle.dump(contig_bin_counts, handle, protocol=pickle.HIGHEST_PROTOCOL)


@pytest.fixture
def contig_lens():
    return [["contig1", 1000], ["contig2", 5000], ["contig3", 10000]]


@pytest.fixture
def contig_bin_counts():
    return np.array(
        [
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 1, 1, 1, 2, 0, 0, 0, 0, 0],
            [1, 1, 1, 2, 1, 0, 5, 5, 10, 5],
        ],
        dtype=np.int32,
    )


@pytest.fixture
def kwarguments(tmp_path):
    return {
        "sample": "sample1",
        "count_file": tmp_path / "counts.pkl",
        "bin_width": 1000,
        "output_file": tmp_path / "tempOutput",
    }


@pytest.fixture
def count_expected_output():
    out_file_contents = [
        "sample1\tcontig1\t1\t2.632e+04\t2.632e+04\t1\t1.887e+05\t1\t1\t1\t0\n",
        "sample1\tcontig2\t6\t1.579e+05\t3.158e+04\t1.2\t2.264e+05\t1.2\t1\t1\t0.2\n",
        "sample1\tcontig3\t31\t8.158e+05\t8.158e+04\t3.1\t5.849e+05\t3.1\t1.5\t0.9\t9.656\n",
    ]
    return out_file_contents


@pytest.fixture
def count_empty_output():
    out_file_contents = [
        "sample1\tcontig1\t0\t0\t0\t0\t0\t0\t0\t0\t0\n",
        "sample1\tcontig2\t0\t0\t0\t0\t0\t0\t0\t0\t0\n",
        "sample1\tcontig3\t0\t0\t0\t0\t0\t0\t0\t0\t0\n",
    ]
    return out_file_contents


def test_calculate_coverage_stats_from_counts(
    contig_lens, contig_bin_counts, kwarguments, count_expected_output
):
    dump_pickle(contig_lens, contig_bin_counts, **kwarguments)

    calculate_coverage_stats_from_counts(**kwarguments)

    with open(kwarguments["output_file"], "r") as output_file:
        actual_output = output_file.readlines()

    assert actual_output == count_expected_output


def test_empty_coverage_stats_from_counts(contig_lens, kwarguments, count_empty_output):
    contig_bin_counts = np.zeros([3, 10], dtype=np.int32)
    dump_pickle(contig_lens, contig_bin_counts, **kwarguments)

    calculate_coverage_stats_from_counts(**kwarguments)

    with open(kwarguments["output_file"], "r") as output_file:
        actual_output = output_file.readlines()

    assert actual_output == count_empty_output
