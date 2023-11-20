import pytest
import pickle
import numpy as np
import tempfile
import os
from unittest.mock import mock_open, patch, call

from koverage.scripts.sampleCoverage import (
    calculate_coverage_stats_from_counts,
    # print_coverage_stats,
)


@pytest.fixture
def count_pickle(tmp_path):
    file_path = tmp_path / "counts.pkl"
    contig_lengths = [
        ["contig1", 1000],
        ["contig2", 5000],
        ["contig3", 10000],
    ]
    contig_bin_counts = np.array(
        [
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 1, 1, 1, 2, 0, 0, 0, 0, 0],
            [1, 1, 1, 2, 1, 0, 5, 5, 10, 5],
        ],
        dtype=np.int32,
    )
    with open(file_path, "wb") as handle:
        pickle.dump(contig_lengths, handle, protocol=pickle.HIGHEST_PROTOCOL)
        pickle.dump(contig_bin_counts, handle, protocol=pickle.HIGHEST_PROTOCOL)
    kwargs = {
        "sample": "sample1",
        "count_file": file_path,
        "bin_width": 1000,
        "output_file": "tempOutput",
    }
    return file_path, kwargs


@pytest.fixture
def count_expected_output():
    out_file_contents = [
        "sample1\tcontig1\t1\t2.632e+04\t2.632e+04\t1\t1.887e+05\t1\t1\t1\tnan\n",
        "sample1\tcontig2\t6\t1.579e+05\t3.158e+04\t1.2\t2.264e+05\t1.2\t1\t1\t0.2\n",
        "sample1\tcontig3\t31\t8.158e+05\t8.158e+04\t3.1\t5.849e+05\t3.1\t1.5\t0.9\t9.656\n"
    ]
    return out_file_contents


def test_calculate_coverage_stats_from_counts(count_pickle, count_expected_output):
    pickle_file, kwargs = count_pickle
    with tempfile.TemporaryDirectory() as temp_dir:
        kwargs["output_file"] = os.path.join(temp_dir, kwargs["output_file"])
        calculate_coverage_stats_from_counts(**kwargs)

        with open(kwargs["output_file"], 'r') as output_file:
            actual_output = output_file.readlines()

        assert actual_output == count_expected_output






# def test_print_coverage_stats(kwargs):
#     expected_output = [
#         "sample1\tcontig1\t500\t500\t500\t0.5\t0.3243\t1.5\t1.2\t0.8\t0.1\n",
#         "sample1\tcontig2\t1000\t1000\t666.7\t0.6667\t0.4324\t2.5\t2.2\t0.9\t0.2\n",
#         "sample1\tcontig3\t750\t750\t375\t0.375\t0.2432\t3.5\t3.2\t0.7\t0.3\n",
#     ]
#     with patch("builtins.open", mock_open()) as mock_file:
#         print_coverage_stats(**kwargs)
#         mock_file.assert_called_once_with("output.txt", "w")
#         handle = mock_file()
#         handle.write.assert_has_calls(
#             [call(line) for line in expected_output], any_order=False
#         )
