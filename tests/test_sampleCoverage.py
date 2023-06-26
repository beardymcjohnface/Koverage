import pytest
from unittest.mock import mock_open, patch, call

from koverage.scripts.sampleCoverage import (
    calculate_coverage_stats_from_counts,
    print_coverage_stats,
)


@pytest.fixture
def lib_file(tmp_path):
    file_path = tmp_path / "library_size.txt"
    file_path.write_text("1000000\n")
    return file_path


@pytest.fixture
def count_file(tmp_path):
    file_path = tmp_path / "count_data.txt"
    file_path.write_text(
        (
            "contig1\t1000\t500\t1.5\t1.2\t0.5\t0.2\n"
            "contig2\t1500\t1000\t2.5\t2.2\t0.8\t0.1\n"
            "contig3\t2000\t750\t3.5\t3.2\t0.6\t0.3\n"
        )
    )
    return file_path


def test_calculate_coverage_stats_from_counts(lib_file, count_file):
    expected_counts = {
        "contig1": {
            "count": "500",
            "rpm": 500.0,
            "rpkm": 500.0,
            "rpk": 500.0,
            "mean": "1.5",
            "median": "1.2",
            "hitrate": "0.5",
            "variance": "0.2",
        },
        "contig2": {
            "count": "1000",
            "rpm": 1000.0,
            "rpkm": 666.6666666666666,
            "rpk": 666.66666666666666,
            "mean": "2.5",
            "median": "2.2",
            "hitrate": "0.8",
            "variance": "0.1",
        },
        "contig3": {
            "count": "750",
            "rpm": 750.0,
            "rpkm": 375.0,
            "rpk": 375.0,
            "mean": "3.5",
            "median": "3.2",
            "hitrate": "0.6",
            "variance": "0.3",
        },
    }
    expected_rpkscale = 0.0015416666666666667
    counts, rpkscale = calculate_coverage_stats_from_counts(lib_file, count_file)
    assert counts == expected_counts
    assert round(rpkscale, 5) == round(expected_rpkscale, 5)


@pytest.fixture
def kwargs():
    return {
        "sample": "sample1",
        "counts": {
            "contig1": {
                "count": "500",
                "rpm": 500.0,
                "rpkm": 500.0,
                "rpk": 0.5,
                "mean": "1.5",
                "median": "1.2",
                "hitrate": "0.8",
                "variance": "0.1",
            },
            "contig2": {
                "count": "1000",
                "rpm": 1000.0,
                "rpkm": 666.6666666666666,
                "rpk": 0.6666666666666666,
                "mean": "2.5",
                "median": "2.2",
                "hitrate": "0.9",
                "variance": "0.2",
            },
            "contig3": {
                "count": "750",
                "rpm": 750.0,
                "rpkm": 375.0,
                "rpk": 0.375,
                "mean": "3.5",
                "median": "3.2",
                "hitrate": "0.7",
                "variance": "0.3",
            },
        },
        "rpkscale": 1.5416666666666667,
        "output_file": "output.txt",
    }


def test_print_coverage_stats(kwargs):
    expected_output = [
        "sample1\tcontig1\t500\t500\t500\t0.5\t0.3243\t1.5\t1.2\t0.8\t0.1\n",
        "sample1\tcontig2\t1000\t1000\t666.7\t0.6667\t0.4324\t2.5\t2.2\t0.9\t0.2\n",
        "sample1\tcontig3\t750\t750\t375\t0.375\t0.2432\t3.5\t3.2\t0.7\t0.3\n",
    ]
    with patch("builtins.open", mock_open()) as mock_file:
        print_coverage_stats(**kwargs)
        mock_file.assert_called_once_with("output.txt", "w")
        handle = mock_file()
        handle.write.assert_has_calls(
            [call(line) for line in expected_output], any_order=False
        )
