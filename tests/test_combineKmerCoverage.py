import gzip
import koverage.scripts.combineKmerCoverage as ckc


def test_collect_kmer_coverage_stats(tmp_path):
    input_file = tmp_path / "test_input.txt"
    with gzip.open(input_file, "wt") as f:
        f.write(
            "Sample\tContig\tSum\tMean\tMedian\tHitrate\tVariance\n"
            "sample\tcontig1\t5\t0.25\t0.25\t0.25\t0.9\n"
            "sample\tcontig1\t5\t0.25\t1\t0.5\t0\n"
            "sample\tcontig2\t20\t1.0\t2.5\t1.5\t1.8\n"
        )
    expected_result = {
        "contig1": {"sum": 10, "mean": 0.5, "median": 1.25},
        "contig2": {"sum": 20, "mean": 1.0, "median": 2.5},
    }
    result = ckc.collect_kmer_coverage_stats(input_file)
    assert result == expected_result


def test_print_kmer_coverage():
    output_file = "test_output.txt.gz"
    all_coverage = {
        "contig1": {"sum": 10, "mean": 0.5, "median": 1.25},
        "contig2": {"sum": 20, "mean": 1.0, "median": 2.5},
        "contig3": {"sum": 20, "mean": 1.0, "median": 2.5},
    }
    ckc.print_kmer_coverage(all_coverage, output_file, lines_per_batch=2)
    with gzip.open(output_file, "rt") as file:
        lines = file.readlines()
        assert len(lines) == 4
        assert lines[0] == "Contig\tSum\tMean\tMedian\n"
        assert lines[1] == "contig1\t10\t0.5\t1.25\n"
        assert lines[2] == "contig2\t20\t1\t2.5\n"
        assert lines[3] == "contig3\t20\t1\t2.5\n"
