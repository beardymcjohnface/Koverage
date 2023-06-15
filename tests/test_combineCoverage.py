import koverage.scripts.combineCoverage as cc


def test_collect_coverage_stats(tmp_path):
    input_file = tmp_path / "test_input.txt"
    input_file.write_text(
        "Sample\tContig\tCount\tRPM\tRPKM\tRPK\tTPM\tHitrate\tVariance\n"
        "sample\tcontig1\t5\t0.25\t0.25\t0.25\t0.9\n"
        "sample\tcontig1\t5\t0.25\t1\t0.5\t0\n"
        "sample\tcontig2\t20\t1.0\t2.5\t1.5\t1.8\n"
    )
    expected_result = {
        "contig1": {"count": 10, "rpm": 0.5, "rpkm": 1.25, "rpk": 0.75, "tpm": 0.9},
        "contig2": {"count": 20, "rpm": 1.0, "rpkm": 2.5, "rpk": 1.5, "tpm": 1.8},
    }
    result = cc.collect_coverage_stats(input_file)
    assert result == expected_result


def test_print_sample_coverage():
    output_file = "test_output.txt"
    all_coverage = {
        "contig1": {"count": 10, "rpm": 0.5, "rpkm": 1.25, "rpk": 0.75, "tpm": 0.9},
        "contig2": {"count": 20, "rpm": 1.0, "rpkm": 2.5, "rpk": 1.5, "tpm": 1.8},
    }
    cc.print_sample_coverage(output_file, all_coverage)
    with open(output_file, "r") as file:
        lines = file.readlines()
        assert len(lines) == 3
        assert lines[0] == "Contig\tCount\tRPM\tRPKM\tRPK\tTPM\n"
        assert lines[1] == "contig1\t10\t0.5\t1.25\t0.75\t0.9\n"
        assert lines[2] == "contig2\t20\t1\t2.5\t1.5\t1.8\n"
