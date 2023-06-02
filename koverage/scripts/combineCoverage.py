#!/usr/bin/env python3

import logging


def collect_coverage_stats(input_file):
    """Combine the mapped coverage stats for all samples.

    :param input_file: Text TSV file (Sample\tContig\tCount\tRPM\tRPKM\tRPK\tTPM\tHitrate\tVariance)
    :returns all_coverage: Dictionary of counts for each contig (dict[contigID]["count"/"rpm"/"rpkm"/"rpk"/"tpm"])
    """
    all_coverage = {}
    with open(input_file, "r") as infh:
        infh.readline()
        for line in infh:
            l = line.strip().split()
            try:
                assert(type(all_coverage[l[1]]) is dict)
            except (AssertionError, KeyError):
                all_coverage[l[1]] = {"count": 0, "rpm": 0, "rpkm": 0, "rpk": 0, "tpm": 0}
            all_coverage[l[1]]["count"] += int(l[2])
            all_coverage[l[1]]["rpm"] += float(l[3])
            all_coverage[l[1]]["rpkm"] += float(l[4])
            all_coverage[l[1]]["rpk"] += float(l[5])
            all_coverage[l[1]]["tpm"] += float(l[6])
    return all_coverage


def print_sample_coverage(output_file, all_coverage):
    """Print the combined coverage statistics from collect_coverage_stats().

    :param output_file: Text TSV filepath for writing
    :param all_coverage: Dictionary of counts for each contig (dict[contigID]["count"/"rpm"/"rpkm"/"rpk"/"tpm"])
    :returns: None
    """
    with open(output_file, "w") as outCov:
        outCov.write("Contig\tCount\tRPM\tRPKM\tRPK\tTPM\n")
        for contig in sorted(all_coverage.keys()):
            outCov.write("\t".join([
                contig,
                str(all_coverage[contig]["count"]),
                "{:.{}g}".format(all_coverage[contig]["rpm"], 4),
                "{:.{}g}".format(all_coverage[contig]["rpkm"], 4),
                "{:.{}g}".format(all_coverage[contig]["rpk"], 4),
                "{:.{}g}".format(all_coverage[contig]["tpm"], 4) + "\n"
            ]))


def main(input_file, output_file, log_file):
    logging.basicConfig(filename=log_file, filemode="w", level=logging.DEBUG)
    logging.debug("Collecting combined coverage stats")
    all_coverage = collect_coverage_stats(input_file)
    logging.debug("Printing all sample coverage")
    print_sample_coverage(output_file, all_coverage)


if __name__ == "__main__":
    main(snakemake.input[0], snakemake.output.all_cov, snakemake.log[0])
