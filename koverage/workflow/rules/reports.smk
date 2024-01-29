rule sample_tsv:
    output:
        tsv = os.path.join(dir["out"], "koverage.samples.tsv")
    params:
        sample_dict = samples["reads"]
    run:
        from metasnek import fastq_finder
        fastq_finder.write_samples_tsv(params.sample_dict, output.tsv)


rule coverage_report:
    """Generate html report for koverage"""
    input:
        smpl = os.path.join(dir["result"], "sample_coverage.tsv"),
        all = os.path.join(dir["result"], "all_coverage.tsv")
    output:
        html = os.path.join(dir["result"], "report.html")
    params:
        # pyspy = config["args"]["pyspy"],
        sample_cov_desc = config["report"]["map"]["sample_cov_desc"],
        all_cov_desc = config["report"]["map"]["all_cov_desc"],
        sample_names = samples["names"],
        ref_fasta = config["args"]["ref"],
        max_ctg = config["args"]["report_max_ctg"]
    threads: 1
    log:
        err = os.path.join(dir["log"], "coverage_report.err"),
        # pyspy = os.path.join(dir["log"], "coverage_report.svg")
    benchmark:
        os.path.join(dir["bench"], "coverage_report.txt")
    script:
        os.path.join(dir["scripts"], "koverageReport.py")


rule buildEnv:
    output:
        os.path.join(dir["temp"], "{env}.done")
    conda:
        lambda wildcards: os.path.join(dir["env"], wildcards.env)
    shell:
        "touch {output}"
