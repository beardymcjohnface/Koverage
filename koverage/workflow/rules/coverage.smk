rule sample_coverage:
    """convert raw counts to RPKM, FPKM, TPM, etc values"""
    input:
        counts = os.path.join(dir.temp,"{sample}.counts.tsv"),
        lib = os.path.join(dir.temp,"{sample}.lib"),
        var = os.path.join(dir.temp, "{sample}.variance.tsv")
    output:
        temp(os.path.join(dir.temp,"{sample}.cov.tsv"))
    threads: 1
    log:
        os.path.join(dir.log, "sample_coverage.{sample}.err")
    benchmark:
        os.path.join(dir.bench, "sample_coverage.{sample}.txt")
    script:
        os.path.join(dir.scripts, "sampleCoverage.py")


rule all_sample_coverage:
    """Concatenate the sample coverage TSVs"""
    input:
        expand(os.path.join(dir.temp,"{sample}.cov.tsv"), sample=samples.names)
    output:
        os.path.join(dir.result, "sample_coverage.tsv")
    threads: 1
    log:
        os.path.join(dir.log, "all_sample_coverage.err")
    benchmark:
        os.path.join(dir.bench, "all_sample_coverage.txt")
    shell:
        """
        printf "Sample\tContig\tCount\tRPM\tRPKM\tRPK\tTPM\tHitrate\tVariance\n" > {output} 2> {log}
        cat {input} >> {output} 2> {log}
        """


rule combine_coverage:
    """Combine all sample coverages"""
    input:
        os.path.join(dir.result,"sample_coverage.tsv")
    output:
        all_cov = os.path.join(dir.result, "all_coverage.tsv"),
        sample_sum = os.path.join(dir.result, "sample_summary.tsv"),
        all_sum = os.path.join(dir.result, "all_summary.tsv")
    threads: 1
    log:
        os.path.join(dir.log, "combine_coverage.err")
    benchmark:
        os.path.join(dir.bench, "combine_coverage.txt")
    script:
        os.path.join(dir.scripts, "combineCoverage.py")

