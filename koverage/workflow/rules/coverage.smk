rule read_r1:
    """Read the R1 file"""
    input:
        lambda wildcards: samples.reads[wildcards.sample]["R1"]
    output:
        pipe(os.path.join(dir.temp, "{sample}.R1.fastq"))
    params:
        lambda wildcards: "zcat" if samples.reads[wildcards.sample]["R1"].endswith(".gz") else "cat"
    shell:
        """
        {params} {input} >> {output.fh2}
        """

if config.args.bams:
    rule sam_to_counts_bam:
        """Collect the counts for each contig from the piped SAM output"""
        input:
            os.path.join(dir.temp,"{sample}.sam"),
        output:
            tsv = temp(os.path.join(dir.temp, "{sample}.counts.tsv")),
            sam = pipe(os.path.join(dir.temp, "{sample}.save.sam"))
        params:
            config.args.bams
        script:
            os.path.join(dir.scripts, "samToCounts.py")
else:
    rule sam_to_counts:
        """Collect the counts for each contig from the piped SAM output"""
        input:
            os.path.join(dir.temp,"{sample}.sam"),
        output:
            tsv = temp(os.path.join(dir.temp, "{sample}.counts.tsv"))
        params:
            config.args.bams
        script:
            os.path.join(dir.scripts, "samToCounts.py")

rule save_bam:
    """Sort the reads and save the bam file from the piped SAM output"""
    input:
        os.path.join(dir.temp,"{sample}.save.sam"),
    output:
        os.path.join(dir.out,"{sample}.bam")
    threads:
        config.resources.bam.threads
    resources:
        mem_mb = config.resources.bam.mem_mb,
        time = config.resources.bam.time
    conda:
        os.path.join(dir.env, "samtools.yaml")
    shell:
        """samtools sort -@ {threads} {input} > {output}"""


rule sample_coverage:
    """convert raw counts to RPKM, FPKM, TPM, etc values"""
    input:
        tsv = os.path.join(dir.temp,"{sample}.counts.tsv"),
        r1 = os.path.join(dir.temp,"{sample}.R1.counts")
    output:
        temp(os.path.join(dir.temp,"{sample}.cov.tsv"))
    script:
        os.path.join(dir.scripts, "sampleCoverage.py")


rule combine_coverage:
    """Combine all sample coverages"""
    input:
        expand(os.path.join(dir.temp,"{sample}.cov.tsv"), sample=samples.names)
    output:
        sample_cov = os.path.join(dir.result, "sample_coverage.tsv"),
        all_cov = os.path.join(dir.result, "all_coverage.tsv"),
        sample_sum = os.path.join(dir.result, "sample_summary.tsv"),
        all_sum = os.path.join(dir.result, "all_summary.tsv")
    script:
        os.path.join(dir.scripts, "combineCoverage.py")

