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


rule sam_to_counts:
    """Collect the counts for each contig from the piped SAM output"""
    input:
        os.path.join(dir.temp,"{sample}.sam"),
    output:
        tsv = temp(os.path.join(dir.temp, "{sample}.counts.tsv")),
        sam = pipe(os.path.join(dir.temp, "{sample}.depth.sam"))
    script:
        os.path.join(dir.scripts, "samToCounts.py")


if config.args.bams:
    rule mpileup_save_bam:
        input:
            os.path.join(dir.temp,"{sample}.depth.sam")
        output:
            mp = pipe(os.path.join(dir.temp,"{sample}.mpileup")),
            bm = os.path.join(dir.bam,"{sample}.bam")
        conda:
            os.path.join(dir.env, "samtools.yaml")
        shell:
            """
            samtools view -b {input} \
                | tee {output.bm} \
                | samtools mpileup >> {output.mp}
            """
else:
    rule mpileup:
        input:
            os.path.join(dir.temp,"{sample}.depth.sam")
        output:
            pipe(os.path.join(dir.temp,"{sample}.mpileup"))
        conda:
            os.path.join(dir.env, "samtools.yaml")
        shell:
            """
            samtools view -b {input} \
                | samtools mpileup -Aa - \
                | cut -f1,4 >> {output}
            """


rule mpileup_to_depth:
    """Collect depth histograms for each contig for sample and echo sam output"""
    input:
        os.path.join(dir.temp,"{sample}.mpileup")
    output:
        os.path.join(dir.temp,"{sample}.depth.tsv")
    params:
        config.args.bams
    script:
        os.path.join(dir.scripts, "mpileupToDepth.py")


rule sample_coverage:
    """convert raw counts to RPKM, FPKM, TPM, etc values"""
    input:
        tsv = os.path.join(dir.temp,"{sample}.counts.tsv"),
        r1 = os.path.join(dir.temp,"{sample}.R1.counts")
    output:
        temp(os.path.join(dir.temp,"{sample}.cov.tsv"))
    script:
        os.path.join(dir.scripts, "sampleCoverage.py")


rule all_sample_coverage:
    """Concatenate the sample coverage TSVs"""
    input:
        expand(os.path.join(dir.temp,"{sample}.cov.tsv"), sample=samples.names)
    output:
        os.path.join(dir.result,"sample_coverage.tsv")
    shell:
        """
        printf "Sample\tContig\tRPM\tRPKM\tRPK\tTPM\n" > {output}
        cat {input} >> {output}
        """


rule combine_coverage:
    """Combine all sample coverages"""
    input:
        os.path.join(dir.result,"sample_coverage.tsv")
    output:
        all_cov = os.path.join(dir.result, "all_coverage.tsv"),
        sample_sum = os.path.join(dir.result, "sample_summary.tsv"),
        all_sum = os.path.join(dir.result, "all_summary.tsv")
    script:
        os.path.join(dir.scripts, "combineCoverage.py")
