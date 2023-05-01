rule read_r1:
    """Read the R1 file"""
    input:
        lambda wildcards: samples.reads[wildcards.sample]["R1"]
    output:
        os.path.join(dir.temp, "{sample}.R1.fastq")
        # pipe(os.path.join(dir.temp, "{sample}.R1.fastq"))
    threads:
        config.resources.pipe.cpu
    resources:
        mem_mb = config.resources.pipe.mem_mb,
        time = config.resources.pipe.time_min
    params:
        lambda wildcards: "zcat" if samples.reads[wildcards.sample]["R1"].endswith(".gz") else "cat"
    group:
        "pipejob"
    log:
        os.path.join(dir.log, "read_r1.{sample}.err")
    shell:
        """
        {params} {input} >> {output} 2> {log}
        """


rule sam_to_counts:
    """Collect the counts for each contig from the piped SAM output"""
    input:
        os.path.join(dir.temp,"{sample}.sam"),
    output:
        tsv = temp(os.path.join(dir.temp, "{sample}.counts.tsv")),
        sam = os.path.join(dir.temp, "{sample}.depth.sam")
        # sam = pipe(os.path.join(dir.temp, "{sample}.depth.sam"))
    threads:
        config.resources.pipe.cpu
    resources:
        mem_mb = config.resources.pipe.mem_mb,
        time = config.resources.pipe.time_min
    group:
        "pipejob"
    log:
        os.path.join(dir.log, "sam_to_counts.{sample}.err")
    script:
        os.path.join(dir.scripts, "samToCounts.py")


if config.args.bams:
    rule mpileup_save_bam:
        input:
            os.path.join(dir.temp,"{sample}.depth.sam")
        output:
            mp = os.path.join(dir.temp,"{sample}.depth"),
            # mp = pipe(os.path.join(dir.temp,"{sample}.depth")),
            bm = os.path.join(dir.bam,"{sample}.bam")
        threads:
            config.resources.pipe.cpu
        resources:
            mem_mb=config.resources.pipe.mem_mb,
            time=config.resources.pipe.time_min
        conda:
            os.path.join(dir.env, "minimap.yaml")
        group:
            "pipejob"
        log:
            os.path.join(dir.log,"mpileup.{sample}.err")
        shell:
            """
            {{
            samtools view -b {input} \
                | tee {output.bm} \
                | samtools depth -aa - >> {output.mp};
            }} 2> {log}
            """
else:
    rule mpileup:
        input:
            os.path.join(dir.temp,"{sample}.depth.sam")
        output:
            os.path.join(dir.temp,"{sample}.depth")
            # pipe(os.path.join(dir.temp,"{sample}.depth"))
        threads:
            config.resources.pipe.cpu
        resources:
            mem_mb=config.resources.pipe.mem_mb,
            time=config.resources.pipe.time_min
        conda:
            os.path.join(dir.env, "minimap.yaml")
        group:
            "pipejob"
        log:
            os.path.join(dir.log,"mpileup.{sample}.err")
        shell:
            """
            {{
            samtools view -b {input} \
                | samtools depth -aa - >> {output};
            }} 2> {log}
            """


rule mpileup_to_depth:
    """Collect depth histograms for each contig for sample and echo sam output"""
    input:
        os.path.join(dir.temp,"{sample}.depth")
    output:
        hist = temp(os.path.join(dir.temp,"{sample}.depth.tsv")), # todo: keep or delete?
        kurt = temp(os.path.join(dir.temp, "{sample}.kurtosis.tsv"))
    threads:
        config.resources.pipe.cpu
    resources:
        mem_mb = config.resources.pipe.mem_mb,
        time = config.resources.pipe.time_min
    params:
        max_depth = config.args.max_depth
    group:
        "pipejob"
    log:
        os.path.join(dir.log, "mpileup_to_depth.{sample}.err")
    script:
        os.path.join(dir.scripts, "mpileupToDepth.py")


rule sample_coverage:
    """convert raw counts to RPKM, FPKM, TPM, etc values"""
    input:
        tsv = os.path.join(dir.temp,"{sample}.counts.tsv"),
        r1 = os.path.join(dir.temp,"{sample}.R1.count"),
        kurt = os.path.join(dir.temp, "{sample}.kurtosis.tsv")
    output:
        temp(os.path.join(dir.temp,"{sample}.cov.tsv"))
    threads: 1
    log:
        os.path.join(dir.log, "sample_coverage.{sample}.err")
    script:
        os.path.join(dir.scripts, "sampleCoverage.py")


rule all_sample_coverage:
    """Concatenate the sample coverage TSVs"""
    input:
        expand(os.path.join(dir.temp,"{sample}.cov.tsv"), sample=samples.names)
    output:
        os.path.join(dir.result,"sample_coverage.tsv")
    threads: 1
    log:
        os.path.join(dir.log, "all_sample_coverage.err")
    shell:
        """
        printf "Sample\tContig\tRPM\tRPKM\tRPK\tTPM\tKurtosis\n" > {output} 2> {log}
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
    script:
        os.path.join(dir.scripts, "combineCoverage.py")
