rule idx_ref:
    """Map and collect the raw read counts for each sample"""
    input:
        config.args.ref
    output:
        config.args.ref + '.idx'
    threads:
        config.resources.map.cpu
    resources:
        mem_mb = config.resources.map.mem_mb,
        time = config.resources.map.time_min
    conda:
        os.path.join(dir.env, "minimap.yaml")
    benchmark:
        os.path.join(dir.bench, "idx_ref.txt")
    log:
        os.path.join(dir.log, "idx_ref.err")
    shell:
        """
        minimap2 -t {threads} -d {output} {input} 2> {log}
        """


rule raw_coverage:
    """Map and collect the raw read counts for each sample
    
    lib: single line, total number of reads
    counts: "contig\tcontig_len\tcount\tmean\tmedian\thitrate\tvariance
    """
    input:
        ref = config.args.ref + ".idx",
        r1=lambda wildcards: samples.reads[wildcards.sample]["R1"],
    output:
        lib = temp(os.path.join(dir.temp, "{sample}.lib")),
        counts = temp(os.path.join(dir.temp, "{sample}.counts.tsv")),
    threads:
        config.resources.map.cpu
    resources:
        mem_mb = config.resources.map.mem_mb,
        time = config.resources.map.time_min
    params:
        r2 = lambda wildcards: samples.reads[wildcards.sample]["R2"],
        pafs = config.args.pafs,
        paf_dir = dir.paf,
        max_depth = config.args.max_depth,
        bin_width = config.args.bin_width,
        minimap = config.args.minimap,
        pyspy = config.args.pyspy
    conda:
        os.path.join(dir.env, "minimap.yaml")
    benchmark:
        os.path.join(dir.bench, "raw_coverage.{sample}.txt")
    log:
        err = os.path.join(dir.log, "raw_coverage.{sample}.err"),
        pyspy = os.path.join(dir.log, "raw_coverage.{sample}.svg")
    script:
        os.path.join(dir.scripts, "minimapWrapper.py")


rule sample_coverage:
    """convert raw counts to coverage values
    
    output: sample\tcontig\tCount\tRPM\tRPKM\tRPK\tTPM\tMean\tMedian\tHitrate\tVariance
    """
    input:
        counts = os.path.join(dir.temp,"{sample}.counts.tsv"),
        lib = os.path.join(dir.temp,"{sample}.lib"),
    output:
        temp(os.path.join(dir.temp,"{sample}.cov.tsv"))
    params:
        pyspy = config.args.pyspy
    threads: 1
    log:
        err =os.path.join(dir.log, "sample_coverage.{sample}.err"),
        pyspy = os.path.join(dir.log, "sample_coverage.{sample}.svg")
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
        err = os.path.join(dir.log, "all_sample_coverage.err"),
        pyspy = os.path.join(dir.log, "all_sample_coverage.err")
    benchmark:
        os.path.join(dir.bench, "all_sample_coverage.txt")
    shell:
        """
        printf "Sample\tContig\tCount\tRPM\tRPKM\tRPK\tTPM\tMean\tMedian\tHitrate\tVariance\n" > {output} 2> {log}
        cat {input} >> {output} 2> {log}
        """


rule combine_coverage:
    """Combine all sample coverages"""
    input:
        os.path.join(dir.result,"sample_coverage.tsv")
    output:
        all_cov = os.path.join(dir.result, "all_coverage.tsv"),
        # sample_sum = os.path.join(dir.result, "sample_summary.tsv"),
        # all_sum = os.path.join(dir.result, "all_summary.tsv")
    params:
        pyspy = config.args.pyspy
    threads: 1
    log:
        err = os.path.join(dir.log, "combine_coverage.err"),
        pyspy = os.path.join(dir.log, "combine_coverage.svg")
    benchmark:
        os.path.join(dir.bench, "combine_coverage.txt")
    script:
        os.path.join(dir.scripts, "combineCoverage.py")

