
rule raw_coverage:
    """Map and collect the raw read counts for each sample"""
    input:
        assembly = config.args.assembly,
        r1=lambda wildcards: samples.reads[wildcards.sample]["R1"],
        r2=lambda wildcards: samples.reads[wildcards.sample]["R2"]
    output:
        lib = temp(os.path.join(dir.temp, "{sample}.lib")),
        var = temp(os.path.join(dir.temp, "{sample}.variance.tsv")),
        counts = temp(os.path.join(dir.temp, "{sample}.counts.tsv")),
        bamfile = os.path.join(dir.bam,"{sample}.bam"),
    threads:
        config.resources.map.cpu
    resources:
        mem_mb = config.resources.map.mem_mb,
        time = config.resources.map.time_min
    params:
        bams = config.args.bams,
        max_depth = config.args.max_depth,
        bin_width = config.args.bin_width
    conda:
        os.path.join(dir.env, "minimap.yaml")
    benchmark:
        os.path.join(dir.bench, "raw_coverage.{sample}.txt")
    log:
        os.path.join(dir.log, "minimap2.{sample}.err")
    script:
        os.path.join(dir.scripts, "minimapWrapperPaf.py")
