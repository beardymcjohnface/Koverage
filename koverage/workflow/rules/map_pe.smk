rule raw_coverage: # todo: convert the awk bit to another named pipe?
    """Map and collect the raw read counts for each sample"""
    input:
        assembly = config.args.assembly,
        r1=lambda wildcards: samples.reads[wildcards.sample]["R1"],
        r2=lambda wildcards: samples.reads[wildcards.sample]["R2"],
    output:
        sam = pipe(os.path.join(dir.temp, "{sample}.sam")),
        r1 = temp(os.path.join(dir.temp, "{sample}.R1.counts")),
    threads:
        config.resources.map.threads
    resources:
        mem_mb = config.resources.map.mem_mb,
        time = config.resources.map.time_min
    params:
        minimap = "-ax sr --secondary=no",
    conda:
        os.path.join(dir.env, "minimap.yaml")
    log:
        os.path.join(dir.log, "{sample}.raw_coverage.err")
    shell:
        """
        minimap2 -t {threads} {params.minimap} {input.assembly} \
            <(zcat -f {input.r1} | tee >( wc -l | awk '{{print $1 / 4}}' > {output.r1})) \
            <(zcat -f {input.r2})) \
            2> {log} > {output.sam}
        """
