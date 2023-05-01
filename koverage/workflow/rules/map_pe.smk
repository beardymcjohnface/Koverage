import re

rule raw_coverage:
    """Map and collect the raw read counts for each sample"""
    input:
        assembly = config.args.assembly,
        r1=os.path.join(dir.temp, "{sample}.R1.fastq"),
        r2=lambda wildcards: samples.reads[wildcards.sample]["R2"]
    output:
        sam = pipe(os.path.join(dir.temp, "{sample}.sam")),
        r1 = os.path.join(dir.temp, "{sample}.R1.count")
    threads:
        config.resources.map.cpu
    resources:
        mem_mb = config.resources.map.mem_mb,
        time = config.resources.map.time_min
    params:
        minimap = "-ax sr --secondary=no",
        div = lambda wildcards: "4" if re.match(".*fastq$|.*fastq.gz$", samples.reads[wildcards.sample]["R1"]) else "2"
    conda:
        os.path.join(dir.env, "minimap.yaml")
    log:
        os.path.join(dir.log, "{sample}.minimap2.err")
    group:
        "pipejob"
    shell:
        """
        minimap2 -t {threads} {params.minimap} {input.assembly} \
            <( cat {input.r1} | tee >( wc -l | awk '{{ print $1 / {params.div} }}' > {output.r1} ) \
            {input.r2} 2> {log} \
            | samtools sort -@ {threads} >> {output.sam}
        """
