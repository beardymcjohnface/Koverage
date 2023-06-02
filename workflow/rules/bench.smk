rule bench_map_pe:
    input:
        ref = config.args.ref,
        r1=lambda wildcards: samples.reads[wildcards.sample]["R1"],
    output:
        bam = os.path.join(dir.temp, "{sample}.bam"),
        bai = os.path.join(dir.temp, "{sample}.bam.bai")
    params:
        r2 = lambda wildcards: samples.reads[wildcards.sample]["R2"]
    threads:
        config.resources.map.cpu
    resources:
        mem_mb = config.resources.map.mem_mb,
        time = config.resources.map.time_min
    conda:
        os.path.join(dir.env, "minimap.yaml")
    benchmark:
        os.path.join(dir.bench, "bench_map_pe.{sample}.txt")
    log:
        os.path.join(dir.log, "bench_map_pe.{sample}.err")
    shell:
        """
        minimap2 -t {threads} -ax sr --secondary=no {input.ref} {input.r1} {params.r2} \
            | samtools sort -@ {threads} - > {output.bam}
        samtools index {output.bam}
        """


rule bench_depth:
    input:
        os.path.join(dir.temp, "{sample}.bam")
    output:
        os.path.join(dir.temp, "{sample}.depth")
    threads:
        config.resources.map.cpu
    resources:
        mem_mb = config.resources.map.mem_mb,
        time = config.resources.map.time_min
    conda:
        os.path.join(dir.env, "minimap.yaml")
    benchmark:
        os.path.join(dir.bench, "bench_depth.{sample}.txt")
    log:
        os.path.join(dir.log, "bench_depth.{sample}.err")
    shell:
        """
        samtools depth -aa -@ {threads} {input} > {output}
        """


rule bench_bam2counts:
    input:
        os.path.join(dir.temp, "{sample}.bam")
    output:
        os.path.join(dir.temp, "{sample}.cov")
    conda:
        os.path.join(dir.env, "coverm.yaml")
    shell:
        """
        coverm contig -b {input} -m count -m rpkm -m tpm -m covered_bases -m variance > {output}
        """


rule bench_combine:
    input:
        expand(os.path.join(dir.temp, "{sample}.cov"), sample=samples.names)
    output:
        os.path.join(dir.result, "sample_bench_coverage.tsv")
    params:
        samples = samples.names,
        dir = dir.temp
    run:
        with open(output[0], "w") as outfh:
            outfh.write("Sample\tContig\tCount\tRPKM\tTPM\tCovered_bases\tVariance\n")
            for sample in params.samples:
                with open(os.path.join(params.dir, f"{sample}.cov"), "r") as infh:
                    for line in infh:
                        if not line.startswith("Contig\t"):
                            outfh.write(f"{sample}\t{line}")