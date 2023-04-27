rule sam_to_counts:
    """Collect the counts for each contig from the piped SAM output"""
    input:
        os.path.join(dir.temp,"{sample}.sam"),
    output:
        temp(os.path.join(dir.temp, "{sample}.counts.tsv"))
    run:
        pass


rule save_bam:
    """Sort the reads and save the bam file from the piped SAM output"""
    input:
        os.path.join(dir.temp,"{sample}.sam"),
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
    run:
        with open(input.r1, 'r') as f:
            lib = int(f.readline().strip()) / 1000000
        with open(output[0], 'w') as o:
            with open(input.tsv, 'r') as t:
                for line in t:
                    l = line.strip().split()
                    rpm = int(l[2]) / lib
                    rpkm = rpm / ( int(l[1]) / 1000 )
                    o.write(f"{l[0]}\t{rpkm}\n")


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


rule save_bam:
    """Save the BAM file"""
    input:
        os.path.join(dir.temp,"{sample}.bam")
    output:
        os.path.join(dir.out, "{sample}.bam")
    shell:
        "mv {input} {output}"
