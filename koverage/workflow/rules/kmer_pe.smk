rule jellyfish_db:
    """Calculate a jellyfish database of the reads"""
    input:
        r1=lambda wildcards: samples.reads[wildcards.sample]["R1"],
        r2=lambda wildcards: samples.reads[wildcards.sample]["R2"]
    output:
        os.path.join(dir.temp, "{sample}." + str(config.args.kmer_size) + "mer"),
    threads:
        config.resources.map.cpu
    resources:
        mem_mb = config.resources.map.mem_mb,
        time = config.resources.map.time_min
    params:
        kmer = config.args.kmer_size,
        jf = config.params.jellyfish,
        cat = lambda wildcards: "zcat" if samples.reads[wildcards.sample]["R1"].endswith(".gz") else "cat"
    conda:
        os.path.join(dir.env, "jellyfish.yaml")
    benchmark:
        os.path.join(dir.bench, "jellyfish_db.{sample}.txt")
    log:
        os.path.join(dir.log, "jellyfish.{sample}.err")
    shell:
        """
        jellyfish count \
            -m {params.kmer} \
            -t {threads} \
            -o {output} \
            {params.jf} \
            <({params.cat} {input.r1} {input.r2})
        """


rule ref_kmer_prep:
    """Sample kmers for a reference fasta"""
    input:
        config.args.assembly
    output:
        config.refkmers
    threads:
        config.resources.map.cpu
    benchmark:
        os.path.join(dir.bench, "ref_kmer_prep.txt")
    log:
        os.path.join(dir.log, "ref_kmer_prep.err")
    script:
        os.path.join(dir.scripts, "refSampleKmer.py")
