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
        config.args.ref
    output:
        config.refkmers
    threads:
        config.resources.map.cpu
    resources:
        mem_mb = config.resources.map.mem_mb,
        time = config.resources.map.time_min
    params:
        ksize = config.args.kmer_size,
        kspace = config.args.kmer_sample,
        kmin = config.args.kmer_min,
        kmax = config.args.kmer_max,
    benchmark:
        os.path.join(dir.bench, "ref_kmer_prep.txt")
    log:
        os.path.join(dir.log, "ref_kmer_prep.err")
    script:
        os.path.join(dir.scripts, "refSampleKmer.py")


rule kmer_screen:
    """Screen jellyfish database for ref kmers"""
    input:
        ref = config.refkmers,
        db = os.path.join(dir.temp, "{sample}." + str(config.args.kmer_size) + "mer")
    output:
        os.path.join(dir.temp, "{sample}." + str(config.args.kmer_size) + "mer.kcov.zst")
    threads:
        config.resources.jf.cpu
    resources:
        mem_mb = config.resources.jf.mem_mb,
        time = config.resources.jf.time_min
    conda:
        os.path.join(dir.env,"jellyfish.yaml")
    benchmark:
        os.path.join(dir.bench, "kmer_screen.{sample}.txt")
    log:
        os.path.join(dir.log, "kmer_screen.{sample}.err")
    script:
        os.path.join(dir.scripts, "kmerScreen.py")
