rule jellyfish_db:
    """Calculate a jellyfish database of the reads"""
    input:
        r1=lambda wildcards: samples.reads[wildcards.sample]["R1"],
    output:
        os.path.join(dir.temp, "{sample}." + str(config.args.kmer_size) + "mer"),
    threads:
        config.resources.map.cpu
    resources:
        mem_mb = config.resources.map.mem_mb,
        time = config.resources.map.time_min
    params:
        r2 = lambda wildcards: samples.reads[wildcards.sample]["R2"],
        kmer = config.args.kmer_size,
        jf = config.params.jellyfish,
        cat = lambda wildcards: "gunzip -c" if samples.reads[wildcards.sample]["R1"].endswith(".gz") else "cat"
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
            <({params.cat} {input.r1} {params.r2})
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
        temp(os.path.join(dir.temp, "{sample}." + str(config.args.kmer_size) + "mer.kcov.zst"))
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


rule all_sample_kmer_coverage:
    """Concatenate the sample coverage TSVs"""
    input:
        expand(os.path.join(dir.temp, "{sample}." + str(config.args.kmer_size) + "mer.kcov.zst"), sample=samples.names)
    output:
        os.path.join(dir.result, "sample_kmer_coverage." + str(config.args.kmer_size) + "mer.tsv.gz")
    threads: 1
    conda:
        os.path.join(dir.env, "zstd.yaml")
    log:
        os.path.join(dir.log, "all_sample_kmer_coverage.err")
    benchmark:
        os.path.join(dir.bench, "all_sample_kmer_coverage.txt")
    shell:
        """
        {{
            printf "Sample\tContig\tMean\tMedian\tHitrate\tVariance\n" 2> {log};
            zstdcat {input} 2> {log};
        }} | gzip -1 - > {output}
        """


rule combine_kmer_coverage:
    """Combine all sample kmer coverages"""
    input:
        os.path.join(dir.result, "sample_kmer_coverage." + str(config.args.kmer_size) + "mer.tsv.gz")
    output:
        all_cov = os.path.join(dir.result, "all_kmer_coverage.tsv"),
        # sample_sum = os.path.join(dir.result, "sample_summary.tsv"),
        # all_sum = os.path.join(dir.result, "all_summary.tsv")
    threads: 1
    log:
        os.path.join(dir.log, "combine_kmer_coverage.err")
    benchmark:
        os.path.join(dir.bench, "combine_kmer_coverage.txt")
    script:
        os.path.join(dir.scripts, "combineKmerCoverage.py")
