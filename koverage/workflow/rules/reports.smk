rule sample_tsv:
    output:
        tsv = os.path.join(dir.out, "koverage.samples.tsv")
    params:
        sample_dict = samples.reads
    run:
        from metasnek import fastq_finder
        fastq_finder.write_samples_tsv(params.sample_dict,output.tsv)