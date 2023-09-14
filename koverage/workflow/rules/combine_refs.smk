rule combine_fastas:
    """Combine the multiple ref fastas"""
    output:
        os.path.join(dir["temp"],"concatenated_refs.fasta")
    params:
        ref_files = references
    localrule:
        True
    run:
        fasta_finder.combine_fastas(params.ref_files, output[0])
