rule combine_fastas:
    """Combine the multiple ref fastas"""
    output:
        os.path.join(dir["temp"],"concatenated_refs.fasta")
    localrule:
        True
    run:
        fasta_finder.combine_fastas(references, output[0])
