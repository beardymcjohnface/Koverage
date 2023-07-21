---
title: "Koverage: Read-coverage analysis for massive (meta)genomics datasets"
tags:
  - Python
  - Snakemake
  - Genomics
  - Metagenomics
authors:
  - name: Michael J. Roach
    orcid: 0000-0003-1488-5148
    affiliation: 1
    corresponding: true
  - name: Brad Hart
    orcid: 
    affiliation: 2
  - name: Sarah Beecroft
    orcid: 
    affiliation: 3
  - name: Bhavya Papudeshi
    orcid:
    affiliation: 1
  - name: Laura Inglis
    orcid: 
    affiliation: 1
  - name: Susie Grigson
    orcid:
    affiliation: 1
  - name: Vijini Mallawaarachchi
    orcid:
    affiliation: 1
  - name: George Bouras
    orcid:
    affiliation: "4,5"
  - name: Robert A. Edwards
    orcid:
    affiliation: 1
affiliations:
  - name: "Flinders Accelerator for Microbiome Exploration, Flinders University, Adelaide, SA, Australia"
    index: 1
  - name: "Health and Biomedical Innovation, Clinical and Health Sciences, University of South Australia, SA, Australia"
    index: 2
  - name: "Harry Perkins Institute of Medical Research, Perth, WA, Australia"
    index: 3
  - name: "Adelaide Medical School, Faculty of Health and Medical Sciences, The University of Adelaide, Adelaide, South Australia 5005, Australia"
    index: 4
  - name: "The Department of Surgery – Otolaryngology Head and Neck Surgery, Central Adelaide Local Health Network, Adelaide, South Australia 5000, Australia"
    index: 5
date:
bibliography: paper.bib
---

# Summary

Calculating the read-coverage of sequencing reads to a reference genome is a routine task with many applications. Some 
examples include identifying duplication or deletion events in a draft genome assembly, identifying related contigs for 
binning in metagenome assemblies, or analysing taxonomic compositions of metagenome samples. Calculating read-coverage 
information typically involves several read and write operations of the sequencing data. This is not a problem for small
datasets. However, this can be a significant bottleneck when analysing a large number of samples, or when screening very
large reference sequence files. Koverage is designed to reduce the I/O burden as much as possible to enable maximum 
scalability for large sample sizes. It also includes a kmer-based coverage method that significantly reduces the 
computational complexity of screening large reference genomes. Koverage is a Snakemake [@snakemake] based pipeline, 
providing out-of-the-box support for HPC and cloud environments. It utilises the Snaketool [@snaketool] command line 
interface and is available to install via PIP or Conda for maximum ease-of-use. The source code and documentation is 
available at [https://github.com/beardymcjohnface/Koverage](https://github.com/beardymcjohnface/Koverage).


# Statement of need

With the current state of sequencing technologies, it is trivial to generate terabytes of sequencing data for
hundreds or even thousands of samples. Furthermore, databases such as the Sequence Read Archive (SRA) and the European
Nucleotide Archive (ENV), containing nearly 100 petabytes combined of sequencing data, are constantly being mined and
reanalysed in bioinformatics analyses. Computational inefficiencies at such scales can translate into thousands of
dollars worth of service units, while memory and I/O bottlenecks can lead to under-utilisation of CPUs. In more severe
cases, I/O heavy processes in large parallel batches can result in significantly impaired performance, especially for
HPC clusters with a shared scratch space of spinning disk hard drives.

While there are tools for performing coverage calculations, they are not well optimised for deployment at large scales,
or when analysing large reference files. A typical approach may require several complete read and write operations of
the sequencing data in order to generate the coverage statistics. Furthermore, mapping to very large reference sequence
files can require large amounts of memory, or alternatively, aligning reads in chunks and coalescing these chunked
alignments at the end, resulting in even more I/O operations. Some solutions involve moving I/O operations into memory,
for instance via tempfs. However, whether this is a feasable option is highly system-dependent and will nevertheless
exacerbate any existing memory bottlenecks.

Koverage addresses the I/O bottleneck of large datasets by eliminating the sorting, reading, and writing of intermediate
alignment files. Koverage also includes a kmer-based implementation to eliminate memory bottlenecks that may arise from
screening large reference files. 

# Implementation

Koverage is written in Snakemake [@snakemake] and Python, and uses the Snaketool [@snaketool] command line interface 
(CLI). The Snaketool CLI will take the user input command line arguments and Koverage's default configuration to 
build a runtime config file. It will then build the Snakemake command and run the pipeline. Any unrecognised command 
line arguments are assumed to be Snakemake args and are added to the Snakemake command. For cluster or cloud execution, 
users are encouraged to generate a Snakemake profile for their chosen deployment, and Koverage has been designed to be 
compatible with Snakemake's Cookiecutter [@cookiecutter] template profiles. The only required inputs are the reference 
FASTA-format file (--ref), and the sample reads (--reads).

# Sample parsing

Koverage will parse sample reads (--reads) using MetaSnek fastq_finder [@metasnek]. Users supply either a directory 
containing their sequencing reads, or a tab-separated values (TSV) file listing their sample names and corresponding 
sequencing read filepaths. If users supply a directory to --reads, sample names and read file pairs will be inferred 
from the file names. If users supply a TSV file, sample names and filepaths will simply be read from the file. More 
information and examples are available at [https://gist.github.com/beardymcjohnface/bb161ba04ae1042299f48a4849e917c8](https://gist.github.com/beardymcjohnface/bb161ba04ae1042299f48a4849e917c8)

# Mapping-based coverage

This is the default method for calculating coverage statistics. Reads are mapped sample-by-sample to the reference 
genome using Minimap2 [@minimap]. The minimap2 alignments are parsed in real-time by a wrapper script that collects the
counts per contig and total counts per sample. Koverage also uses the read mapping coordinates to collect read counts 
for 'windows' or 'bins' along the contig. This allows for a fast approximation of the coverage of each contig by at 
least one read (hitrate), and of the evenness of coverage (variance) for each contig. Following mapping, the final 
counts, mean, median, hitrate, and variance are written to a TSV file. A second script calculates the Reads Per Million
(RPM), Reads Per Kilobase Million (RPKM), Reads Per Kilobase (RPK), and Transcripts Per Million (TPM) like so:

Method | Calculation
--- | ---
RPM | N/(T/1,000,000)
RPKM | (N/(T/1,000,000))/L
RPK | N/L
TPM | RPK / (R/1,000,000)

Where:
 - N = number of reads mapped to the contig
 - T = Total number of mapped reads for that sample
 - L = length of contig in kilobases
 - R = sum of RPK values for that sample

Lasty, the coverage from all samples are collated, and a summary of the coverage for 
each contig by all samples is calculated.

# Kmer-based coverage



# CoverM wrapper



# Acknowledgments



# References

