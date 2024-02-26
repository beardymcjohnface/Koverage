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
    affiliation: "1,2"
    corresponding: true
  - name: Bradley J. Hart
    orcid: 0000-0001-8110-2460
    affiliation: 3
  - name: Sarah J. Beecroft
    orcid: 0000-0002-3935-2279
    affiliation: 4
  - name: Bhavya Papudeshi
    orcid: 0000-0001-5359-3100
    affiliation: 1
  - name: Laura K. Inglis
    orcid: 0000-0001-7919-8563
    affiliation: 1
  - name: Susanna R. Grigson
    orcid: 0000-0003-4738-3451
    affiliation: 1
  - name: Vijini Mallawaarachchi
    orcid: 0000-0002-2651-8719
    affiliation: 1
  - name: George Bouras
    orcid: 0000-0002-5885-4186
    affiliation: "5,6"
  - name: Robert A. Edwards
    orcid: 0000-0001-8383-8949
    affiliation: 1
affiliations:
  - name: "Flinders Accelerator for Microbiome Exploration, Flinders University, Adelaide, SA, Australia"
    index: 1
  - name: "Adelaide Centre for Epigenetics and the South Australian Immunogenomics Cancer Institute, Faculty of Health and Medical Sciences, The University of Adelaide, Adelaide, SA, Australia"
    index: 2
  - name: "Health and Biomedical Innovation, Clinical and Health Sciences, University of South Australia, SA, Australia"
    index: 3
  - name: "Pawsey Supercomputing Research Centre, Kensington, WA, Australia"
    index: 4
  - name: "Adelaide Medical School, Faculty of Health and Medical Sciences, The University of Adelaide, Adelaide, SA, Australia"
    index: 5
  - name: "The Department of Surgery – Otolaryngology Head and Neck Surgery, Central Adelaide Local Health Network, Adelaide, SA, Australia"
    index: 6
date:
bibliography: paper.bib
---

# Summary

Genomes of organisms are constructed by assembling sequencing from whole genome sequencing (WGS). It is useful to determine sequencing read-coverage of the genome assembly, for instance identifying duplication or deletion events, related contigs for binning metagenomes [@metacoag;@graphbin2], or analysing taxonomic compositions of metagenomes [@condiga]. Although calculating read-coverage is a routine task, it typically involves several complete read and write operations (I/O operations). This is not a problem for small datasets, but can be a significant bottleneck for very large datasets. Koverage reduces I/O burden as much as possible to enable maximum scalability. Koverage includes a kmer-based method that significantly reduces the computational complexity for very large reference genomes. Koverage uses Snakemake [@snakemake], providing out-of-the-box support for HPC and cloud environments. It utilises the Snaketool [@snaketool] command line interface, and is installable with PIP or Conda for maximum ease of use. Source code and documentation are available at [https://github.com/beardymcjohnface/Koverage](https://github.com/beardymcjohnface/Koverage).


# Statement of need

With the current state of sequencing technologies, it is trivial to generate terabytes of sequencing data for hundreds or even thousands of samples. Databases such as the Sequence Read Archive (SRA) and the European Nucleotide Archive (ENA), containing nearly 100 petabytes combined of sequencing data, are constantly being mined and reanalysed in bioinformatics analyses. Memory and I/O bottlenecks lead to under-utilisation of CPUs, and computational inefficiencies at such scales waste thousands of dollars in compute costs. I/O heavy processes in large parallel batches can result in significantly impaired performance. This is especially true for HPC clusters with shared file storage, or for cloud environments using cost-efficient bucket storage.

While there are existing tools for performing coverage calculations, they are not optimised for deployment at large scales, or when analysing large reference files. They require several complete I/O operations of the sequencing data in order to generate coverage statistics. Mapping to very large genomes requires large amounts of memory, or alternatively, aligning reads in chunks creating more I/O operations. Moving I/O operations into memory, for example via `tempfs` may alleviate I/O bottlenecks. However, this is highly system-dependent and will exacerbate memory bottlenecks. 

Koverage addresses these I/O bottlenecks by eliminating the sorting, reading, and writing of intermediate alignments. Koverage includes a kmer-based implementation to eliminate memory bottlenecks from screening large genomes. Koverage can be utilised as is, but has also been incorporated into the metagenomics pipelines Hecatomb [@hecatomb], Phables [@phables], and Reneo [@reneo].

# Implementation

Koverage is written in Snakemake [@snakemake] and Python, and uses the Snaketool [@snaketool] command line interface (CLI). Snaketool takes the user input command line arguments to build a runtime config file, generate the Snakemake command, and run the pipeline. Any unrecognised arguments are assumed to be Snakemake arguments and are added to the Snakemake command. For cluster- or cloud-based execution, users are encouraged to generate a Snakemake profile, and Koverage is compatible with Snakemake's Cookiecutter [@cookiecutter] template profiles. The only required inputs are the reference FASTA-format file (`--ref`), and the sample reads (`--reads`).

# Sample parsing
Koverage will parse reads (`--reads`) using MetaSnek `fastq_finder` [@metasnek]. Users supply either a directory of sequencing reads, or a tab-separated values (TSV) file of sample names and corresponding read filepaths. For a directory, sample names and read file pairs will be inferred from the file names. For a TSV file, sample names and filepaths are read from the file. More information and examples are available at [https://gist.github.com/beardymcjohnface/bb161ba04ae1042299f48a4849e917c8](https://gist.github.com/beardymcjohnface/bb161ba04ae1042299f48a4849e917c8)

# Mapping-based coverage

This is the default method for calculating coverage statistics. Reads are mapped sample-by-sample to the reference genome using Minimap2 [@minimap]. The alignments are piped to a script that collects the counts per contig and total counts per sample. Koverage uses mapping coordinates to collect read counts for _bins_ or _windows_ along each contig. This allows for a fast approximation of the coverage of each contig by at least one read (hitrate), and of the evenness of coverage (variance) for each contig. The final counts, mean, median, hitrate, and variance are written to a Python pickle. A second script calculates the Reads Per Million (RPM), Reads Per Kilobase Million (RPKM), Reads Per Kilobase (RPK), and Transcripts Per Million (TPM) like so:

__RPM__ = $\frac{10^6 \times N}{T}$

__RPKM__ = $\frac{ 10^6 \times N}{T \times L}$

__RPK__ = $\frac{N}{L}$

__TPM__ = $\frac{10^6 \times RPK}{R}$

Where:

 - N = number of reads mapped to the contig
 - T = Total number of mapped reads for that sample
 - L = length of contig in kilobases
 - R = sum of all RPK values for that sample

To generate fast estimations for mean, median, hitrate, and variance, Koverage first collects the counts of the start coordinates of mapped reads within _bins_ (or _windows_) across each contig (\autoref{fig:counts}). The user can customise bin width (default 100 bp); mean and median counts are comparable to read-depth when the binwidth is equal to the library's read length. Variance is calculated directly as the standard variance of the bin counts. The hitrate is calculated as the fraction of bins greater than zero.

![Windowed-coverage counts. Counts of start coordinates of mapped reads are collected for each bin across a contig. The counts array is used to calculate estimates for coverage hitrate and variance.\label{fig:counts}](fig1.png){ width=100% }

The coverage from all samples are collated, and a summary for each contig coverage by all samples is calculated. A summary HTML report is generated which includes interactive graphs and tables for both the per sample coverage, and the combined coverage from all samples. We utilized Datapane [@datapane] to embed a combined bar and line chart from Plotly [@plotly] and an interactive table displaying the results.

# Kmer-based coverage

Mapping to very large reference genomes can place considerable strain on computer resources. Koverage offers a kmer-based approach to estimating coverage. First, the reference genome is processed and kmers are sampled evenly across each contig. The user can customise kmer size, sampling interval, and minimum and maximum number of kmers to sample per contig. Jellyfish [@jellyfish] databases are created for each sample. The sampled reference kmers are queried against the sample kmer database. The kmer counts, and a kmer count array is created for each contig. The sum, mean, and median are calculated directly from the count array. Hitrate is calculated as the number of kmer counts > 0 divided by the total number of kmers queried. Variance is highly sensitive to large outliers, and kmer counts are especially prone to large outliers. Therefore, variance is calculated as the standard variance of the lowest 95 % of kmer counts.

# CoverM wrapper

Koverage includes a wrapper for the popular CoverM [@coverm] tool. CoverM can parse aligned and sorted reads in BAM format. It can also align reads with minimap2, saving the sorted alignments in a temporary filesystem (tempfs), and then process the aligned and sorted reads from tempfs. When a large enough tempfs is available, this method of running CoverM is extremely fast. However, if the tempfs is insufficient for storing the alignments, they are instead written to and read from regular disk storage which can be a significant I/O bottleneck. This wrapper in Koverage will generate alignments with Minimap2, sort and save them in BAM format with SAMtools [@samtools], and run CoverM on the resulting BAM file. CoverM is currently not available for macOS and as such, this wrapper will only run on Linux systems.

# Benchmarks

We tested Koverage's methods on the Pawsey Supercomputing Research Centre's Setonix HPC (commissioned in 2023) [@setonix] using a small coral metagenome dataset [@coral] consisting of 34 samples, a 360 Mbp metagenome assembly, and 9.1 GB of sequencing reads. This represents a typical metagenomics application in optimal conditions. Table 1 shows that CoverM is slightly faster than the default mapping-based method in spite of the extra read and write operations.

__Table 1: Coral metagenome benchmarks with high performance I/O__
 
| Method   | Runtime (HH:MM:SS) | CPU Walltime (HH:MM:SS) | Mean load (%) | Peak memory (Gb) |
|----------|--------------------|-------------------------|---------------|------------------|
| Map      | 00:40:34           | 01:49:38                | 270           | 4.6              |
| Kmer     | 02:20:58           | 00:51:40                | 37            | 4.2              |
| CoverM   | 00:31:49           | 01:12:17                | 227           | 7.4              |

We repeated the above benchmarking with Koverage directly reading and writing to Pawsey's S3 network bucket storage mounted using s3fs-fuse. Unlike the local scratch partition, this represents a scenario with a significant I/O bottleneck. Table 2 shows that while all methods are slower, Koverage's mapping and kmer methods perform much faster than the CoverM wrapper. The poor performance of the CoverM wrapper is entirely the result of generating the alignment BAM files, accounting for 85% of the overall runtime, rather than CoverM itself. 

__Table 2: Coral metagenome benchmarks with bottlenecked I/O__

| Method   | Runtime (HH:MM:SS) | CPU Walltime (HH:MM:SS) | Mean load (%) | Peak memory (Gb) |
|----------|--------------------|-------------------------|---------------|------------------|
| Map      | 03:34:15           | 01:49:01                | 50            | 4.6              |
| Kmer     | 03:18:33           | 01:13:53                | 14            | 4.6              |
| CoverM   | 11:13:39           | 01:32:10                | 22            | 7.3              |

# Acknowledgments

This work was supported by resources provided by the Pawsey Supercomputing Research Centre with funding from the Australian Government and the Government of Western Australia. The support provided by Flinders University for HPC research resources is acknowledged. This work was supported by an award from NIH NIDDK RC2DK116713 and an award from the Australian Research Council DP220102915. 

# References
