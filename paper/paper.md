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
  - name: Brad Hart
    orcid: 0000-0001-8110-2460
    affiliation: 3
  - name: Sarah Beecroft
    orcid: 0000-0002-3935-2279
    affiliation: 4
  - name: Bhavya Papudeshi
    orcid: 0000-0001-5359-3100
    affiliation: 1
  - name: Laura Inglis
    orcid: 0000-0001-7919-8563
    affiliation: 1
  - name: Susanna Grigson
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
  - name: "Adelaide Centre for Epigenetics and the South Australian Immunogenomics Cancer Institute, Faculty of Health and Medical Sciences, The University of Adelaide, Adelaide, South Australia, Australia"
    index: 2
  - name: "Health and Biomedical Innovation, Clinical and Health Sciences, University of South Australia, SA, Australia"
    index: 3
  - name: "Pawsey Supercomputing Research Centre, Kensington, WA, Australia"
    index: 4
  - name: "Adelaide Medical School, Faculty of Health and Medical Sciences, The University of Adelaide, Adelaide, South Australia 5005, Australia"
    index: 5
  - name: "The Department of Surgery – Otolaryngology Head and Neck Surgery, Central Adelaide Local Health Network, Adelaide, South Australia 5000, Australia"
    index: 6
date:
bibliography: paper.bib
---

# Summary

Genomes of organisms are constructed by assembling short fragments (called sequencing reads) that are the resulting data outputs of whole genome sequencing (WGS). It is useful to determine the read-coverage of sequencing reads in the resulting genome assembly for many reasons, such as identifying duplication or deletion events, identifying related contigs for binning in metagenome assemblies [@metacoag;@graphbin2], or analysing taxonomic compositions of metagenomic samples [@condiga]. Although calculating the read-coverage of sequencing reads to a reference genome is a routine task, it typically involves several read and write operations (I/O operations) of the sequencing data. Although this is not a problem for small datasets, it can be a significant bottleneck when analysing a large number of samples, or when screening very large reference sequence files. Koverage is designed to reduce the I/O burden as much as possible to enable maximum scalability for large sample sizes. Koverage also includes a kmer-based coverage method that significantly reduces the computational complexity of screening large reference genomes such as the human genome. Koverage is a Snakemake [@snakemake] based pipeline, providing out-of-the-box support for HPC and cloud environments. It utilises the Snaketool [@snaketool] command line interface and is available to install via PIP or Conda for maximum ease-of-use. The source code and documentation is available at [https://github.com/beardymcjohnface/Koverage](https://github.com/beardymcjohnface/Koverage).


# Statement of need

With the current state of sequencing technologies, it is trivial to generate terabytes of sequencing data for hundreds or even thousands of samples. Databases such as the Sequence Read Archive (SRA) and the European Nucleotide Archive (ENA), containing nearly 100 petabytes combined of sequencing data, are constantly being mined and reanalysed in bioinformatics analyses. Computational inefficiencies at such scales waste thousands of dollars in compute costs and contribute to excess CO2 emissions. Memory and I/O bottlenecks can lead to under-utilisation of CPUs and exacerbate these inefficiencies. In severe cases, I/O heavy processes in large parallel batches can result in significantly impaired performance. This is especially true for HPC clusters with a shared scratch space of spinning disk hard drives, or for cloud-based analyses using cost-efficient network file systems or bucket storage.

While there are existing tools for performing coverage calculations, they are not optimised for deployment at large scales, or when analysing large reference files. This typically requires several complete I/O operations of the sequencing data in order to generate the coverage statistics. Furthermore, mapping to very large reference sequence files can require large amounts of memory, or alternatively, aligning reads in chunks and merging these chunked alignments at the end, resulting in even more I/O operations. Some proposed solutions involve moving I/O operations into memory, for example via `tempfs`. However, whether leveraging memory instead of I/O is a feasible option is highly system-dependent and will exacerbate any existing memory bottlenecks. 

Koverage addresses the I/O bottleneck of large datasets by eliminating the sorting, reading, and writing of intermediate alignment files. Koverage also includes a kmer-based implementation to eliminate memory bottlenecks that may arise from screening large reference files. Koverage can be utilised as is, but has also been incorporated into the metagenomics pipelines Hecatomb [@Hecatomb], Phables [@Phables], and Reneo [@Reneo].

# Implementation

Koverage is written in Snakemake [@snakemake] and Python, and uses the Snaketool [@snaketool] command line interface (CLI). The Snaketool CLI will take the user input command line arguments and Koverage's default configuration to build a runtime config file, build the Snakemake command and run the pipeline. Any unrecognised command line arguments are assumed to be Snakemake arguments and are added to the Snakemake command. For cluster- or cloud-based execution, users are encouraged to generate a Snakemake profile for their chosen deployment. Koverage has been designed to be compatible with Snakemake's Cookiecutter [@cookiecutter] template profiles. The only required inputs are the reference FASTA-format file (`--ref`), and the sample reads (`--reads`).

# Sample parsing
Koverage will parse reads (`--reads`) using MetaSnek `fastq_finder` [@metasnek]. Users supply either a directory containing their sequencing reads, or a tab-separated values (TSV) file listing their sample names and corresponding sequencing read filepaths. If users supply a directory to `--reads`, sample names and read file pairs will be inferred from the file names. If users supply a TSV file, sample names and filepaths will simply be read from the file. More information and examples are available at [https://gist.github.com/beardymcjohnface/bb161ba04ae1042299f48a4849e917c8](https://gist.github.com/beardymcjohnface/bb161ba04ae1042299f48a4849e917c8)

# Mapping-based coverage

This is the default method for calculating coverage statistics. Reads are mapped sample-by-sample to the reference genome using Minimap2 [@minimap]. The minimap2 alignments are parsed in real-time by a script that collects the counts per contig and total counts per sample. Koverage also uses the read mapping coordinates to collect read counts for `_bins_` or `_windows_` along the contig. This allows for a fast approximation of the coverage of each contig by at least one read (hitrate), and of the evenness of coverage (variance) for each contig. Following mapping, the final counts, mean, median, hitrate, and variance are written to a TSV file. A second script calculates the Reads Per Million (RPM), Reads Per Kilobase Million (RPKM), Reads Per Kilobase (RPK), and Transcripts Per Million (TPM) like so:

__RPM__ = $\frac{10^6 \times N}{T}$

__RPKM__ = $\frac{ 10^6 \times N}{T \times L}$

__RPK__ = $\frac{N}{L}$

__TPM__ = $\frac{10^6 \times RPK}{R}$

Where:

 - N = number of reads mapped to the contig
 - T = Total number of mapped reads for that sample
 - L = length of contig in kilobases
 - R = sum of all RPK values for that sample

As mentioned, Koverage uses a fast estimation for mean, median, hitrate, and variance. It estimates these values by first collecting the counts of the start coordinates of mapped reads within _bins_ (or _windows_) across each contig (\autoref{fig:counts}). The user can customise the bin width (default 100 bp); mean and median counts are comparable to read-depth when the binwidth is equal to the library's read length. Variance is calculated directly as the standard variance of the bin counts. The hitrate is calculated as the fraction of bins greater than zero.

![Windowed-coverage counts. Counts of start coordinates of mapped reads are collected for each bin across a contig. The counts array is used to calculate estimates for coverage hitrate and variance.\label{fig:counts}](fig1.png){ width=100% }

Lastly, the coverage from all samples are collated, and a summary of the coverage for each contig by all samples is calculated. A summary HTML report is then generated which includes interactive graphs and tables for both the per sample coverge, and the combined coverage from all samples. In the HTML report, we utilized Datapane [@datapane] to embed both a combined bar and line chart from Plotly [@plotly] and an interactive table displaying the results. This visualization represents the reads that have been mapped to each contig within the given reference sequence. The visualization is organized into two distinct tabs: one showcasing the individual read files with their associated mapping, and the other illustrating the combined read files with their respective mapping.

# Kmer-based coverage

Mapping to very large reference genomes can place considerable strain on computer resources. As an alternative, Koverage offers a kmer-based approach to estimating coverage across contigs. First, the reference genome is processed and kmers are sampled evenly across each contig. The user can customise the kmer size, sampling interval, and minimum and maximum number of kmers to sample for each contig. Jellyfish [@jellyfish] databases are then created for each sample. Koverage will initiate an interactive Jellyfish session for each sample's kmer database. The kmers that were sampled from each reference contig are queried against the sample kmer database and the kmer counts, and a kmer count array is created for each contig. The sum, mean, and median are calculated directly from the count array, and the hitrate is calculated as the number of kmer counts > 0 divided by the total number of kmers queried. As variance is highly sensitive to large outliers, and kmer counts are especially prone to large outliers for repetitive sequences, the variance is calculated as the standard variance of the lowest 95 % of kmer counts.

# CoverM wrapper

Koverage includes a wrapper for the popular CoverM [@coverm] tool. CoverM can parse aligned and sorted reads in BAM format. It can also align reads with minimap2, saving the sorted alignments in a temporary filesystem (tempfs), and then process the aligned and sorted reads from tempfs. When a large enough tempfs is available, this method of running CoverM is extremely fast. However, if the tempfs is insufficient for storing the alignments, they are instead written to and read from regular disk storage which can be a significant I/O bottleneck. This wrapper in Koverage will use Minimap2 to generate alignments, sort them and save them in BAM format with SamTools [@samtools], and then run CoverM on the resulting BAM file. While this is not the fastest method for running CoverM, it is convenient for users wishing to retain the sorted alignments in BAM format, and for automated running over many samples with a combined output summary file. CoverM is currently not available for MacOS and as such, this wrapper will only run on Linux systems.

# Benchmarks

We tested Koverage's methods on the Pawsey Supercomputing Research Centre's Setonix HPC (commissioned in 2023) [@setonix] using a small coral metagenome dataset [@coral] consisting of 68 samples, a 360 Mbp metagenome assembly, and 9.1 GB of sequencing reads. This represents a typical metagenomics application in optimal conditions. Table 1 shows that the CoverM wrapper is slightly faster than the default mapping-based method in spite of the extra read and write operations.

__Table 1: Coral metagenome benchmarks with high performance I/O__
 
| Method   | Runtime (HH:MM:SS) | CPU Walltime (HH:MM:SS) | Mean load (%) | Peak memory (Gb) |
|----------|--------------------|-------------------------|---------------|------------------|
| Map      | 00:40:34           | 01:49:38                | 270           | 4.6              |
| Kmer     | 02:20:58           | 00:51:40                | 37            | 4.2              |
| CoverM   | 00:31:49           | 01:12:17                | 227           | 7.4              |

We repeated the above benchmarking with Koverage directly reading and writing to Pawsey's S3 network bucket storage mounted using s3fs-fuse. Unlike Setonix's high performance local scratch partition, this represents a scenario with a significant I/O bottleneck. Table 2 shows that while all methods are slower, Koverage's mapping and kmer methods perform much faster than the CoverM wrapper. The poor performance of the CoverM wrapper is entirely the result of generating the alignment BAM files, which accounted for 85% of the overall runtime, rather than CoverM itself. 

__Table 2: Coral metagenome benchmarks with bottlenecked I/O__

| Method | Runtime (HH:MM:SS) | CPU Walltime (HH:MM:SS) | Mean load (%) | Peak memory (Gb) |
|--------|--------------------|-------------------------|---------------|------------------|
| Map    | 03:34:15           | 01:49:01                | 50            | 4.6              |
| Kmer   | 03:18:33           | 01:13:53                | 14            | 4.6              |
| CoverM | 11:13:39           | 01:32:10                | 22            | 7.3              |

# Acknowledgments

This work was supported by resources provided by the Pawsey Supercomputing Research Centre with funding from the Australian Government and the Government of Western Australia. The support provided by Flinders University for HPC research resources is acknowledged. This work was supported by an award from NIH NIDDK RC2DK116713 and an award from the Australian Research Council DP220102915. 

# References