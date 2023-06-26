![](koverage.png)


[![](https://img.shields.io/static/v1?label=CLI&message=Snaketool&color=blueviolet)](https://github.com/beardymcjohnface/Snaketool)
[![](https://img.shields.io/static/v1?label=Licence&message=MIT&color=black)](https://opensource.org/license/mit/)
[![](https://img.shields.io/static/v1?label=Install%20with&message=PIP&color=success)](https://pypi.org/project/koverage/)
[![](https://github.com/beardymcjohnface/Koverage/actions/workflows/py-app.yaml/badge.svg)](https://github.com/beardymcjohnface/Koverage/actions/workflows/py-app.yaml/)
[![Documentation Status](https://readthedocs.org/projects/koverage/badge/?version=latest)](https://koverage.readthedocs.io/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/beardymcjohnface/Koverage/branch/main/graph/badge.svg?token=17P2ZEL44U)](https://codecov.io/gh/beardymcjohnface/Koverage)

---

Quickly get coverage statistics given reads and an assembly.

# Motivation

While there are tools that will calculate read-coverage statistics, they do not scale particularly well for large 
datasets, large sample numbers, or large reference FASTAs.
Koverage is designed to place minimal burden on I/O and RAM to allow for maximum scalability.

# Install

Koverage is still in development, but is available on PyPI.
Easy install: 

```shell
pip install koverage
```

Developer install:

```shell
git clone https://github.com/beardymcjohnface/Koverage.git
cd Koverage
pip install -e .
```

# Usage

Get coverage statistics from mapped reads (default method).

```shell
koverage run --reads readDir --ref assembly.fasta
```

Get coverage statistics using kmers (scales much better for very large reference FASTAs).

```shell
koverage run --reads readDir --ref assembly.fasta kmer
```

Any unrecognised commands are passed onto Snakemake.
Run Koverage on a HPC using a Snakemake profile.

```shell
koverage run --reads readDir --ref assembly.fasta --profile mySlurmProfile
```

# Test

You can test the methods like so.

```shell
# test default method
koverage test

# test all methods
koverage test map kmer bench
```

# Coverage methods

## Mapping-based (default)

```shell
koverage run ...
# or 
koverage run ... map
```

This method will map reads using minimap2 and use the mapping coordinates to calculate coverage.
This method is suitable for most applications.

## Kmer-based

```shell
koverage run ... kmer
```

This method calculates [Jellyfish](https://github.com/gmarcais/Jellyfish) databases of the sequencing reads.
It samples kmers from all reference contigs and queries them from the Jellyfish DBs to calculate coverage statistics.
This method is exceptionally fast for very large reference genomes.

## CoverM

```shell
koverage run ... bench
```

We've included a wrapper for [CoverM](https://github.com/wwood/CoverM) which you may find useful.
The wrapper manually runs minimap2 and then invokes CoverM on the sorted BAM file. 
It then combines the output from all samples like the other methods.
If you have a large tempfs/ you'll probably find it faster to run CoverM directly on your reads.

# Outputs

## Mapping-based

Default output files using fast estimations for mean, median, hitrate, and variance.

<details>
    <summary><b>sample_coverage.tsv</b></summary>
Per sample and per contig counts.

Column | description
--- | ---
Sample | Sample name derived from read file name
Contig | Contig ID from assembly FASTA
Count | Raw mapped read count
RPM | Reads per million
RPKM | Reads per kilobase million
RPK | Reads per kilobase
TPM | Transcripts per million
Mean | _Estimated_ mean read depth
Median | _Estimated_ median read depth
Hitrate | _Estimated_ fraction of contig with depth > 0
Variance | _Estimated_ read depth variance

</details>

<br>

<details>
    <summary><b>all_coverage.tsv</b></summary>
Per contig counts (all samples).

Column | description
--- | ---
Contig | Contig ID from assembly FASTA
Count | Raw mapped read count
RPM | Reads per million
RPKM | Reads per kilobase million
RPK | Reads per kilobase
TPM | Transcripts per million

</details>

_(more outputs to come, watch this space)_
    
## Kmer-based

Outputs for kmer-based coverage metrics.
Kmer outputs are gzipped as it is anticipated that this method will be used with very large reference FASTA files.

<details>
    <summary><b>sample_kmer_coverage.NNmer.tsv.gz</b></summary>
Per sample and contig kmer coverage.

Column | description
--- | ---
Sample | Sample name derived from read file name
Contig | Contig ID from assembly FASTA
Sum | Sum of sampled kmer depths
Mean | Mean sampled kmer depth
Median | Median sampled kmer depth
Hitrate | Fraction of kmers with depth > 0
Variance | Variance of lowest 95 % of sampled kmer depths

</details>

<br>

<details>
    <summary><b>all_kmer_coverage.NNmer.tsv.gz</b></summary>
Contig kmer coverage (all samples).

Column | description
--- | ---
Contig | Contig ID from assembly FASTA
Sum | Sum of sampled kmer depths
Mean | Mean sampled kmer depth
Median | Median sampled kmer depth

</details>

_(more outputs to come, watch this space)_    
