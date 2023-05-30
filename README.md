![](koverage.png)


[![](https://img.shields.io/static/v1?label=CLI&message=Snaketool&color=blueviolet&style=for-the-badge)](https://github.com/beardymcjohnface/Snaketool)
[![](https://img.shields.io/static/v1?label=Licence&message=MIT&color=black&style=for-the-badge)](https://opensource.org/license/mit/)
![](https://img.shields.io/static/v1?label=Install%20with&message=PIP&color=success&style=for-the-badge)

---

Quickly get coverage statistics given reads and an assembly.

# Motivation

While there are tools that will calculate read-coverage statistics, they do not scale particularly well for large 
datasets, large sample numbers, or large reference FASTAs.
Koverage is designed to place minimal burden on I/O and RAM to allow for maximum scalability.

# Install

Koverage is still in development and is not yet packaged on Bioconda or PyPI.
Setup.py may be missing some packages but just install them with pip.

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

# Outputs

## Mapping-based


<details>
    <summary><b>sample_coverage.tsv</b></summary>

Column | description
--- | ---
Sample | Sample name derived from read file name
Contig | Contig ID from assembly FASTA
Count | Raw mapped read count
RPM | Reads per million
RPKM | Reads per kilobase million
RPK | Reads per kilobase
TPM | Transcripts per million
Variance | _Estimated_ read depth variance

</details>

_(more outputs to come, watch this space)_
    
## Kmer-based

<details>
    <summary><b>sample_kmer_coverage.NNmer.tsv.gz</b></summary>

Column | description
--- | ---
Sample | Sample name derived from read file name
Contig | Contig ID from assembly FASTA
Mean | Mean sampled kmer depth
Median | Median sampled kmer depth
Hitrate | Fraction of kmers with depth > 0
Variance | Variance of sampled kmer depths

</details>

_(more outputs to come, watch this space)_    
