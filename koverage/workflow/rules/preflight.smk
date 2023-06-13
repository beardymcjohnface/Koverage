import attrmap as ap
import attrmap.utils as au
import glob
import os

from metasnek import fastq_finder


# Concatenate Snakemake's own log file with the master log file
def copy_log_file():
    files = glob.glob(os.path.join(".snakemake", "log", "*.snakemake.log"))
    if not files:
        return None
    current_log = max(files, key=os.path.getmtime)
    shell("cat " + current_log + " >> " + config.args.log)

onsuccess:
    copy_log_file()

onerror:
    copy_log_file()


# DIRECTORIES
dir = ap.AttrMap()
dir.base = os.path.join(workflow.basedir, "..")
dir.env = os.path.join(workflow.basedir, "envs")
dir.scripts = os.path.join(dir.base, "scripts")

try:
    assert(ap.utils.to_dict(config.args)["output"]) is not None
    dir.out = config.args.output
except (KeyError, AssertionError):
    dir.out = "koverage.out"

dir.temp = os.path.join(dir.out, "temp")
dir.log = os.path.join(dir.out, "logs")
dir.paf = os.path.join(dir.out, "pafs")
dir.hist = os.path.join(dir.out, "histograms")
dir.result = os.path.join(dir.out, "results")
dir.bench = os.path.join(dir.out, "benchmarks")


config.refkmers = os.path.join(dir.temp, os.path.basename(config.args.ref) + "." + str(config.args.kmer_size) + "mer.zst")
config.samplekmers = os.path.join(dir.result, "sample_kmer_coverage." + str(config.args.kmer_size) + "mer.tsv.gz")
config.allkmers = os.path.join(dir.result, "all_kmer_coverage." + str(config.args.kmer_size) + "mer.tsv.gz")


# PARSE SAMPLES
samples = ap.AttrMap()
samples.reads = fastq_finder.parse_samples_to_dictionary(config.args.reads)
samples.names = list(ap.utils.get_keys(samples.reads))
samples = au.convert_state(samples, read_only=True)
fastq_finder.write_samples_tsv(samples.reads, os.path.join(dir.out, "samples.tsv"))

# TARGETS
targets = ap.AttrMap()

if config.args.pafs:
    targets.pafs = expand(os.path.join(dir.paf,"{sample}.paf.zst"), sample=samples.names)
else:
    targets.pafs = []

targets.coverage = [
    os.path.join(dir.result, "sample_coverage.tsv"),
    os.path.join(dir.result, "all_coverage.tsv"),
    # os.path.join(dir.result, "sample_summary.tsv"),
    # os.path.join(dir.result, "all_summary.tsv")
]

targets.kmercov = [
    config.samplekmers,
    config.allkmers
]

targets.benchmark = [
    os.path.join(dir.result, "sample_bench_coverage.tsv")
]

