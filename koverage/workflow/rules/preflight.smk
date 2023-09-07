import glob
import os

from metasnek import fastq_finder, fasta_finder


# Concatenate Snakemake's own log file with the master log file
def copy_log_file():
    files = glob.glob(os.path.join(".snakemake", "log", "*.snakemake.log"))
    if not files:
        return None
    current_log = max(files, key=os.path.getmtime)
    shell("cat " + current_log + " >> " + config["args"]["log"])

onsuccess:
    copy_log_file()

onerror:
    copy_log_file()


# DIRECTORIES
dir = dict()
dir["base"] = os.path.join(workflow.basedir, "..")
dir["env"] = os.path.join(workflow.basedir, "envs")
dir["scripts"] = os.path.join(dir["base"], "scripts")

try:
    assert(config["args"]["output"]) is not None
    dir["out"] = config["args"]["output"]
except (KeyError, AssertionError):
    dir["out"] = "koverage.out"

dir["temp"] = os.path.join(dir["out"], "temp")
dir["log"] = os.path.join(dir["out"], "logs")
dir["paf"] = os.path.join(dir["out"], "pafs")
dir["hist"] = os.path.join(dir["out"], "histograms")
dir["result"] = os.path.join(dir["out"], "results")
dir["bench"] = os.path.join(dir["out"], "benchmarks")


config["refkmers"] = os.path.join(dir["temp"], os.path.basename(config["args"]["ref"]) + "." + str(config["args"]["kmer_size"]) + "mer.zst")
config["samplekmers"] = os.path.join(dir["result"], "sample_kmer_coverage." + str(config["args"]["kmer_size"]) + "mer.tsv.gz")
config["allkmers"] = os.path.join(dir["result"], "all_kmer_coverage." + str(config["args"]["kmer_size"]) + "mer.tsv.gz")


# PARSE SAMPLES
samples = dict()
samples["reads"] = fastq_finder.parse_samples_to_dictionary(config["args"]["reads"])
samples["names"] = list(samples["reads"].keys())


# PARSE REF(S)
references = fasta_finder.parse_fastas(config["args"]["ref"])
if len(references) > 1:
    config["args"]["ref"] = os.path.join(dir["temp"], "concatenated_refs.fasta")


# TARGETS
targets = dict()

if config["args"]["pafs"]:
    targets["pafs"] = expand(os.path.join(dir["paf"],"{sample}.paf.zst"), sample=samples["names"])
else:
    targets["pafs"] = []

targets["coverage"] = [
    os.path.join(dir["result"], "sample_coverage.tsv"),
    os.path.join(dir["result"], "all_coverage.tsv"),
]

if config["args"]["report"]:
    targets["coverage"].append(os.path.join(dir["result"], "report.html"))

targets["kmercov"] = [
    config["samplekmers"],
    config["allkmers"]
]

targets["coverm"] = [
    os.path.join(dir["result"], "sample_coverm_coverage.tsv")
]

targets["reports"] = [
    os.path.join(dir["out"], "koverage.samples.tsv")
]


# Add targets for pre-building the environments
targets["envs"] = []

for filename in os.listdir(dir["env"]):
    if filename.endswith(".yaml") or filename.endswith(".yml"):
        targets["envs"].append(os.path.join(dir["temp"], filename + ".done"))
