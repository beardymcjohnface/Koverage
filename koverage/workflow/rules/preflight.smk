import attrmap as ap


# DIRECTORIES
dir = ap.AttrMap()
dir.env = os.path.join(workflow.basedir, "envs")
dir.scripts = os.path.join(workflow.basedir, "scripts")

try:
    assert(ap.utils.to_dict(config.args)["output"]) is not None
    dir.out = config.args.output
except (KeyError, AssertionError):
    dir.out = "koverage.out"

dir.temp = os.path.join(dir.out, "temp")
dir.log = os.path.join(dir.out, "logs")
dir.bam = os.path.join(dir.out, "bams")
dir.hist = os.path.join(dir.out, "histograms")
dir.result = os.path.join(dir.out, "results")


# PARSE SAMPLES
include: config.modules[config.args.library]["preprocessing"]

samples = ap.AttrMap()
samples.reads = parseSamples(config.args.reads)
samples.names = list(ap.utils.get_keys(samples.reads))
samples = au.convert_state(samples, read_only=True)


# LIBRARY SPECIFIC RULES
include: config.modules[config.args.library]["mapping"]


# TARGETS
targets = ap.AttrMap()

if config.args.bams:
    targets.bams = expand(os.path.join(dir.bam,"{sample}.bam"), sample=samples.names)
else:
    targets.bams = []

if config.args.histograms:
    targets.bams = expand(os.path.join(dir.hist,"{sample}.png"), sample=samples.names)

targets.coverage = [
    os.path.join(dir.result, "sample_coverage.tsv"),
    os.path.join(dir.result, "all_coverage.tsv"),
    os.path.join(dir.result, "sample_summary.tsv"),
    os.path.join(dir.result, "all_summary.tsv")
]
