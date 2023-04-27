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

targets.


