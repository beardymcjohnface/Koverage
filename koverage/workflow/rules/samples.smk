"""
Function for parsing the 'Reads' config and identifying samples and read files
"""

def samplesFromDirectory(dir):
    """Parse samples from a directory"""
    outDict = {}
    samples, extensions = glob_wildcards(os.path.join(dir,"{sample}_R1{extensions}"))
    if not extensions:
        sys.stderr.write("\n"
                         "    FATAL: We could not parse the sequence file names from the specified directory.\n"
                         "    We are expecting {sample}_R1{extension}, and so your files should contain the \n"
                         "    characters '_R1' in the fwd reads and optional '_R2' in the rev reads. \n"
                         "    Alternatively you can specify a 2- or 3-column TSV file instead to declare the sample\n"
                         "    names and corresponding R1/R2 files. e.g. \n"
                         "    sample1\tpath/to/reads/sample1.1.fastq.gz\tpath/to/reads/sample1.2.fastq.gz\n"
                         "    sample2\tpath/to/reads/sample2.1.fastq.gz\tpath/to/reads/sample2.2.fastq.gz\n"
                         "    ..."
                         "    See https://hecatomb.readthedocs.io/en/latest/usage/#read-directory for more info\n"
                         "\n")
        sys.exit(1)
    else:
        for sample in samples:
            outDict[sample] = {}
            R1 = os.path.join(dir,f"{sample}_R1{extensions[0]}")
            R2 = os.path.join(dir,f"{sample}_R2{extensions[0]}")
            if os.path.isfile(R1):
                outDict[sample]["R1"] = R1
            else:
                sys.stderr.write(f"    FATAL: Error globbing files. {R1} does not exist.\n")
                sys.exit(1)
            if os.path.isfile(R2):
                outDict[sample]["R2"] = R2
            else:
                sys.stderr.write(f"    INFO: no R2 file found for {sample}\n")
                outDict[sample]["R2"] = str()
    return outDict


def samplesFromTsv(tsvFile):
    """Read samples and files from a TSV"""
    outDict = {}
    with open(tsvFile,"r") as tsv:
        for line in tsv:
            l = line.strip().split("\t")
            if len(l) >= 2:
                outDict[l[0]] = {}
                if os.path.isfile(l[1]):
                    outDict[l[0]]["R1"] = l[1]
                else:
                    sys.stderr.write(f"    FATAL: file {l[1]} for sample {l[0]} not found.")
                    sys.exit(1)
                if len(l) == 3 and l[2] != str():
                    if os.path.isfile(l[2]):
                        outDict[l[0]]["R2"] = l[2]
                    else:
                        sys.stderr.write(f"    FATAL: file {l[2]} for sample {l[0]} not found.")
                        sys.exit(1)
                else:
                    sys.stderr.write(f"    INFO: no R2 file specified for {l[0]}\n")
                    outDict[l[0]]["R2"] = str()
    return outDict


def parseSamples(readFileDir):
    """Parse samples from either a directory or TSV file"""
    if os.path.isdir(readFileDir):
        sampleDict = samplesFromDirectory(readFileDir)
    elif os.path.isfile(readFileDir):
        sampleDict = samplesFromTsv(readFileDir)
    else:
        sys.stderr.write("\n"
                         f"    FATAL: {readFileDir} is neither a file nor directory.\n"
                         "    See https://hecatomb.readthedocs.io/en/latest/usage/#read-directory for more info\n"
                         "\n")
        sys.exit(1)
    if len(sampleDict.keys()) == 0:
        sys.stderr.write("\n"
                         "    FATAL: We could not detect any samples at all.\n"
                         "    See https://hecatomb.readthedocs.io/en/latest/usage/#read-directory for more info\n"
                         "\n")
        sys.exit(1)
    return sampleDict


def writeSamplesTsv(dict, outfh):
    """Write the samples to a TSV file"""
    with open(outfh, "w") as out:
        for sample in dict.keys():
            out.write(f"{sample}\t{dict[sample]['R1']}\t{dict[sample]['R2']}\n")
    return None
