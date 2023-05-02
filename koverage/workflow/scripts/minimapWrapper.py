#!/usr/bin/env python3

import subprocess
import threading
import queue
from scipy.stats import kurtosis
import sys


### TODO: CAPTURE AND CHECK ERROR CODES FOR SYSTEM COMMANDS
### TODO: CAPTURE STDERR AND SAVE TO LOG FOR SYSTEM COMMANDS
### TODO: Other teadious crap that Snakemake usually handles?

# """For testing as a standalone script"""
# import attrmap as ap
# # test inputs
# snakemake.input.assembly = "ref.fa"
# snakemake.input.r1 = "test.r1.fastq.gz"
# snakemake.input.r2 = "test.r2.fastq.gz"
# snakemake.threads = "4"
# snakemake.params.bams = True
# # test outputs
# snakemake = ap.AttrMap()
# snakemake.params.max_depth = 300
# snakemake.output.hist = "test.hist"
# snakemake.output.kurt = "test.kurtosis"
# snakemake.output.bam = "test.bam"
# snakemake.output.count = "test.count"
# snakemake.output.lib = "test.lib"



def read_stools_sort_bam(pipe, q1, q2, q3):
    """Capture the samtools sort output and add it to three queues for samtools depth, read counts, and bam"""
    for line in iter(pipe.stdout.readline, b''):
        line = line.decode()
        q1.put(line)
        q2.put(line)
        q3.put(line)
    for q in [q1, q2, q3]:
        q.put(None)


def read_stools_sort(pipe, q1, q2):
    """Capture the samtools sort output and add it to two queues for samtools depth and read counts"""
    for line in iter(pipe.stdout.readline, b''):
        line = line.decode()
        q1.put(line)
        q2.put(line)
    for q in [q1, q2]:
        q.put(None)


"""Samtools depth processing functions"""
def initDepth():
    d = dict()
    for i in range(snakemake.params.max_depth):
        d[i] = 0
    return d


def dumpContig(ctg, d, fh):
    for depth in sorted(d.keys()):
        fh.write(f"{ctg}\t{depth}\t{str(d[depth])}\n")


def dumpKurtosis(ctg, l, fh):
    kurt = str(kurtosis(l))
    fh.write(f"{ctg}\t{kurt}\n")


def p3_read_depth(pipe):
    """Read from samtools depth and calc depth histogram and kurtosis for each contig"""
    outhist = open(snakemake.output.hist, 'w')
    outkurt = open(snakemake.output.kurt, 'w')
    curdepth = initDepth()
    curcontig = str()
    curkurt = list()
    for line in iter(pipe.stdout.readline, b''):
        line = line.decode()
        l = line.strip().split()
        try:
            l[1] = int(l[1])
            l[2] = int(l[2])
        except IndexError:
            sys.stderr.write(f"P3 line: {line}")
        if not curcontig == l[0]:
            if curcontig:
                dumpContig(curcontig, curdepth, outhist)
                dumpKurtosis(curcontig, curkurt, outkurt)
                curdepth = initDepth()
                curkurt = list()
            curcontig = l[0]
        curkurt.append(l[2])
        if l[2] > snakemake.params.max_depth:
            l[2] = snakemake.params.max_depth
        curdepth[l[2]] += 1
    dumpContig(curcontig, curdepth, outhist)
    dumpKurtosis(curcontig, curkurt, outkurt)
    outhist.close()
    outkurt.close()


def q1_stools_depth(q1):
    """Echo the samtools sort output to samtools depth and spawn depth parser (p3_read_depth)"""

    # start samtools depth
    depthcmd = ["samtools", "depth", "-aa", "-"]
    p3 = subprocess.Popen(depthcmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE)

    # Start depth parser thread
    thread_read_depth = threading.Thread(target=p3_read_depth, args=(p3,))
    thread_read_depth.start()

    # read from samtools sort and print to samtools depth
    while True:
        line = q1.get()
        if line is None:
            break
        p3.stdin.write(line.encode())
        p3.stdin.flush()

    # close stdin
    p3.stdin.close()

    # Wait for depth parser to finish
    thread_read_depth.join()

    # Close stdout
    p3.stdout.close()




def q3_save_bam(q3):
    """Save the bam file"""
    viewcmd = ["samtools", "view", "-b", "-h", "-o", snakemake.output.bam, "-"]
    pb = subprocess.Popen(viewcmd, stdin=subprocess.PIPE)
    while True:
        line = q3.get()
        # sys.stderr.write(f"Q3: {line}")
        if line is None:
            break
        pb.stdin.write(line.encode())
    pb.stdin.close()


def q2_read_counts(q2):
    """Parse the samtools sort output and get read counts and contig lens"""
    ctglen = dict()
    ctgcnt = dict()
    rcnt = 0
    while True:
        line = q2.get()
        # sys.stderr.write(f"Q2: {line}")
        if line is None:
            break
        l = line.strip().split()
        if line.startswith("@"):            # SAM header
            # sys.stderr.write(f"Starts with @\n")
            if line.startswith("@SQ"):
                # sys.stderr.write(f"Starts with @SQ\n")
                c = l[1].split(":")
                n = l[2].split(":")
                ctglen[c[1]] = n[1]
                ctgcnt[c[1]] = 0
        else:                               # SAM body
            rcnt += 1
            if not l[2] == "*":
                ctgcnt[l[2]] += 1
    with open(snakemake.output.ctgcnts, 'w') as outfh:
        for c in ctgcnt.keys():
            outfh.write(f"{c}\t{ctglen[c]}\t{ctgcnt[c]}\n")
    with open(snakemake.output.lib, 'w') as outfh:
        outfh.write(f"{str(rcnt)}\n")


# minimap2 command
mm2cmd = [
    "minimap2",
    "-t",
    str(snakemake.threads),
    "-ax",
    "sr",
    "--secondary=no",
    snakemake.input.assembly,
    snakemake.input.r1,
    snakemake.input.r2
]


# samtools sort command
sortcmd = ["samtools", "sort", "-@", str(snakemake.threads), "-O", "SAM"]


# Start minimap
p1 = subprocess.Popen(mm2cmd, stdout=subprocess.PIPE)


# Start samtools sort
p2 = subprocess.Popen(sortcmd, stdin=p1.stdout, stdout=subprocess.PIPE)


# Create queues for samtools depth and read counts from minimap|samtools sort
queue1 = queue.Queue()
queue2 = queue.Queue()


if snakemake.params.bams:
    queue3 = queue.Queue()
    # Capture samtools sort -> queue1-3
    thread_reader = threading.Thread(target=read_stools_sort_bam, args=(p2, queue1, queue2, queue3))
    thread_reader.start()
    # Read from q3 and save bam file
    thread_parser3 = threading.Thread(target=q3_save_bam, args=(queue3,))
    thread_parser3.start()
else:
    # Capture samtools sort -> queue1-2
    thread_reader = threading.Thread(target=read_stools_sort, args=(p2, queue1, queue2))
    thread_reader.start()


# Read from q1 to samtools depth and calc histograms
thread_parser1 = threading.Thread(target=q1_stools_depth, args=(queue1,))
thread_parser1.start()


# Read from q2 and get read counts
thread_parser2 = threading.Thread(target=q2_read_counts, args=(queue2,))
thread_parser2.start()


# wait for workers to finish
thread_parser1.join()
thread_parser2.join()

if snakemake.params.bams:
    thread_parser3.join()

# close pipe (probably not needed)
p2.stdout.close()
