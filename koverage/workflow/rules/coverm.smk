rule coverm_map_pe:
    input:
        ref = config["args"]["ref"],
        r1=lambda wildcards: samples["reads"][wildcards.sample]["R1"],
    output:
        bam = os.path.join(dir["temp"], "{sample}.bam"),
        bai = os.path.join(dir["temp"], "{sample}.bam.bai")
    params:
        r2 = lambda wildcards: samples["reads"][wildcards.sample]["R2"] if samples["reads"][wildcards.sample]["R2"] else "",
    threads:
        resources["med"]["cpu"]
    resources:
        mem_mb = resources["med"]["mem"],
        mem = str(resources["med"]["mem"]) + "MB",
        time = resources["med"]["time"]
    conda:
        os.path.join(dir["env"], "minimap.yaml")
    benchmark:
        os.path.join(dir["bench"], "coverm_map_pe.{sample}.txt")
    log:
        os.path.join(dir["log"], "coverm_map_pe.{sample}.err")
    shell:
        ("{{ "
        "minimap2 "
            "-t {threads} "
            "-ax sr "
            "--secondary=no "
            "{input.ref} "
            "{input.r1} "
            "{params.r2} "
        "| samtools sort "
            "-T {wildcards.sample} "
            "-@ {threads} - "
            "> {output.bam}; "
        "samtools index "
            "{output.bam}; "
        "}} 2> {log}")


rule coverm_bam2counts:
    input:
        os.path.join(dir["temp"], "{sample}.bam")
    output:
        os.path.join(dir["temp"], "{sample}.cov")
    params:
        params = config["params"]["coverm"]
    conda:
        os.path.join(dir["env"], "coverm.yaml")
    benchmark:
        os.path.join(dir["bench"],"coverm_bam2counts.{sample}.txt")
    log:
        os.path.join(dir["log"], "coverm_bam2counts.{sample}.err")
    shell:
        ("coverm contig "
            "-b {input} " 
            "{params.params} "
            "> {output} "
            "2> {log} ")


rule coverm_combine:
    input:
        expand(os.path.join(dir["temp"], "{sample}.cov"), sample=samples["names"])
    output:
        os.path.join(dir["result"], "sample_coverm_coverage.tsv")
    params:
        samples = samples["names"],
        dir = dir["temp"]
    run:
        with open(output[0], "w") as outfh:
            with open(input[0], "r") as infh:
                header = infh.readline().split("\t")
                header = [' '.join(element.split()[1:]) for element in header]
                header[0] = "Contig"
                header = '\t'.join(header)
            outfh.write("Sample\t" + header + "\n")
            for sample in params.samples:
                with open(os.path.join(params.dir, sample + ".cov"), "r") as infh:
                    infh.readline()
                    for line in infh:
                        outfh.write(sample + "\t" + line)


rule reneo_coverage:
    input:
        expand(os.path.join(dir["temp"],"{sample}.cov"),sample=samples["names"])
    output:
        os.path.join(dir["result"], "reneo.coverage.tsv")
    threads:
        resources["ram"]["cpu"]
    resources:
        mem_mb = resources["ram"]["mem"],
        mem = str(resources["ram"]["mem"]) + "MB",
        time = resources["ram"]["time"]
    benchmark:
        os.path.join(dir["bench"],"reneo_coverage.txt")
    log:
        os.path.join(dir["log"], "reneo_coverage.err")
    shell:
        ("for i in {input}; do "
            "tail -n+2 $i; "
        "done | "
            "awk -F '\t' '{{ sum[$1] += $2 }} END {{ for (key in sum) print key, sum[key] }}' "
            "> {output} 2> {log}; ")
