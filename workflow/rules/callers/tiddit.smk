rule tiddit:
    input:
        ref=config["reference"],
        ref_index=multiext(config["reference"], ".pac", ".fai"),
        bam=config["sample_path"] + "{sample}.bam",
    output:
        vcf=config["output"] + "/callers/tiddit/{sample}.vcf",
    params:
        tmp=directory(config["output"] + "/callers/tiddit/{sample}"),
    benchmark:
        config["output"] + "/benchmarks/tiddit/{sample}.tsv"
    log:
        config["output"] + "/logs/tiddit/{sample}.log",
    threads: 40
    resources:
        walltime_h=24,
        mem_mb=8000,
    shadow: "minimal"
    envmodules:
        "tiddit/" + config["tiddit_version"],
    shell:
        """
        tiddit --sv --bam {input.bam} -o {params.tmp} --ref {input.ref} --threads {threads} > {log} 2>&1
        """
