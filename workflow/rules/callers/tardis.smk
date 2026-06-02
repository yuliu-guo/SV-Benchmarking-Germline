
rule tardis:
    input:
        ref=config["reference"],
        ref_index=multiext(config["reference"], ".pac", ".fai"),
        sonic_file=config["sonic"],
        bam=config["sample_path"] + "{sample}.bam",
    output:
        vcf=config["output"] + "/callers/tardis/{sample}.vcf" if exists(config["sonic"]) else [],
    params:
        out=config["output"] + "/callers/tardis/{sample}",
    resources:
        walltime_h=80,
        mem_mb=80000,
        mem_gb=80,
    benchmark:
        config["output"] + "/benchmarks/tardis/{sample}.tsv"
    threads: 40
    shadow:
        "shallow"
    log:
        config["output"] + "/logs/tardis/{sample}.log",
    envmodules:
        "TARDIS/" + config["tardis_version"],
    shell:
        """
        tardis -i {input.bam} --ref {input.ref} --sonic {input.sonic_file} --out {params.out}  > {log} 2>&1
        """
