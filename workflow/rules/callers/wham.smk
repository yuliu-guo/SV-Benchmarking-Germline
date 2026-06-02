rule wham:
    input:
        bam=config["sample_path"] + "{sample}.bam",
        ref=config["reference"],
        ref_index=multiext(config["reference"], ".pac", ".fai"),
    output:
        vcf=config["output"] + "/callers/wham/{sample}.vcf",
    log:
        config["output"] + "/logs/wham/{sample}.log",
    threads: 40
    resources:
        mem_mb=80000,
        walltime_h=10,
    benchmark:
        config["output"] + "/benchmarks/wham/{sample}.tsv"
    envmodules:
        "wham/" + config["wham_version"],
    shell:
        """
        mkdir -p $(dirname {output.vcf})
        whamg -a {input.ref} -f {input.bam} > {output.vcf} 
        """
