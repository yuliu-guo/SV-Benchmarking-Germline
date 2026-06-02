rule delly:
    input:
        bam=config["sample_path"] + "{sample}.bam",
        ref_index=multiext(config["reference"], ".pac", ".fai"),
        ref=config["reference"],
    output:
        vcf=config["output"] + "/callers/delly/{sample}.vcf",
    resources:
        walltime_h=15,
        mem_mb=80000,
    benchmark:
        config["output"] + "/benchmarks/delly/{sample}.tsv"
    log:
        config["output"] + "/logs/delly/{sample}.log",
    threads: 40
    shell:
        """
        ml Delly/{config[delly_version]}
        export OMP_NUM_THREADS={threads}

        delly call -g {input.ref} -o {output.vcf} {input.bam} > {log} 2>&1

        """
