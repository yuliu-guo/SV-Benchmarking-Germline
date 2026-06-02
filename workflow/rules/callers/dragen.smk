# just sym link the results for easy reference
rule dragen:
    input:
        bam=config["sample_path"] + "{sample}.bam",
        sv_vcf=config["dragen_dir"] + "{sample}.sv.vcf.gz",
        sv_tbi=config["dragen_dir"] + "{sample}.sv.vcf.gz.tbi",
        cnv_vcf=config["dragen_dir"] + "{sample}.cnv.vcf.gz",
        cnv_tbi=config["dragen_dir"] + "{sample}.cnv.vcf.gz.tbi",
    output:
        sv_vcf=config["output"] + "/callers/dragen/{sample}.sv.vcf.gz",
        sv_tbi=config["output"] + "/callers/dragen/{sample}.sv.vcf.gz.tbi",
        cnv_vcf=config["output"] + "/callers/dragen/{sample}.cnv.vcf.gz",
        cnv_tbi=config["output"] + "/callers/dragen/{sample}.cnv.vcf.gz.tbi",
    params:
        dragen_dir=config["dragen_dir"],
    resources:
        walltime_h=1,
        mem_mb=1000,
    wildcard_constraints:
        sample="|".join(SAMPLES),
    threads: 1
    shell:
        """
                ln -s {params.dragen_dir}/{wildcards.sample}.sv.vcf.gz {output.sv_vcf}
                ln -s {params.dragen_dir}/{wildcards.sample}.sv.vcf.gz.tbi {output.sv_tbi}
                ln -s {params.dragen_dir}/{wildcards.sample}.cnv.vcf.gz {output.cnv_vcf}
                ln -s {params.dragen_dir}/{wildcards.sample}.cnv.vcf.gz.tbi {output.cnv_tbi}
                """


rule dragen_combine:
    input:
        sv_vcf=config["output"] + "/callers/dragen/{sample}.sv.vcf.gz",
        sv_tbi=config["output"] + "/callers/dragen/{sample}.sv.vcf.gz.tbi",
        cnv_vcf=config["output"] + "/callers/dragen/{sample}.cnv.vcf.gz",
        cnv_tbi=config["output"] + "/callers/dragen/{sample}.cnv.vcf.gz.tbi",
    output:
        vcf=config["output"] + "/callers/dragen/{sample}.vcf.gz",
        vcf_tbi=config["output"] + "/callers/dragen/{sample}.vcf.gz.tbi",
    wildcard_constraints:
        sample="|".join(SAMPLES),
    resources:
        walltime_h=2,
        mem_mb=16000,
    log:
        config["output"] + "/logs/dragen/{sample}-combine.log",
    threads: 1
    envmodules:
        "bcftools/" + config["bcftools_version"],
    shell:
        """
        bcftools concat {input.sv_vcf} {input.cnv_vcf} -o {output.vcf} -Oz --write-index="tbi" -aD > {log} 2>&1
        """
