## 03. Truvari ##
SVTYPES = config["svtypes"]


rule split_svtypes:
    input:
        vcf=config["output"] + "/truvari/{type}/{caller}/{sample}/tp-base.vcf.gz",
        tbi=config["output"] + "/truvari/{type}/{caller}/{sample}/tp-base.vcf.gz.tbi",
    output:
        vcf=config["output"]
        + "/truvari-split/{type}/{caller}/{sample}/{svtype}.tp-base.vcf.gz",
    log:
        config["output"] + "/logs/split_svtypes/{type}_{caller}_{sample}_{svtype}.log",
    resources:
        mem_mb=2000,
        walltime_h=1,
    threads: 1
    envmodules:
        "bcftools/" + config["bcftools_version"],
    shell:
        """
        bcftools view -i 'SVTYPE=="{wildcards.svtype}"' {input.vcf} -Oz -o {output.vcf}
        bcftools index --tbi {output.vcf}
        """


rule truvari:
    input:
        test_vcf=config["output"] + "/callers/{caller}/{sample}.final.vcf.gz",
        sub_bed=config["truvari_bed"], 
        test_vcf_tbi=config["output"] + "/callers/{caller}/{sample}.final.vcf.gz.tbi",
        true_vcf=config["sample_path"] + "{sample}.vcf.gz",
        ref=config["reference"],
    output:
        fp=config["output"] + "/truvari/{type}/{caller}/{sample}/fp.vcf.gz",
        fn=config["output"] + "/truvari/{type}/{caller}/{sample}/fn.vcf.gz",
        tp_base=config["output"] + "/truvari/{type}/{caller}/{sample}/tp-base.vcf.gz",
        tp_base_tbi=config["output"]
        + "/truvari/{type}/{caller}/{sample}/tp-base.vcf.gz.tbi",
        summary=config["output"] + "/truvari/{type}/{caller}/{sample}/summary.json",
    log:
        config["output"] + "/logs/truvari/{type}/{caller}/{sample}.log",
    params:
        name="{type}",
        prefix=config["output"] + "/truvari/{type}/{caller}/{sample}",
        setting=lambda wildcards: config["truvari_runs"][wildcards.type],
    resources:
        mem_mb=8000,
        walltime_h=4,
    threads: 1
    envmodules:
        "truvari/" + config["truvari_version"],
    shell:
        """
        (mkdir -p {params.prefix};
        rm -rf {params.prefix}/temp;
        rm -rf {params.prefix}/phab_bench;
        truvari bench -b {input.true_vcf}  -f {input.ref}  --includebed {input.sub_bed} --extend 500 -c {input.test_vcf} -o {params.prefix}/temp {params.setting};
        mv --force {params.prefix}/temp/* {params.prefix}/.;) > {log} 2>&1
        """


# use only deletions on callers that only call deletions
use rule truvari as truvari_deletion with:
    input:
        test_vcf=config["output"] + "/callers/{caller}/{sample}.vcf.gz",
        test_vcf_tbi=config["output"] + "/callers/{caller}/{sample}.vcf.gz.tbi",
        sub_bed=config["truvari_bed"], 
        true_vcf=config["sample_path"] + "{sample}_del.vcf.gz",
        ref=config["reference"],
    wildcard_constraints:
        caller="popdel",


rule truvari_consistency:
    input:
        vcfs=expand(
            rules.split_svtypes.output,
            type="{type}",
            caller=CALLERS,
            sample="{sample}",
            svtype="{svtype}",
            output=config["output"],
        ),
    output:
        json=config["output"]
        + "/truvari-consistency/{type}/{sample}/{svtype}/consistency.json",
    log:
        config["output"] + "/logs/truvari_consistency/{type}_{sample}_{svtype}.log",
    threads: 1
    resources:
        mem_mb=4000,
        walltime_h=2,
    envmodules:
        "truvari/" + config["truvari_version"],
    shell:
        """
        truvari consistency -j {input.vcfs} > {output.json} 2> {log}
        """
