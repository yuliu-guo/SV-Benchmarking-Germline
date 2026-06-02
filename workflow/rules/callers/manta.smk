
def get_scratch_directory(
    subdirs: List[str], force_workdir_scratch: bool = False
) -> str:
    """
    Auxiliary function that returns the scratch directory location
    for a given rule.

    Arguments
    ---------
    subdirs: List[str]
        List of subdirectory names.
    force_workdir_scratch: bool
        If true, force the scratch directory inside of the
        working directory.
    """
    basedir = Path(pipeline_workdir).joinpath("scratch")
    if subdirs:
        for element in subdirs:
            basedir = basedir.joinpath(element)
    return basedir.as_posix()
rule manta_calling:
    input:
        ref=config["reference"],
        ref_index=multiext(config["reference"], ".pac", ".fai"),
        cram=config["sample_path"] + "{sample}.cram",
        crai=config["sample_path"] + "{sample}.crai",
    output:
        vcf_diploid=config["output"] + "/callers/manta/{sample}.vcf.gz",
        tbi_diploid=config["output"] + "/callers/manta/{sample}.vcf.gz.tbi",
    log:
        config["output"] + "/logs/manta/{sample}.log",
    params:
        tmp1=get_scratch_directory(["manta_tmp1"], False),
        tmp2=get_scratch_directory(["manta_tmp2"], False),
        output=config["output"] + "/callers/manta/{sample}",
        sample="{sample}",
    threads: 40
    benchmark:
        config["output"] + "/benchmarks/manta/{sample}.tsv"
    envmodules:
        "Anaconda2/" + config["anaconda2_version"],
        "manta/" + config["manta_version"],
    resources:
        mem_gb=80,
        mem_mb=80000,
        walltime_h=20,
    shell:
        "( mkdir -p {params.tmp1}; "
        " ln -rs {input.cram} {params.tmp1}/{params.sample}.cram; "
        " ln -rs {input.crai} {params.tmp1}/{params.sample}.cram.crai; "
        " configManta.py "
        "     --bam {params.tmp1}/{params.sample}.cram "
        "     --referenceFasta {input.ref} "
        "     --runDir {params.tmp2}; "
        " {params.tmp2}/runWorkflow.py "
        "     --mode local "
        "     --jobs {threads} "
        "     --memGb {resources.mem_gb}; "
        " mkdir -p $(dirname {output.vcf_diploid}); "
        " mv {params.tmp2}/results/variants/diploidSV.vcf.gz {output.vcf_diploid}; "
        " mv {params.tmp2}/results/variants/diploidSV.vcf.gz.tbi {output.tbi_diploid}; "

        " rm -rf {params.tmp1} {params.tmp2}; "
        ") 2> {log}"
