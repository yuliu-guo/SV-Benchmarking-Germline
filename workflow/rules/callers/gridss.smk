localrules: 
    gridss_rename 

rule gridss:
    input:
        ref=config["reference"],
        ref_index=multiext(config["reference"], ".pac", ".fai"),
        bam=config["sample_path"] + "{sample}.bam",
    output:
        vcf=config["output"] + "/callers/gridss/raw_{sample}.vcf",
    log:
        config["output"] + "/logs/gridss/{sample}.log",
    threads: 40
    resources:
        mem_mb=80*1024,
        mem_gb=80,
        walltime_h=10,
    benchmark:
        config["output"] + "/benchmarks/gridss/{sample}.tsv"
    shadow:
        "shallow"
    shell:
        """
                module load GRIDSS/2.13.2-foss-2024a-Java-17
                mkdir -p $(dirname {output.vcf})
                ( gridss -r {input.ref} -o {output.vcf} --jvmheap 40g -t {threads} {input.bam} ) > {log} 2>&1 
        """

rule gridss_annotate:
    input: 
        vcfs = expand(config["output"] + "/callers/gridss/raw_{sample}.vcf", sample=SAMPLES),
    output: 
        annotated_vcfs = expand(config["output"] + "/callers/gridss/raw_{sample}_annotated.vcf", sample=SAMPLES),
    threads: 1
    resources: 
        mem_mb = 20000,
        mem_gb = 20,
        walltime_h=10
    params:
        # script = os.path.abspath("workflow/scripts/annotate_gridss.R"),
        script = 'scripts/annotate_gridss.R',
        directory = os.path.abspath(config["output"] + "/callers/gridss/"),
        singularity = '/singularity/rbase_4.5.1-with-bio.sif'
    log:
        config["output"] + "/logs/gridss/annotate.log",
    shell:
        """
        singularity exec --bind {params.directory} --bind {params.script} -e {params.singularity} Rscript {params.script} {params.directory} > {log} 2>&1 
        """

rule gridss_rename: 
    input: 
        config["output"] + "/callers/gridss/raw_{sample}_annotated.vcf"
    output:
        config["output"] + "/callers/gridss/{sample}.vcf"
    shell: 
        """
        cp {input} {output}
        """
