rule smoove:
    input:
        ref=config["reference"],
        ref_index=multiext(config["reference"], ".pac", ".fai",".bwt"),
        cram=config["sample_path"] + "{sample}.cram",
        cram_index =config["sample_path"] + "{sample}.crai",
    output:
        temp_vcf=directory(config["output"] + "/callers/smoove/{sample}-smoove.vcf.gz"),
        vcf=config["output"] + "/callers/smoove/{sample}.vcf.gz",
    params:
        directory = subpath(output.vcf,parent=True),
        sample=subpath(output.vcf, basename=True, strip_suffix=".vcf.gz"),
    threads: 40
    resources:
        walltime_h=24,
        mem_mb=80000,
    benchmark:
        config["output"] + "/benchmarks/smoove/{sample}.tsv"
    log:
        config["output"] + "/logs/smoove/{sample}.log",
    singularity:
        "smoove_0.2.8-svtyper_0.7.1-python_2.7__v0.1__.sif"
 #   envmodules:
 #       "Anaconda2/" + config["anaconda2_version"],
 #       "samblaster/" + config["samblaster_version"],
 #       "sambamba/" + config["sambamba_version"],
 #       "bcftools/" + config["bcftools_version"],
 #       "samtools/" + config["samtools_version"],
 #       "gsort/" + config["gsort_version"],
 #       "LUMPY/" + config["lumpy_version"],
 #       "smoove/" + config["smoove_version"],
    shell:
        """
        (mkdir -p {params.directory}
        smoove call --name {params.sample} --genotype  --fasta {input.ref} -x {input.cram} -p {threads} 
        mv {params.sample}-smoove.vcf.gz {output.vcf}) > {log} 2>&1 
        """

