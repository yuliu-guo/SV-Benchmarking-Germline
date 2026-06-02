##############
### LUMPY ###
##############

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

rule variants_cnv_discordant_reads:
    input:
        ref=config["reference"],
        cram=config["sample_path"] + "{sample}.cram",
    output:
        bam=config["sample_path"] + "{sample}.discordants.bam",
        bai=config["sample_path"] + "{sample}.discordants.bai",
    params:
        tmp1=get_scratch_directory(["discordant_read_tmp1"], False),
        tmp2=get_scratch_directory(["discordant_read_tmp2"], False),
    log:
        config["output"] + "/logs/discordant_reads/{sample}.log",
    benchmark:
        config["output"] + "/benchmarks/lumpy-discordant_reads/{sample}.tsv"
    threads: 2
    resources:
        mem_mb=80000,
        walltime_h=10,
    envmodules:
        "samtools/" + config["samtools_version"],
    shell:
        """
        (
        ml samtools/1.10

         echo "Creating temporary directories..."

         mkdir -p {params.tmp1}

         mkdir -p {params.tmp2}



         echo "Filtering and sorting discordant reads..."

         samtools view -T {input.ref} -h -u -F 1294 {input.cram} \
             | samtools sort -T {params.tmp2}/tmp -o {params.tmp1}/sorted.bam -



         echo "Indexing the sorted BAM file..."

         samtools index -@ {threads} {params.tmp1}/sorted.bam



         echo "Creating output directory..."

         mkdir -p $(dirname {output.bam})



         echo "Copying results..."

         cp {params.tmp1}/sorted.bam {output.bam}

         cp {params.tmp1}/sorted.bam.bai {output.bai}



         echo "Cleaning up temporary files..."

         rm -fr {params.tmp1}

         rm -fr {params.tmp2}

        ) > {log} 2>&1 

        """


rule variants_cnv_split_reads:
    input:
        ref=config["reference"],
        cram=config["sample_path"] + "{sample}.cram",
    output:
        bam=config["sample_path"] + "{sample}.splits.bam",
        bai=config["sample_path"] + "{sample}.splits.bai",
    log:
        config["output"] + "/logs/split_reads/{sample}.log",
    benchmark:
        config["output"] + "/benchmarks/lumpy-split_reads/{sample}.tsv"
    threads: 2
    params:
        tmp1=get_scratch_directory(["split_reads_tmp1"], False),
        tmp2=get_scratch_directory(["split_reads_tmp2"], False),
    resources:
        mem_mb=80000,
        walltime_h=20,
    envmodules:
        "Anaconda2/" + config["anaconda2_version"],
        "samblaster/" + config["samblaster_version"],
        "sambamba/" + config["sambamba_version"],
        "lumpy/" + config["lumpy_version"],
        "samtools/" + config["samtools_version"],
    shell:
        """
        echo "Starting split_reads for {wildcards.sample}" > {log}
        (
        ml samtools/1.10
        echo "Creating temporary directories... {params.tmp1} and {params.tmp2}"
        mkdir -p {params.tmp1}
        mkdir -p {params.tmp2}
        echo "Running samtools view..."
        samtools view -T {input.ref} -h {input.cram} \
            | extractSplitReads_BwaMem -i stdin \
            | samtools sort -T {params.tmp2}/tmp -o {params.tmp1}/sorted.bam -
        echo "Indexing BAM file..."
        samtools index -@ {threads} {params.tmp1}/sorted.bam {params.tmp1}/sorted.bai
        echo "Creating output directory..."
        mkdir -p $(dirname {output.bam})
        echo "Copying results..."
        cp {params.tmp1}/sorted.bam {output.bam}
        cp {params.tmp1}/sorted.bai {output.bai}
        echo "Cleaning up..."
        rm -fr {params.tmp1}
        rm -fr {params.tmp2}
        ) > {log} 2>&1 
        """


rule create_dict:
    input:
        config["reference"],
    output:
        config["human_ref_dict"],
    log:
        config["output"] + "/logs/create_dict.log",
    threads: 1
    resources:
        mem_mb=40000,
        walltime_h=2,
    envmodules:
        "picard-tools/" + config["picard-tools_version"],
        "htslib/" + config["htslib_version"],
        "bcftools/" + config["bcftools_version"],
        "java/" + config["java_version"],
        "gatk/" + config["gatk_version"],
    shell:
        """

        picard CreateSequenceDictionary \ 

                R={input} \ 

                O={output} > {log} 2>&1 

        """


lumpyexpress = os.path.join(workflow.current_basedir, "scripts", "lumpyexpress")


rule lumpy_calling:
    input:
        ref=config["reference"],
        ref_index=multiext(config["reference"], ".pac", ".fai"),
        cram=config["sample_path"] + "{sample}.cram",
        discordant_bam=rules.variants_cnv_discordant_reads.output.bam,
        splitter_bam=rules.variants_cnv_split_reads.output.bam,
    output:
        vcf=config["output"] + "/callers/lumpy/{sample}_tmp.vcf",
    log:
        config["output"] + "/logs/lumpy/calling-{sample}.log",
    benchmark:
        config["output"] + "/benchmarks/lumpy-calling/{sample}.tsv"
    params:
        tmpdir=get_scratch_directory(["lumpy_calling"], False),
    threads: 40
    group:
        "lumpy"
    resources:
        mem_mb=80000,
        walltime_h=16,
    envmodules:
        "Anaconda2/" + config["anaconda2_version"],
        "samblaster/" + config["samblaster_version"],
        "sambamba/" + config["sambamba_version"],
        "samtools/" + config["samtools_version"],
        "LUMPY/" + config["lumpy_version"],
    shell:
        # Using modified version of lumpyexpress (which does not fail when calling samtools view)
        " ({lumpyexpress} -v "
        "     -B {input.cram} "
        "     -S {input.splitter_bam} "
        "     -D {input.discordant_bam} "
        "     -R {input.ref} "
        "     -T {params.tmpdir} "
        "     -o {output.vcf}; "
        ") > {log} 2>&1 "


# TODO find svtyper-sso
rule lumpy_genotypes:
    input:
        cram=config["sample_path"] + "{sample}.cram",
        vcf=rules.lumpy_calling.output.vcf,
    output:
        vcf=config["output"] + "/callers/lumpy/{sample}_geno.vcf",
    log:
        config["output"] + "/logs/lumpy/genotypes-{sample}.log",
    benchmark:
        config["output"] + "/benchmarks/lumpy-genotypes/{sample}.tsv"
    threads: 40
    group:
        "lumpy"
    resources:
        mem_mb=80000,
        walltime_h=16,
    shell:
        # Call genotypes
        " svtyper-sso "
        "     --core {threads} "
        "     -i {input.vcf} "
        "     -B {input.cram} "
        "     > {output.vcf}; "
        ") > {log} 2>&1 "


rule lumpy_postprocess:
    input:
        ref_dict=config["human_ref_dict"],
        vcf=rules.lumpy_genotypes.output.vcf,
    output:
        sorted=config["output"] + "/callers/lumpy/{sample}.lumpy_sorted.vcf",
        header=config["output"] + "/callers/lumpy/{sample}.lumpy.header",
        vcf=config["output"] + "/callers/lumpy/{sample}.vcf",
    log:
        config["output"] + "/logs/lumpy/postprocess-{sample}.log",
    threads: 2
    group:
        "postprocess"
    resources:
        mem_gb=16,
        mem_mb=16000,
        walltime_h=4,
    envmodules:
        "HTSlib/" + config["htslib_version"],
        "bcftools/" + config["bcftools_version"],
        "Java/" + config["java_version"],
        "GATK/" + config["gatk_version"],
    shell:
        "(gatk --java-options '-Xmx{resources.mem_gb}g' SortVcf "
        "     --INPUT {input.vcf} "
        "     --SEQUENCE_DICTIONARY {input.ref_dict} "
        "     --COMPRESSION_LEVEL 0 "
        "     --QUIET true "
        "     --CREATE_INDEX false "
        "     --OUTPUT {output.sorted}; "

        " INFO='##SVCommandLine=<ID=lumpy,"
        'Version="{config[lumpy_version]}",'
        'Date="\'$(date +"%b %d, %Y %r %Z")\'">\'; '
        ' bcftools view -h {output.sorted} | sed "\\$i $INFO" > {output.header}; '
        " bcftools reheader -h {output.header} -o {output.vcf} {output.sorted}; "
        ") > {log} 2>&1 "
