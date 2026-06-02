
rule cnvpytor_calling:
    input:
        ref=config["reference"],
        ref_index=multiext(config["reference"], ".pac", ".fai"),
        cram=config["sample_path"] + "{sample}.cram",
    output:
        root=temp(config["output"] + "/callers/cnvpytor/{sample}.pytor"),
        tsv=temp(config["output"] + "/callers/cnvpytor/{sample}.tsv"),
    log:
        config["output"] + "/logs/cnvpytor/calling-{sample}.log",
    benchmark:
        config["output"] + "/benchmarks/cnvpytor-calling/{sample}.tsv"
    params:
        bin_size=100,
    threads: 40
    resources:
        mem_gb=80,
        mem_mb=80000,
        walltime_h=500,
    envmodules:
        "cnvpytor/" + config["cnvpytor_version"],
    shell:
        "( cnvpytor -root {output.root} -T {input.ref} -rd {input.cram}; "
        "echo === STEP 1 COMPLETE ; "
        " cnvpytor -root {output.root} -his {params.bin_size}; "
        "echo === STEP 2 COMPLETE ; "
        " cnvpytor -root {output.root} -partition {params.bin_size}; "
        "echo === STEP 3 COMPLETE ; "
        " cnvpytor -root {output.root} -call {params.bin_size} | sed -r 's/\tinf/\t9.9e+15/g' > {output.tsv};"
        "echo === STEP 4 COMPLETE ; "
        ") > {log} 2>&1"


rule cnvpytor_postprocess:
    input:
        ref_dict=config["human_ref_dict"],
        ref=config["reference"],
        ref_index=multiext(config["reference"], ".pac", ".fai"),
        tsv=rules.cnvpytor_calling.output.tsv,
    output:
        vcf_tmp1=temp(
            config["output"] + "/temp/cnvpytor/postprocess/{sample}.cnvpytor_tmp1.vcf"
        ),
        vcf_tmp2=temp(
            config["output"]
            + "/temp/cnvpytor/postprocess/{sample}.cnvpytor_tmp2.vcf.gz"
        ),
        tbi_tmp2=temp(
            config["output"]
            + "/temp/cnvpytor/postprocess/{sample}.cnvpytor_tmp2.vcf.gz.tbi"
        ),
        header=temp(
            config["output"] + "/temp/cnvpytor/postprocess/{sample}.cnvpytor.header"
        ),
        vcf=config["output"] + "/callers/cnvpytor/{sample}.vcf.gz",
        tbi=config["output"] + "/callers/cnvpytor/{sample}.vcf.gz.tbi",
    resources:
        mem_gb=20,
        mem_mb=20000,
        walltime_h=10,
    log:
        config["output"] + "/logs/cnvpytor/postprocess-{sample}.log",
    benchmark:
        config["output"] + "/benchmarks/cnvpytor-postprocess/{sample}.tsv"
    shell:
        "(module load BioPerl/{config[perl_version]}  cnvnator/{config[cnvnator_version]}  cnvpytor/{config[cnvpytor_version]}  HTSlib/{config[htslib_version]}  "
        "     bcftools/{config[bcftools_version]} Java/{config[java_version]}  "
        "     GATK/{config[gatk_version]} ; "


        " if [[ '{config[cnvnator_version]}' < '0.3.4' ]]; then "
        "     cnvnator2VCF.pl -prefix {wildcards.sample} -reference hg38 {input.tsv} | perl -F'[\\t]' -an -e '$F[3]=chr(78) if(!m/^#/); print join(chr(9),@F)'; "
        " else "
        "     cnvnator2VCF.pl -prefix {wildcards.sample} -reference hg38 {input.tsv}; "
        " fi | sed 's/##FORMAT=<ID=PE,Number=1,Type=Integer/##FORMAT=<ID=PE,Number=1,Type=Float/' > {output.vcf_tmp1}; "

        " gatk --java-options '-Xmx{resources.mem_gb}g' SortVcf "
        "      --INPUT {output.vcf_tmp1} "
        "      --SEQUENCE_DICTIONARY {input.ref_dict} "
        "      --COMPRESSION_LEVEL 9 "
        "      --CREATE_INDEX true "
        "      -O {output.vcf_tmp2}; "

        ' CALLCMD=$(echo "{rules.cnvpytor_calling.rule.shellcmd}" | tr -d "\'"); '
        " INFO='##SVCommandLine=<ID=cnvpytor,"
        'CommandLine="$CALLCMD",'
        'Version="{config[cnvpytor_version]}",'
        'Date="\'$(date +"%b %d, %Y %r %Z")\'">\'; '

        ' bcftools view -h {output.vcf_tmp2} | sed "\\$i $INFO" > {output.header}; '

        " bcftools reheader -s <(echo {wildcards.sample}) -h {output.header} {output.vcf_tmp2} | bcftools view -O z > {output.vcf}; "

        " tabix -p vcf {output.vcf}; "
        ") 2> {log}"
