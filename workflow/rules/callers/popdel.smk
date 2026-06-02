rule clean_bam_for_popdel:
    input:
        bam=config["sample_path"] + "{sample}.bam",
    output:
        header=temp(config["output"] + "/callers/popdel/profiles/{sample}.header"),
        reheader=temp(config["output"] + "/callers/popdel/profiles/{sample}.reheader"),
        bam=config["output"] + "/callers/popdel/profiles/{sample}.reheader.bam",
        bai=config["output"] + "/callers/popdel/profiles/{sample}.reheader.bam.bai",
    params:
        script="workflow/scripts/sanitize_sq.awk",
    log:
        "logs/popdel-clean/{sample}.log",
    benchmark:
        "benchmarks/popdel-clean/{sample}.txt"
    threads: 2
    resources:
        mem_mb=16000,
        walltime_h=6,
    envmodules:
        "samtools/" + config["samtools_version"],
    shell:
        """

        (

            samtools view -H {input.bam} > {output.header}

            ./{params.script} {output.header} > {output.reheader}

            samtools reheader {output.reheader} {input.bam} > {output.bam}

            samtools index -@ {threads} {output.bam}

        ) > {log} 2>&1  

        """


rule popdel_interval:
    input:
        nonzero_bed=config["sample_path"] + "{sample}_nonzero.bed",
    output:
        interval=config["output"] + "/callers/popdel/profiles/{sample}.interval",
    log:
        profile=config["output"] + "/logs/popdel/{sample}-interval.log",
    resources:
        walltime_h=24,
        mem_mb=20000,
    benchmark:
        config["output"] + "/benchmarks/popdel-interval/{sample}.tsv"
    threads: 1
    shell:
        """

        awk '{print $1 ":" $2+1 "-" $3}' {input.nonzero_bed} > {output.interval}

        """


rule popdel_profile:
    input:
        ref=config["reference"],
        bam=rules.clean_bam_for_popdel.output.bam,
    output:
        profile=config["output"] + "/callers/popdel/profiles/{sample}.profile",
        profile_txt=config["output"] + "/callers/popdel/profiles/{sample}.profile.txt",
    log:
        profile=config["output"] + "/logs/popdel/{sample}-profile.log",
    resources:
        walltime_h=24,
        mem_mb=80000,
    benchmark:
        config["output"] + "/benchmarks/popdel-profile/{sample}.tsv"
    threads: 1
    envmodules:
        "HTSlib/" + config["htslib_version"],
        "popdel/" + config["popdel_version"],
        "bcftools/" + config["bcftools_version"],
    shell:
        """

        popdel profile {input.bam} -o {output.profile}  &> {log.profile}

        realpath {output.profile} > {output.profile_txt}

        """


rule popdel_filter:
    input:
        profile=config["output"] + "/callers/popdel/profiles/{sample}.profile",
    output:
        raw_profile=temp(
            config["output"] + "/callers/popdel/profiles/{sample}.raw.profile"
        ),
        profile=config["output"] + "/callers/popdel/profiles/{sample}.filtered.profile",
        profile_txt=config["output"]
        + "/callers/popdel/profiles/{sample}.filtered.profile.txt",
    log:
        profile=config["output"] + "/logs/popdel/{sample}-profile-filter.log",
    resources:
        walltime_h=24,
        mem_mb=80000,
    benchmark:
        config["output"] + "/benchmarks/popdel-profile/{sample}.tsv"
    threads: 1
    envmodules:
        "HTSlib/" + config["htslib_version"],
        "popdel/" + config["popdel_version"],
        "bcftools/" + config["bcftools_version"],
    shell:
        """

        popdel view {input.profile} -o {output.raw_profile} 

        grep -va 'HLA' {output.raw_profile} > {output.profile}

        realpath {output.profile} > {output.profile_txt}

        """


rule popdel_call:
    input:
        profile=rules.popdel_profile.output.profile,
        profile_txt=config["output"] + "/callers/popdel/profiles/{sample}.profile.txt",
    output:
        vcf=config["output"] + "/callers/popdel/{sample}.vcf",
    log:
        call=config["output"] + "/logs/popdel-call/{sample}-call.log",
    resources:
        walltime_h=24,
        mem_mb=80000,
    benchmark:
        config["output"] + "/benchmarks/popdel/{sample}.tsv"
    threads: 40
    envmodules:
        "HTSlib/" + config["htslib_version"],
        "popdel/" + config["popdel_version"],
        "bcftools/" + config["bcftools_version"],
    shell:
        """
        popdel call --out {output.vcf} {input.profile_txt} &> {log.call} 
        """
