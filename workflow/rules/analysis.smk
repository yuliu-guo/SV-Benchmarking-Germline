# generates fig 1a
rule summary_analysis:
    input:
        expand(
            rules.truvari.output.summary,
            type=config["truvari_runs"],
            sample=SAMPLES,
            caller=CALLERS,
            output=config["output"],
        ),
        script="workflow/r-scripts/01-truvari-summary.R",
    output:
        config["output"] + "/truvari/stat_summary.png",
    log:
        config["output"] + "/logs/truvari-summary.log",
    resources:
        mem_mb=8000,
        walltime_h=1,
    params: 
        results_dir = config["output"]
    envmodules:
        "R-bundle-Bioconductor/3.22-foss-2025a-R-4.5.1",
        "R-bundle-PREDICT/1.0.1-foss-2025a"
    shell:
        """
        Rscript {input.script} -i {params.results_dir} > {log} 2>&1
        """


# generates fig 1b (SVLEN) 
rule vcf_analysis:
    input:
        script="workflow/r-scripts/02-vcf-investigation.R",
        data=rules.summary_analysis.output,
    output:
        config["output"] + "/truvari/type-ignored/fig1_count_split_TP.Rds",
    log:
        config["output"] + "/logs/truvari-summary.log",
    params: 
        results_dir = config["output"]
    resources:
        mem_mb=16000,
        walltime_h=2,
    envmodules:
        "R-bundle-DALYCA/1.0.6-foss-2024a",
        "R-bundle-Bioconductor/3.22-foss-2025a-R-4.5.1",
        "R-bundle-PREDICT/1.0.1-foss-2025a"
    shell:
        """
        Rscript {input.script} -i {params.results_dir} > {log} 2>&1
        """

# RUN 03 OFF SERVER FOR GT COMPLILATION 

# RUN 04 MANUALLY - COMPARES TWO RUNS OF THIS PIPELINE

# 05 upsetr analysis 
rule upset_analysis:
    input:
        json=expand(
            rules.truvari_consistency.output.json,
            type="type-ignored",
            sample=SAMPLES,
            svtype=SVTYPES,
            output=config["output"],
        ),
        script="workflow/r-scripts/05-upsetr_const.R",
        data=rules.summary_analysis.output,
    output:
        config["output"] + "/truvari/type-ignored/_upset-degreesort.png",
    log:
        config["output"] + "/logs/truvari-upset.log",
    resources:
        mem_mb=8000,
        walltime_h=1,
    envmodules:
        "R-bundle-DALYCA/1.0.6-foss-2024a",
        "R-bundle-Bioconductor/3.22-foss-2025a-R-4.5.1",
        "R-bundle-PREDICT/1.0.1-foss-2025a", 
    params: 
        results_dir = config["output"]
    shell:
        """
        Rscript {input.script} > {log} 2>&1
        """

# 06 benchmark plots 
rule benchmark_analysis:
    input:
        expand(
            config["output"] + "/benchmarks/{caller}/{sample}.tsv",
            sample=SAMPLES,
            caller=CALLERS,
            output=config["output"],
        ),
        script="workflow/r-scripts/06-benchmarks.R",
    output:
        config["output"] + "/benchmarks/caller-summary.csv",
    log:
        config["output"] + "/logs/benchmark-summary.log",
    resources:
        mem_mb=8000,
        walltime_h=1,
    envmodules:
        "R-bundle-Bioconductor/3.22-foss-2025a-R-4.5.1",
        "R-bundle-PREDICT/1.0.1-foss-2025a"
    shell:
        """

        R -f {input.script} > {log} 2>&1

        """
