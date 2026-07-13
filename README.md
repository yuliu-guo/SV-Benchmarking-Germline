# Germline SV Benchmarking Pipeline


A Snakemake workflow for running multiple structural‐variant (SV) callers on short-reads data and benchmarking their results against a single truth set with **[truvari](https://github.com/spiralgenetics/truvari)**.

An equivalent **[omnibenchmark](https://docs.omnibenchmark.org/latest/)** pipeline (with a subset of callers) can be found at: https://github.com/cphgeno/OB_SV_BENCHMARK_GERMLINE_SHORTREADS-main. This pipeline has been implemented with omnibenchmark v.0.3.2.  

---

##  Features

* Parallel execution on a compute cluster or locally.
* Modular caller rules—easy to add new SV callers.
* Automatic organization of outputs (VCFs, benchmark metrics, truvari stats).
* R-based downstream analysis.

---

##  Quick Start

```bash
snakemake --configfile config/pedigree.yaml --profile default
```

Or specifying the resources
```bash
snakemake --configfile config/pedigree.yaml --jobs 10 --cluster "qsub -W group_list=group_name -A group_name -l walltime={resources.walltime_h}:00:00 -l mem={resources.mem_mb}MB -l nodes=1:ppn={threads} -M email@address.com" --keep-going --groups truvari=trugroup --group-components trugroup=10

```

*Jobs are submitted to the cluster using the `default` profile.*

---

##  Project Layout

```
project/
├── config/
|   ├── graph.yaml
│   ├── tool_versions.yaml        
│   └── pedigree.yaml          
├── minigraph-cactus-snakemake-1.1/   # seperate snakemake pipeline for graph-aligned BAM files    
├── workflow/
│   ├── Snakefile            # Entry point
│   ├── profiles/default/    # Cluster profile 
│   ├── rules/
│   │   ├── callers/         # One .smk per SV caller
│   │   │   ├── cnvpytor.smk
│   │   │   ├── dragen.smk
│   │   │   ├── gridss.smk
│   │   │   ├── manta.smk
│   │   │   ├── popdel.smk
│   │   │   ├── svaba.smk
│   │   │   ├── tiddit.smk
│   │   │   ├── delly.smk
│   │   │   ├── dysgu.smk
│   │   │   ├── lumpy.smk
│   │   │   ├── octopus.smk
│   │   │   ├── smoove.smk
│   │   │   ├── tardis.smk
│   │   │   └── wham.smk
│   │   ├── truvari.smk
│   │   └── analysis.smk
│   ├── r-scripts/           # R scripts for analysis
│   └── scripts/             # Helper scripts
└── results/                 # Created at runtime

```

---

##  Configuration

Edit `config/config.yaml` to set:
Replace every /PATH/TO/... placeholder with the actual locations on your system.

* `sample_path`: Read-write directory containing the input BAMs (`[ID].bam`).
Symlinking your BAMs into this directory is recommended. 
* `output`: Directory where all pipeline results will be written.
* `CALLERS`: List of SV callers to run (names must match the rules in workflow/rules/callers/).
* `truth_set`: Full path to the truth VCF used for benchmarking.


---

##  Results

The pipeline creates:

* **benchmarks/** – Snakemake `--benchmark-extended` outputs
* **callers/** – VCFs and intermediate files for each caller/sample
* **truvari/** – Per-run comparison outputs and statistics

> Some R scripts may run better locally; see `workflow/r-scripts/` for details.

---

##  Modifying the Pipeline

### Adding a New Caller

1. Create `workflow/rules/callers/{caller}.smk`

   * Ensure a rule outputs:
     `config["output"]/callers/{caller}/{sample}.vcf.gz`
2. Add `{caller}` to the `CALLERS` list in `config/config.yaml`.

### Adding or Removing Samples

* Place new `[ID].bam` files in `sample_path`.
* Add or remove sample names in `config/config.yaml`.

### Changing the Truth Set

* Update `truth_set` in the YAML.
* Adjust the `samples` list and optionally set a new output directory.




