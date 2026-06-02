#!/bin/bash
#PBS -W group_list=group_name
#PBS -A group_name
#PBS -m ae
#PBS -M email@address.com
#PBS -l nodes=1:ppn=1
#PBS -l mem=5GB
#PBS -l walltime=96:00:00
#PBS -j oe
#PBS -o $PBS_JOBID.log

#ml purge
ml load snakemake/7.25.0-gfbf-2023a

cd /PATH/TO/minigraph-cactus-snakemake-1.1

snakemake --configfile minigraph-cactus-config-hprc_v1.0_pped.yaml \
        --snakefile minigraph-cactus.smk \
        --keep-going \
        --cores 80 \
        --latency-wait 30 \
        --printshellcmds \
        --rerun-triggers mtime \
        --rerun-incomplete #--dry-run
        
