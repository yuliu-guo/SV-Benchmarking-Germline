## Minigraph-Cactus 5-step pipeline Snakemake implementation
## Adapted from pipeline by Gregg Thomas, Harvard Informatics (https://github.com/harvardinformatics/cactus-snakemake)
## And from the official MC pangenome pipeline guide (https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md)
## By: Alban Obel (alosla), September 2025

import os
from pathlib import Path
from typing import List

def get_scratch_directory(subdirs: List[str], force_workdir_scratch: bool=False) -> str:
	'''
	Auxiliary function that returns the scratch directory location
	for a given rule.

	Arguments
	---------
	subdirs: List[str]
		List of subdirectory names.
	force_workdir_scratch: bool
		If true, force the scratch directory inside of the
		working directory.
	'''
	jobid = os.environ.get('PBS_JOBID', None)
	basedir = None
	if not jobid or force_workdir_scratch:
		basedir = Path(RESULTS_DIR).joinpath('scratch')
	else:
		basedir = Path('/scratch').joinpath(jobid)
	if subdirs:
		for element in subdirs:
			basedir = basedir.joinpath(element)
	return basedir.as_posix()


def retrieve_input(target_file: str) -> str:
	'''
	Check first for input in local scratch, then in results folder.
	'''
	location = None
	if os.path.isfile(os.path.join(SCRATCH_RESULTS_DIR, target_file)):
		location = os.path.join(SCRATCH_RESULTS_DIR, target_file)
	elif os.path.isfile(os.path.join(RESULTS_DIR, target_file)):
		location = os.path.join(SCRATCH_RESULTS_DIR, target_file)
	else:
		location = os.path.join(RESULTS_DIR, target_file)
	return location


PIPELINE_VERSION = "1.1"
PREFIX = config['prefix']
INPUT_SAMPLES = config['input_file']

RESULTS_DIR = os.path.join(config['output_dir'], f'{PREFIX}_MC-smk-{PIPELINE_VERSION}')
#JOIN_FINAL_RESULTS = os.path.join(RESULTS_DIR, '/final_results')

SCRATCH_BASE_DIR = get_scratch_directory([''], False)
SCRATCH_RESULTS_DIR = os.path.join(SCRATCH_BASE_DIR, f'{PREFIX}_MC-smk-{PIPELINE_VERSION}')
SCRATCH_INPUT_DIR = os.path.join(SCRATCH_RESULTS_DIR,'fasta_in')


RETRIEVED_RESULTS_DIR = get_scratch_directory(['retrieved_results'], False)

LOG_LEVEL = "INFO" ## CRITICAL ERROR WARNING INFO DEBUG TRACE

CHROMOSOMES = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM", "chrOther", "chrEBV"]

# Save outputs including intermediary results from following rules
ARCHIVE_RESULTS_RULES = ["minigraph_build", "graphmap", "graphmap_split", "align", "join"] 

#EXISTING_RESULTS = os.listdir(RESULTS_DIR)
EXISTING_RESULTS = [] 
file_list = [f'{PREFIX}.sv.gfa', f'{PREFIX}.paf'] + expand("chrom_alignments/{chrom}.{align_type}", chrom=CHROMOSOMES, align_type=['vg', 'hal'])
EXISTING_DIRS = []
dir_list = ['split_chroms', 'chrom_alignments', 'join_results']
for file in file_list:
	if os.path.isfile(os.path.join(RESULTS_DIR, file)):
		EXISTING_RESULTS.append(file)
for dir in dir_list:
	if os.path.isdir(os.path.join(RESULTS_DIR, dir)):
		EXISTING_DIRS.append(dir)

ruleorder: retrieve_results > minigraph_build > graphmap > graphmap_split > align > join

rule all:
	input:
		#os.path.join(RESULTS_DIR, f'{PREFIX}.sv.gfa'),
		#os.path.join(RESULTS_DIR, f'{PREFIX}.paf'),
		os.path.join(RESULTS_DIR, "join_results", f"{PREFIX}.gbz")
		

rule prepare_inputs:
	input:
		sample_list = os.path.join(config['input_files_dir'], INPUT_SAMPLES)
	output:
		local_sample_list = os.path.join(RESULTS_DIR, INPUT_SAMPLES),
		scratch_sample_list = os.path.join(RESULTS_DIR, f'scratchPaths.{INPUT_SAMPLES}'),
		temp_sample_list = os.path.join(SCRATCH_RESULTS_DIR, f'scratchPaths.{INPUT_SAMPLES}'),
		fasta_input_dir = directory(SCRATCH_INPUT_DIR)
	params:
		scratch_results = SCRATCH_RESULTS_DIR,
		input_samples = INPUT_SAMPLES
	threads: 10
	resources:
		mem_gb = 200,
		walltime_h = 1
	benchmark: RESULTS_DIR + '/benchmarks/prepare_inputs.tsv'
	log: RESULTS_DIR + '/logs/prepare_inputs.log'
	shell:
		'''
		(
		mkdir -p {output.fasta_input_dir}

		cp {input.sample_list} {output.local_sample_list}

		printf '' > {params.scratch_results}/{params.input_samples}
		while IFS=$'\t' read -r sample path
		do
			filename=$(basename $path)
			cp $path {output.fasta_input_dir}
			printf "$sample\t/scratch_input/$filename\n" >> {output.temp_sample_list}
		done < {output.local_sample_list}

		cp {output.temp_sample_list} {output.scratch_sample_list} 

		) > {log} 2>&1
		'''


## Takes previously generated intermediary files from the RESULTS_DIR and copies to scratch disk of the current pipeline job.
rule retrieve_results:
	input:
		os.path.join(RESULTS_DIR, "MINIGRAPH_BUILD_COMPLETE.flag")
	output:
		retrieved_results = expand(SCRATCH_RESULTS_DIR + f'/{{results}}', results=EXISTING_RESULTS),
		retrieved_dirs = directory(expand(SCRATCH_RESULTS_DIR + f'/{{dirs}}', dirs=EXISTING_DIRS))
	params:
		sv_gfa = os.path.join(RESULTS_DIR, f'{PREFIX}.sv.gfa'),
		scratch_sv_gfa = os.path.join(SCRATCH_RESULTS_DIR, f'{PREFIX}.sv.gfa'),
		scratch = SCRATCH_BASE_DIR,
		results_dir = RESULTS_DIR,
		scratch_results_dir = SCRATCH_RESULTS_DIR,
		existing_results = EXISTING_RESULTS,
		existing_dirs = EXISTING_DIRS
	threads: 10
	resources:
		mem_gb = 200,
		walltime_h = 1
	log: RESULTS_DIR + '/logs/retrieve_results.log'
	shell:
		'''
		(

		for file in {params.existing_results}; do
			cp -pr {params.results_dir}/$file {params.scratch_results_dir}
		done

		for dir in {params.existing_dirs}; do
			cp -pr {params.results_dir}/$dir {params.scratch_results_dir}
		done
		
		touch {params.results_dir}/RESULTS_RETRIEVED.flag

		) > {log} 2>&1
		'''



rule minigraph_build:
	input:
		scratch_sample_list = os.path.join(SCRATCH_RESULTS_DIR, f'scratchPaths.{INPUT_SAMPLES}')
	output:
		sv_gfa = os.path.join(RESULTS_DIR, f'{PREFIX}.sv.gfa'),
		scratch_sv_gfa = os.path.join(SCRATCH_RESULTS_DIR, f'{PREFIX}.sv.gfa')
		#scratch_sv_gfa = os.path.join(SCRATCH_RESULTS_DIR, f'{PREFIX}.sv.gfa')
	params:
		tmp = get_scratch_directory(['minigraph_tmp'], False),
		jobstore = get_scratch_directory(['minigraph_js'], False),
		scratch_input = SCRATCH_INPUT_DIR,
		scratch_results = SCRATCH_RESULTS_DIR,
		results_dir = RESULTS_DIR,
		singularity_version = config['singularity_version'],
		cactus_sif = config['cactus_2.9.9'],
		bind_string = RESULTS_DIR + "," + SCRATCH_INPUT_DIR + ":/scratch_input," + SCRATCH_BASE_DIR,
		file_prefix = PREFIX,
		reference = config['reference_samples'],
		log_level = LOG_LEVEL,
		configfile = config['cactus-minigraph-config']
	threads: 60
	resources:
		mem_gb = 500,
		walltime_h = 80
	benchmark: RESULTS_DIR + '/benchmarks/minigraph_build.tsv'
	log: RESULTS_DIR + '/logs/minigraph_build.log'
	shell:
		'''
		(
		ml load singularity/{params.singularity_version}

		rm -rf -v {params.jobstore}
		mkdir -p {params.tmp}

		singularity exec --bind {params.bind_string} -e -H {params.scratch_results} \
			{params.cactus_sif} cactus-minigraph \
			{params.jobstore} \
			{input.scratch_sample_list} \
			{output.scratch_sv_gfa} \
			--workDir {params.tmp} \
			--reference {params.reference} \
			--mgMemory {resources.mem_gb}G \
			--logLevel {params.log_level} \
			--configFile {params.configfile}

		cp {output.scratch_sv_gfa} {output.sv_gfa}

		touch {params.results_dir}/MINIGRAPH_BUILD_COMPLETE.flag

		rm -r {params.tmp}

		) > {log} 2>&1
		'''

rule graphmap:
	input:
		sv_gfa = os.path.join(SCRATCH_RESULTS_DIR, f'{PREFIX}.sv.gfa'),
		scratch_sample_list = os.path.join(SCRATCH_RESULTS_DIR, f'scratchPaths.{INPUT_SAMPLES}'),
		fasta_input_dir = directory(SCRATCH_INPUT_DIR)
	output:
		paf = os.path.join(RESULTS_DIR, f'{PREFIX}.paf'),
		scratch_paf = os.path.join(SCRATCH_RESULTS_DIR, f'{PREFIX}.paf'),
		fasta = os.path.join(RESULTS_DIR, f'{PREFIX}.sv.gfa.fa'),
		scratch_fasta = os.path.join(SCRATCH_RESULTS_DIR, f'{PREFIX}.sv.gfa.fa'),
		gaf = os.path.join(RESULTS_DIR, f'{PREFIX}.gaf.gz'),
		scratch_gaf = os.path.join(SCRATCH_RESULTS_DIR, f'{PREFIX}.gaf.gz'),
		paf_filter_log = os.path.join(RESULTS_DIR, f'{PREFIX}.paf.filter.log'),
		scratch_paf_filter_log = os.path.join(SCRATCH_RESULTS_DIR, f'{PREFIX}.paf.filter.log'),
		paf_unfiltered = os.path.join(RESULTS_DIR, f'{PREFIX}.paf.unfiltered.gz'),
		scratch_paf_unfiltered = os.path.join(SCRATCH_RESULTS_DIR, f'{PREFIX}.paf.unfiltered.gz')
	params:
		tmp = get_scratch_directory(['graphmap_tmp'], False),
		jobstore = get_scratch_directory(['graphmap_js'], False),
		scratch_input = SCRATCH_INPUT_DIR,
		scratch_results = SCRATCH_RESULTS_DIR,
		singularity_version = config['singularity_version'],
		cactus_sif = config['cactus_2.9.9'],
		bind_string = RESULTS_DIR + "," + SCRATCH_INPUT_DIR + ":/scratch_input," + SCRATCH_BASE_DIR,
		file_prefix = PREFIX,
		reference = config['reference_samples'],
		log_level = LOG_LEVEL,
		configfile = config['cactus-minigraph-config']
	threads: 60
	resources:
		mem_gb = 500,
		walltime_h = 40
	benchmark: RESULTS_DIR + '/benchmarks/graphmap.tsv'
	log: RESULTS_DIR + '/logs/graphmap.log'
	shell:
		'''
		(
		ml load singularity/{params.singularity_version}

		rm -rf -v {params.jobstore}
		mkdir -p {params.tmp}

		singularity exec --bind {params.bind_string} -e -H {params.scratch_results} \
			{params.cactus_sif} cactus-graphmap \
			{params.jobstore} \
			{input.scratch_sample_list} \
			{input.sv_gfa} \
			{output.scratch_paf} \
			--workDir {params.tmp} \
			--reference {params.reference} \
			--outputFasta {output.scratch_fasta} \
			--mapCores $(({threads}/4)) \
			--logLevel {params.log_level} \
			--configFile {params.configfile}
		## mapCores defaults to graphmap cpu = 6

		cp {output.scratch_paf} {output.paf}
		cp {output.scratch_fasta} {output.fasta}
		cp {output.scratch_gaf} {output.gaf}
		cp {output.scratch_paf_filter_log} {output.paf_filter_log}
		cp {output.scratch_paf_unfiltered} {output.paf_unfiltered}

		rm -r {params.tmp}

		) > {log} 2>&1
		'''



rule graphmap_split:
	input:
		sv_gfa = os.path.join(SCRATCH_RESULTS_DIR, f'{PREFIX}.sv.gfa'),
		paf = os.path.join(SCRATCH_RESULTS_DIR, f'{PREFIX}.paf'),
		scratch_sample_list = os.path.join(SCRATCH_RESULTS_DIR, f'scratchPaths.{INPUT_SAMPLES}'),
		fasta_input_dir = directory(SCRATCH_INPUT_DIR)
	output:
		chromfile = os.path.join(RESULTS_DIR, 'split_chroms/chromfile.sort.txt'),
		scratch_chromfile = os.path.join(SCRATCH_RESULTS_DIR, 'split_chroms/chromfile.sort.txt'),
		contig_sizes = os.path.join(RESULTS_DIR, "split_chroms/contig_sizes.tsv"),
		scratch_contig_sizes = os.path.join(SCRATCH_RESULTS_DIR, "split_chroms/contig_sizes.tsv"),
		minigraph_split_log = os.path.join(RESULTS_DIR, "split_chroms/minigraph.split.log"),
		scratch_minigraph_split_log = os.path.join(SCRATCH_RESULTS_DIR, "split_chroms/minigraph.split.log"),
		chrom_seqfiles = expand(os.path.join(SCRATCH_RESULTS_DIR, "split_chroms/seqfiles/{chrom}.seqfile"), chrom=CHROMOSOMES),
		chrom_pafs = expand(os.path.join(SCRATCH_RESULTS_DIR, "split_chroms/{chrom}/{chrom}.paf"), chrom=CHROMOSOMES)
	params:
		tmp = get_scratch_directory(['split_tmp'], False),
		jobstore = get_scratch_directory(['split_js'], False),
		scratch_results = SCRATCH_RESULTS_DIR,
		results_dir = RESULTS_DIR,
		singularity_version = config['singularity_version'],
		cactus_sif = config['cactus_2.9.9'],
		bind_string = RESULTS_DIR + "," + SCRATCH_INPUT_DIR + ":/scratch_input," + SCRATCH_BASE_DIR,
		reference = config['reference_samples'],
		log_level = LOG_LEVEL
	threads: 60
	resources:
		mem_gb = 500,
		walltime_h = 40
	benchmark: RESULTS_DIR + '/benchmarks/graphmap_split.tsv'
	log: RESULTS_DIR + '/logs/graphmap_split.log'
	shell:
		'''
		(
		ml load singularity/{params.singularity_version}

		rm -rf -v {params.jobstore}
		mkdir -p {params.tmp}

		singularity exec --bind {params.bind_string} -e -H {params.scratch_results} \
			{params.cactus_sif} cactus-graphmap-split \
			{params.jobstore} \
			{input.scratch_sample_list} \
			{input.sv_gfa} \
			{input.paf} \
			--outDir {params.scratch_results}/split_chroms \
			--workDir {params.tmp} \
			--reference {params.reference} \
			--logLevel {params.log_level}
		
		sed 's/chrX/chr23/;s/chrY/chr24/;s/chrM/chr25/;s/chrOther/chr26/;s/chrEBV/chr27/' < {params.scratch_results}/split_chroms/chromfile.txt |
		sort -V |
		sed 's/chr23/chrX/;s/chr24/chrY/;s/chr25/chrM/;s/chr26/chrOther/;s/chr27/chrEBV/' |
		sed 's#{params.scratch_results}##g' > {params.scratch_results}/split_chroms/chromfile.sort.txt

		for seqfile in {params.scratch_results}/split_chroms/seqfiles/*.seqfile
		do
			sed -i 's#file://{params.scratch_results}##g' $seqfile
		done

		cp -r {params.scratch_results}/split_chroms {params.results_dir}

		rm -r {params.tmp}

		) > {log} 2>&1
		'''



rule align:
	input:
		chrom_seqfile = os.path.join(SCRATCH_RESULTS_DIR, "split_chroms/seqfiles/{chrom}.seqfile"),
		chrom_paf = os.path.join(SCRATCH_RESULTS_DIR, "split_chroms/{chrom}/{chrom}.paf")
	output:
		chrom_hal = os.path.join(RESULTS_DIR, "chrom_alignments/{chrom}.hal"),
		scratch_chrom_hal = os.path.join(SCRATCH_RESULTS_DIR, "chrom_alignments/{chrom}.hal"),
		chrom_vg = os.path.join(RESULTS_DIR, "chrom_alignments/{chrom}.vg"),
		scratch_chrom_vg = os.path.join(SCRATCH_RESULTS_DIR, "chrom_alignments/{chrom}.vg")
	params:
		tmp = get_scratch_directory(['align_tmp_{chrom}'], False),
		jobstore = get_scratch_directory(['align_js_{chrom}'], False),
		alignment_output = os.path.join(SCRATCH_RESULTS_DIR, 'chrom_alignments'),
		scratch_results = SCRATCH_RESULTS_DIR,
		results_dir = RESULTS_DIR,
		singularity_version = config['singularity_version'],
		cactus_sif = config['cactus_2.9.9'],
		bind_string = RESULTS_DIR + "," + SCRATCH_INPUT_DIR + ":/scratch_input," + SCRATCH_RESULTS_DIR + "/split_chroms:/split_chroms," + SCRATCH_BASE_DIR,
		reference = config['reference_samples'],
		log_level = LOG_LEVEL,
		configfile = config['cactus-minigraph-config']
	threads: 10
	resources:
		mem_gb = 200,
		walltime_h = 24
	benchmark: RESULTS_DIR + '/benchmarks/align_{chrom}.tsv'
	log: RESULTS_DIR + '/logs/align_{chrom}.log'
	shell:
		'''
		(
		ml load singularity/{params.singularity_version}

		rm -rf -v {params.jobstore}
		mkdir -p {params.tmp}

		singularity exec --bind {params.bind_string} -e -H {params.scratch_results} \
			{params.cactus_sif} cactus-align \
			{params.jobstore} \
			{input.chrom_seqfile} \
			{input.chrom_paf} \
			{output.scratch_chrom_hal} \
			--pangenome \
			--outVG \
			--workDir {params.tmp} \
			--reference {params.reference} \
			--logLevel {params.log_level} \
			--configFile {params.configfile}
		
		cp {output.scratch_chrom_hal} {output.chrom_hal}
		cp {output.scratch_chrom_vg} {output.chrom_vg}

		rm -r {params.tmp}

		) > {log} 2>&1
		'''


rule join:
	input:
		chrom_hal_dir = expand(os.path.join(SCRATCH_RESULTS_DIR, "chrom_alignments/{chrom}.hal"), chrom=CHROMOSOMES),
		chrom_vg_dir = expand(os.path.join(SCRATCH_RESULTS_DIR, "chrom_alignments/{chrom}.vg"), chrom=CHROMOSOMES)
	output:
		final_gbz = os.path.join(RESULTS_DIR, "join_results", f"{PREFIX}.gbz"),
		#final_hal = os.path.join(RESULTS_DIR, "join_results", f"{PREFIX}.full.hal"),
		#final_gfa = os.path.join(RESULTS_DIR, "join_results", f"{PREFIX}.gfa.gz"),
		#final_vcf = os.path.join(RESULTS_DIR, "join_results", f"{PREFIX}.vcf.gz"),
		#final_vcf_index = os.path.join(RESULTS_DIR, "join_results", f"{PREFIX}.vcf.gz.tbi"),
		#final_dist = os.path.join(RESULTS_DIR, "join_results", f"{PREFIX}.dist"),
		#final_min = os.path.join(RESULTS_DIR, "join_results", f"{PREFIX}.min"),
		#final_raw_vcf = os.path.join(RESULTS_DIR, "join_results", f"{PREFIX}.raw.vcf.gz"),
		#final_raw_vcf_index = os.path.join(RESULTS_DIR, "join_results", f"{PREFIX}.raw.vcf.gz.tbi"),
		#final_stats = os.path.join(RESULTS_DIR, "join_results", f"{PREFIX}.stats.tgz")
	params:
		vg_input_dir = os.path.join(SCRATCH_RESULTS_DIR, 'chrom_alignments'),
		tmp = get_scratch_directory(['join_tmp'], False),
		jobstore = get_scratch_directory(['join_js'], False),
		scratch_join_out = os.path.join(SCRATCH_RESULTS_DIR, 'join_results'),
		scratch_results = SCRATCH_RESULTS_DIR,
		results_dir = RESULTS_DIR,
		file_prefix = PREFIX,
		chromosomes = CHROMOSOMES,
		singularity_version = config['singularity_version'],
		cactus_sif = config['cactus_2.9.9'],
		bind_string = RESULTS_DIR + "," + SCRATCH_BASE_DIR,
		reference = config['reference_samples'],
		log_level = LOG_LEVEL,
		configfile = config['cactus-minigraph-config']
	threads: 60
	resources:
		mem_gb = 500,
		walltime_h = 40
	benchmark: RESULTS_DIR + '/benchmarks/join.tsv'
	log: RESULTS_DIR + '/logs/join.log'
	shell:
		'''
		(
		ml load singularity/{params.singularity_version}

		rm -rf -v {params.jobstore}
		mkdir -p {params.tmp}
		mkdir -p {params.scratch_join_out}

		chroms=({params.chromosomes})
		hal_files_dir=()
		for chr in "${{chroms[@]}}"; do
			hal_files_dir+=("{params.vg_input_dir}/${{chr}}.hal")
		done

		vg_files_dir=()
		for chr in "${{chroms[@]}}"; do
			vg_files_dir+=("{params.vg_input_dir}/${{chr}}.vg")
		done

		## The sort order in which the --vg and --hal files are listed below should be considered
		## This order is preserved in reference paths of the graph and downstream in the contig sort order of surjected BAM files (and in VCFs for some variant callers)
		## Here the order is manually specified with the CHROMOSOMES variable
		singularity exec --bind {params.bind_string} -e -H {params.scratch_results} \
			{params.cactus_sif} cactus-graphmap-join \
			{params.jobstore} \
			--hal "${{hal_files_dir[@]}}" \
			--vg "${{vg_files_dir[@]}}" \
			--workDir {params.tmp} \
			--outDir {params.scratch_join_out} \
			--outName {params.file_prefix} \
			--reference {params.reference} \
			--logLevel {params.log_level} \
			--haplo clip \
			--gfa clip \
			--viz \
			--chrom-og \
			--configFile {params.configfile}
		
		cp -r {params.scratch_join_out} {params.results_dir}

		) > {log} 2>&1
		'''