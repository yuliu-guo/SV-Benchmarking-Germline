#!/usr/bin/env python3

import os
import subprocess
import json
import pandas as pd
import logging
from shutil import rmtree
from datetime import datetime
import glob
import pickle
import argparse
import sys



def setup_logging(log_dir=" "):
    """Setup logging with both file and console output."""
    os.makedirs(log_dir, exist_ok=True)
    log_filename = os.path.join(log_dir, f"sv_combination_pedigree_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
    
    # Create a logger
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    
    # Create file handler
    file_handler = logging.FileHandler(log_filename)
    file_handler.setLevel(logging.INFO)
    file_format = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    file_handler.setFormatter(file_format)
    
    # Create console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(file_format)
    
    # Add handlers to logger
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    
    return logger


def run_command(cmd, description="Running command", check=True):
    """Run a command with error handling and logging."""
    logger.info(f"{description}: {cmd}")
    try:
        result = subprocess.run(cmd, shell=True, check=check, 
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               text=True)
        return True, result
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed: {cmd}")
        logger.error(f"Error: {str(e)}")
        logger.error(f"Stdout: {e.stdout}")
        logger.error(f"Stderr: {e.stderr}")
        return False, e
    except Exception as e:
        logger.error(f"Exception while running command: {cmd}")
        logger.error(f"Error: {str(e)}")
        return False, e



def run_truvari(truth_vcf, input_vcf, reference, includebed, output_dir):
    """Run truvari bench and refine on a VCF file and return the results."""
    # Create output directory if it doesn't exist
    if os.path.exists(output_dir):
        logger.info(f"Removing existing output directory: {output_dir}")
        rmtree(output_dir)
    
    logger.info(f"Running truvari bench on {input_vcf}")

    bench_cmd = f"truvari bench -b {truth_vcf} -f {reference} --includebed {includebed} --extend 500 -c {input_vcf} -o {output_dir} --passonly --typeignore --pctseq 0"

    success, result = run_command(bench_cmd, "Running truvari bench")
    if not success:
        return None

    summary_file = f"{output_dir}/summary.json"
    if not os.path.exists(summary_file):
        logger.error(f"Summary file not found: {summary_file}")
        return None
    
    # Parse summary results
    try:
        with open(summary_file, "r") as f:
            summary = json.load(f)
        
        logger.info(f"Results for {input_vcf}: TP={summary['TP-base']}, FP={summary['FP']}, Recall={summary['recall']:.4f}, Precision={summary['precision']:.4f}, F1={summary['f1']:.4f}")
        
        return {
            "TP": summary["TP-base"],
            "FP": summary["FP"],
            "FN": summary["FN"],
            "recall": summary["recall"],
            "precision": summary["precision"],
            "f1": summary["f1"]
        }
    except Exception as e:
        logger.error(f"Error parsing summary file: {e}")
        return None



def get_merge_vcf_sets(toolcomb, vcf_files):
    """Generate the space-separated string of VCF files for the given tool combination."""
    to_return = []
    for tool in toolcomb:
        #to_return.append(f'{vcf_files[tool]}')
        # Uncomment the line below if you want to add tool names as tags
        to_return.append(f'{vcf_files[tool]}:{tool}')
    return ' '.join(to_return)



def merge_vcfs(tools_to_merge, vcf_files, output_vcf, svdb_sif_path):
    """Merge VCF files using SVDB with singularity container."""
    logger.info(f"Merging VCFs for tools: {tools_to_merge} to {output_vcf}")
    
    # Check if output file already exists
    if os.path.exists(output_vcf) and os.path.exists(f"{output_vcf}.tbi"):
        logger.info(f"Merged VCF already exists: {output_vcf}")
        return output_vcf
    
    # Get the VCF files string
    vcfs_to_merge = get_merge_vcf_sets(tools_to_merge, vcf_files)
    priority = ','.join(tools_to_merge)
    
    # Create bind string for singularity - bind all directories containing input/output files
    bind_dirs = set()
    
    # Add directories from VCF files
    for tool in tools_to_merge:
        if tool in vcf_files:
            bind_dirs.add(os.path.dirname(os.path.abspath(vcf_files[tool])))
    
    # Add output directory
    bind_dirs.add(os.path.dirname(os.path.abspath(output_vcf)))
    
    # Create bind string
    bind_string = ' '.join([f"-B {d}" for d in bind_dirs])
    
    # Build the merge command using singularity - pipe directly to bgzip
    merge_cmd = (f"singularity exec {bind_string} "
                f"{svdb_sif_path} svdb --merge --vcf {vcfs_to_merge} "
                f"--no_intra --no_var --overlap 0.6 --bnd_distance 1000 "
                f"--priority {priority} "
                f"| bcftools annotate -x ^FORMAT/GT "   # <-- fix header issue
                f"| bgzip > {output_vcf}")

    success, _ = run_command(merge_cmd, "Running SVDB merge with singularity")
    if not success:
        return None
    
    # Index the compressed VCF
    index_cmd = f"tabix -p vcf {output_vcf}"
    success, _ = run_command(index_cmd, "Indexing merged VCF")
    if not success:
        return None
    
    # Check if the output file exists
    if not os.path.exists(output_vcf):
        logger.error(f"Failed to create merged VCF: {output_vcf}")
        return None
    
    return output_vcf


def clean_temp_files(combined_dir, keep_file=None):
    """Clean up temporary files but keep the specified file."""
    temp_files = glob.glob(f"{combined_dir}/temp_*.vcf.gz*")
    for temp_file in temp_files:
        if keep_file and temp_file.startswith(keep_file):
            continue
        try:
            os.remove(temp_file)
            # Also remove the index file if it exists
            index_file = temp_file + ".tbi"
            if os.path.exists(index_file):
                os.remove(index_file)
            logger.info(f"Removed temporary file: {temp_file}")
        except Exception as e:
            logger.warning(f"Failed to remove temporary file {temp_file}: {e}")


def save_state(state, checkpoint_file):
    """Save the current state to a checkpoint file."""
    try:
        with open(checkpoint_file, 'wb') as f:
            pickle.dump(state, f)
        logger.info(f"State saved to {checkpoint_file}")
        return True
    except Exception as e:
        logger.error(f"Failed to save checkpoint: {e}")
        return False

def load_state(checkpoint_file):
    """Load the state from a checkpoint file."""
    try:
        with open(checkpoint_file, 'rb') as f:
            state = pickle.load(f)
        logger.info(f"State loaded from {checkpoint_file}")
        return state
    except Exception as e:
        logger.error(f"Failed to load checkpoint: {e}")
        return None



def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='SV caller combination evaluation with checkpointing')
    parser.add_argument('--resume', action='store_true', help='Resume from checkpoint')
    parser.add_argument('--checkpoint', default='', help='Checkpoint file path')
    parser.add_argument('--vcf-paths', default='', help='Path to VCF files')
    parser.add_argument('--vcf-paths-octopus', default='', help='Path to Octopus VCF files')
    parser.add_argument('--truthset', default='', help='Path to truth set VCF')
    parser.add_argument('--reference', default='', help='Path to reference genome')
    parser.add_argument('--includebed', default='', help='Path to include BED file')
    parser.add_argument('--results-dir', default='',help='Path to results directory')
    parser.add_argument('--combined-dir', default='',help='Path to combined results directory')
    parser.add_argument('--sample', default='NA12878.final',help='Sample name')
    parser.add_argument('--metric', default='recall', choices=['recall', 'precision', 'f1'], help='Metric to use for selecting the best tool (recall, precision, or f1)')
    parser.add_argument('--svdb-sif', default='', help='Path to SVDB singularity container')
    args = parser.parse_args()

    # Setup directories
    os.makedirs(args.results_dir, exist_ok=True)
    os.makedirs(args.combined_dir, exist_ok=True)

    singularity_tmp = ""
    os.makedirs(singularity_tmp, exist_ok=True)

    os.environ['SINGULARITY_TMPDIR'] = singularity_tmp
    os.environ['APPTAINER_TMPDIR'] = singularity_tmp
    os.environ['TMPDIR'] = singularity_tmp
    
    logger.info(f"Set SINGULARITY_TMPDIR to: {singularity_tmp}")
    # Initialize or load state
    if args.resume and os.path.exists(args.checkpoint):
        state = load_state(args.checkpoint)
        if state is None:
            logger.error("Failed to load checkpoint. Exiting.")
            sys.exit(1)
        
        # Extract state variables
        combination_set = state['combination_set']
        tools_pool = state['tools_pool']
        iteration = state['iteration']
        iterations = state['iterations']
        current_combined_vcf = state['current_combined_vcf']
        current_combined_name = state['current_combined_name']
        vcf_files = state['vcf_files']
        metric = state.get('metric', 'recall')  # Default to recall if not found in saved state
        
        # Override metric if specified in command line
        if args.metric and args.metric != metric:
            logger.warning(f"Overriding saved metric '{metric}' with command line metric '{args.metric}'")
            metric = args.metric
        
        logger.info(f"Resuming from iteration {iteration}")
        logger.info(f"Current combination set: {combination_set}")
        logger.info(f"Remaining tools in pool: {tools_pool}")
        logger.info(f"Using metric: {metric}")
    else:
        # Configuration
        sample = args.sample
        metric = args.metric
        
        logger.info(f"Using metric: {metric} for tool selection")
        
        # Original VCF files
        original_vcf_files = {
            "manta": f"{args.vcf_paths}/manta/{sample}.vcf.gz",
            "dysgu": f"{args.vcf_paths}/dysgu/{sample}.vcf.gz",
            "svaba": f"{args.vcf_paths}/svaba/{sample}.vcf.gz",
            "popdel": f"{args.vcf_paths}/popdel/{sample}.vcf.gz",
            "octopus": f"{args.vcf_paths_octopus}/{sample}.filtered.vcf.gz",
            "tiddit": f"{args.vcf_paths}/tiddit/{sample}.vcf.gz",
            "delly": f"{args.vcf_paths}/delly/{sample}.vcf.gz",
            "lumpy": f"{args.vcf_paths}/lumpy/{sample}.vcf.gz",
            "smoove": f"{args.vcf_paths}/smoove/{sample}.vcf.gz",
            "tardis": f"{args.vcf_paths}/tardis/{sample}.vcf.gz",
            "cnvpytor": f"{args.vcf_paths}/cnvpytor/{sample}.vcf.gz",
            "gridss": f"{args.vcf_paths}/gridss/{sample}.vcf.gz",
            "wham": f"{args.vcf_paths}/wham/{sample}.vcf.gz"
        }


        # Check if all files exist
        all_exist = True
        for tool, path in original_vcf_files.items():
            if not os.path.exists(path):
                logger.error(f"VCF file not found: {path}")
                all_exist = False
        
        if not all_exist:
            logger.error("Some VCF files are missing. Please check paths and try again.")
            sys.exit(1)
        
        # Initialize state variables
        combination_set = []
        tools_pool = list(original_vcf_files.keys())
        iteration = 0
        iterations = []
        current_combined_vcf = None
        current_combined_name = ""
        vcf_files = original_vcf_files.copy()
    
    # Main loop
    try:
        while tools_pool:
            iteration += 1
            logger.info(f"--- Iteration {iteration} ---")
            logger.info(f"Tools in pool: {tools_pool}")
            logger.info(f"Current combination set: {combination_set}")
            
            # Step 1: Evaluate each remaining tool in the pool
            tool_results = {}
            temp_merged_files = {}  # Store the paths to temp merged files for each tool
            all_evaluations_successful = True
            
            for tool in tools_pool:
                logger.info(f"Evaluating {tool}...")
                output_dir = f"{args.results_dir}/{tool}_iter{iteration}"
                
                # For first iteration, evaluate individual tools
                if not combination_set:
                    result = run_truvari(args.truthset, vcf_files[tool], args.reference, args.includebed, output_dir)
                    if result is None:
                        logger.error(f"Evaluation failed for tool: {tool}")
                        all_evaluations_successful = False
                        continue
                    temp_merged_files[tool] = vcf_files[tool]  # Just use the original VCF
                else:
                    # Create a name for the potential new combination (sorted alphabetically)
                    temp_combined_tools = sorted(combination_set + [tool])
                    temp_combined_name = "_".join(temp_combined_tools)
                    
                    # Prepare the merging - we'll use the current combined VCF and the new tool
                    temp_merged_vcf = f"{args.combined_dir}/temp_{temp_combined_name}.vcf.gz"
                    
                    # Use current combined VCF and the new tool
                    tools_to_merge = [tool, current_combined_name]
                    merged_vcf = merge_vcfs(tools_to_merge, vcf_files, temp_merged_vcf, args.svdb_sif)
                    if merged_vcf is None:
                        logger.error(f"Merging failed for tool combination: {tools_to_merge}")
                        all_evaluations_successful = False
                        continue
                    
                    result = run_truvari(args.truthset, temp_merged_vcf, args.reference, args.includebed, output_dir)
                    if result is None:
                        logger.error(f"Evaluation failed for merged VCF: {temp_merged_vcf}")
                        all_evaluations_successful = False
                        continue
                    
                    # Store the temp merged file path for this tool
                    temp_merged_files[tool] = temp_merged_vcf
                
                tool_results[tool] = result
                logger.info(f"  {tool} results: Recall={result['recall']:.4f}, Precision={result['precision']:.4f}, F1={result['f1']:.4f}")
            
            # Check if we have at least one successful evaluation
            if not tool_results:
                logger.error("All evaluations failed. Please check the logs and fix the issues.")
                logger.info("Saving state before exiting...")
                state = {
                    'combination_set': combination_set,
                    'tools_pool': tools_pool,
                    'iteration': iteration - 1,  # We're retrying the current iteration
                    'iterations': iterations,
                    'current_combined_vcf': current_combined_vcf,
                    'current_combined_name': current_combined_name,
                    'vcf_files': vcf_files,
                    'metric': metric
                }
                save_state(state, args.checkpoint)
                sys.exit(1)
            
            # Step 2: Find the tool with highest metric value (recall, precision, or f1)
            best_tool = max(tool_results.items(), key=lambda x: x[1][metric])
            best_tool_name = best_tool[0]
            best_tool_result = best_tool[1]
            
            logger.info(f"Best tool: {best_tool_name} with {metric} {best_tool_result[metric]:.4f}")
            
            # Add best tool to combination set
            combination_set.append(best_tool_name)
            
            # Step 3: Create the final merged VCF by reusing or renaming the temp merged file
            if len(combination_set) > 1:
                # Sort tools alphabetically for consistent naming
                sorted_tools = sorted(combination_set)
                combined_name = "_".join(sorted_tools)
                final_merged_vcf = f"{args.combined_dir}/{combined_name}.vcf.gz"
                
                # Check if we need to rename the temp file or if it's already correctly named
                temp_merged_vcf = temp_merged_files[best_tool_name]
                
                if temp_merged_vcf != final_merged_vcf:
                    # Copy/rename the temp merged file instead of merging again
                    copy_cmd = f"cp {temp_merged_vcf} {final_merged_vcf}"
                    success, _ = run_command(copy_cmd, "Copying temp merged file to final")
                    if not success:
                        logger.error(f"Failed to copy temp merged file to final: {temp_merged_vcf} -> {final_merged_vcf}")
                        # Save checkpoint before exiting
                        state = {
                            'combination_set': combination_set,
                            'tools_pool': tools_pool,
                            'iteration': iteration - 1,  # We're retrying the current iteration
                            'iterations': iterations,
                            'current_combined_vcf': current_combined_vcf,
                            'current_combined_name': current_combined_name,
                            'vcf_files': vcf_files,
                            'metric': metric
                        }
                        save_state(state, args.checkpoint)
                        sys.exit(1)
                    
                    # Copy the index file as well
                    index_copy_cmd = f"cp {temp_merged_vcf}.tbi {final_merged_vcf}.tbi"
                    success, _ = run_command(index_copy_cmd, "Copying temp merged index file to final")
                    if not success:
                        logger.warning(f"Failed to copy temp merged index file. Will recreate it.")
                        index_cmd = f"tabix -p vcf {final_merged_vcf}"
                        success, _ = run_command(index_cmd, "Indexing final merged VCF")
                        if not success:
                            logger.error(f"Failed to index final merged VCF: {final_merged_vcf}")
                            # Save checkpoint before exiting
                            state = {
                                'combination_set': combination_set,
                                'tools_pool': tools_pool,
                                'iteration': iteration - 1,  # We're retrying the current iteration
                                'iterations': iterations,
                                'current_combined_vcf': current_combined_vcf,
                                'current_combined_name': current_combined_name,
                                'vcf_files': vcf_files,
                                'metric': metric
                            }
                            save_state(state, args.checkpoint)
                            sys.exit(1)
                
                # Update the current combined VCF for next iteration
                current_combined_vcf = final_merged_vcf
                current_combined_name = combined_name
                vcf_files[current_combined_name] = current_combined_vcf
                
                # Evaluate the combined set - we already have the results from the evaluation step
                combined_result = tool_results[best_tool_name]
                
                # Store iteration results
                iteration_result = {
                    "iteration": iteration,
                    "tools": combined_name,
                    "tools_list": combination_set.copy(),
                    "TP": combined_result["TP"],
                    "FP": combined_result["FP"],
                    "recall": combined_result["recall"],
                    "precision": combined_result["precision"],
                    "f1": combined_result["f1"]
                }
            else:
                # Just use the best tool's results for the first iteration
                iteration_result = {
                    "iteration": iteration,
                    "tools": best_tool_name,
                    "tools_list": [best_tool_name],
                    "TP": best_tool_result["TP"],
                    "FP": best_tool_result["FP"],
                    "recall": best_tool_result["recall"],
                    "precision": best_tool_result["precision"],
                    "f1": best_tool_result["f1"]
                }
                
                # Set the current combined VCF to the best tool's VCF
                current_combined_vcf = vcf_files[best_tool_name]
                current_combined_name = best_tool_name
            
            iterations.append(iteration_result)
            
            # Save checkpoint after each successful iteration
            state = {
                'combination_set': combination_set,
                'tools_pool': tools_pool,
                'iteration': iteration,
                'iterations': iterations,
                'current_combined_vcf': current_combined_vcf,
                'current_combined_name': current_combined_name,
                'vcf_files': vcf_files,
                'metric': metric
            }
            save_state(state, args.checkpoint)
            
            # Clean up temporary files but keep the one we're using
            clean_temp_files(args.combined_dir)
            
            # Step 5: Remove the best tool from the pool
            tools_pool.remove(best_tool_name)
        
        # Generate final report
        logger.info("--- Final Results ---")
        df = pd.DataFrame(iterations)
        # Remove the tools_list column for display purposes
        display_df = df.drop('tools_list', axis=1) if 'tools_list' in df.columns else df
        logger.info("\n" + display_df.to_string())
        
        # Save results to CSV
        csv_filename = f"sv_combination_results_{metric}.csv"
        display_df.to_csv(csv_filename, index=False)
        logger.info(f"Results saved to {csv_filename}")
        
        # Log best combinations for all metrics regardless of which one was used for selection
        best_f1 = max(iterations, key=lambda x: x["f1"])
        logger.info(f"Best F1 score ({best_f1['f1']:.4f}) achieved with tools: {best_f1['tools']}")
        best_recall = max(iterations, key=lambda x: x["recall"])
        logger.info(f"Best recall ({best_recall['recall']:.4f}) achieved with tools: {best_recall['tools']}")
        best_precision = max(iterations, key=lambda x: x["precision"])
        logger.info(f"Best precision ({best_precision['precision']:.4f}) achieved with tools: {best_precision['tools']}")
        
        logger.info("SV combination evaluation complete")
        
        # Remove checkpoint file if everything completed successfully
        if os.path.exists(args.checkpoint):
            os.remove(args.checkpoint)
            logger.info(f"Removed checkpoint file: {args.checkpoint}")
        
    except KeyboardInterrupt:
        logger.warning("Process interrupted. Saving state before exiting...")
        state = {
            'combination_set': combination_set,
            'tools_pool': tools_pool,
            'iteration': iteration,
            'iterations': iterations,
            'current_combined_vcf': current_combined_vcf,
            'current_combined_name': current_combined_name,
            'vcf_files': vcf_files,
            'metric': metric
        }
        save_state(state, args.checkpoint)
        logger.info(f"State saved to {args.checkpoint}. Run with --resume to continue.")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error during execution: {e}")
        import traceback
        logger.error(traceback.format_exc())
        state = {
            'combination_set': combination_set,
            'tools_pool': tools_pool,
            'iteration': iteration,
            'iterations': iterations,
            'current_combined_vcf': current_combined_vcf,
            'current_combined_name': current_combined_name,
            'vcf_files': vcf_files,
            'metric': metric
        }
        save_state(state, args.checkpoint)
        logger.info(f"State saved to {args.checkpoint}. Run with --resume to continue after fixing the issue.")
        sys.exit(1)

if __name__ == "__main__":
    # Setup logging
    logger = setup_logging()
    main()