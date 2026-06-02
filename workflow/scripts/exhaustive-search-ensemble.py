# module purge
# module load bcftools/1.22
# module load Python/3.12.3-GCCcore-13.3.0
# module load SciPy-bundle/2024.05-gfbf-2024a
# module load matplotlib/3.9.2-gfbf-2024a
# module load scikit-learn/1.5.2-gfbf-2024a


import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from itertools import combinations
from sklearn.tree import DecisionTreeClassifier, plot_tree, _tree
#from cyvcf2 import VCF
from sklearn.metrics import precision_score, recall_score, f1_score
import subprocess

tools = ['Manta','Dysgu', 'Octopus', 'DELLY','Smoove','GRIDSS', 'SvABA', 'Tiddit', 'LUMPY', 'PopDel', 'CNVpytor', 'Tardis', 'Wham']
sample = 'NA12878' 
truvari_res_dir = ''
results_dir = '' 
tp_file = 'tp-base.vcf.gz'
fp_file = 'fp.vcf.gz'
fn_file = 'fn.vcf.gz'

min_support = 2
max_combo_size = 10

os.makedirs(results_dir, exist_ok=True)

def parse_vcf(filepath):
    svs = set()
    try:
        cmd = ['bcftools', 'query', '-f', '%CHROM\t%POS\t%REF\t%ALT\n', filepath]
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        
        count = 0
        for line in process.stdout:
            count += 1
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[2] if fields[2] else "."
            alt = fields[3] if fields[3] else "."
            
            sv_id = (chrom, pos, ref, alt)
            svs.add(sv_id)
        
        process.wait()
        if process.returncode != 0:
            stderr_output = process.stderr.read()
            print(f"bcftools error: {stderr_output}")
        
        print(f"  Parsed {count} records from {os.path.basename(filepath)}")
        
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        import traceback
        traceback.print_exc()
    
    return svs


# ------------------------
# FP Clustering function 
# ------------------------
def cluster_fps(all_fp_svs, window=20):
    """
    Cluster FP SVs across all tools that share the same chromosome
    and are within `window` bp of each other.
    Maps each original (CHROM, POS, REF, ALT) → canonical representative.
    """
    if not all_fp_svs:
        return {}

    sorted_fps = sorted(all_fp_svs, key=lambda x: (x[0], x[1]))
    fp_to_canonical = {}
    current_cluster = [sorted_fps[0]]

    for sv in sorted_fps[1:]:
        chrom, pos = sv[0], sv[1]
        last_chrom, last_pos = current_cluster[-1][0], current_cluster[-1][1]

        if chrom == last_chrom and abs(pos - last_pos) <= window:
            current_cluster.append(sv)
        else:
            canonical = min(current_cluster, key=lambda x: x[1])
            for s in current_cluster:
                fp_to_canonical[s] = canonical
            current_cluster = [sv]

    # flush last cluster
    canonical = min(current_cluster, key=lambda x: x[1])
    for s in current_cluster:
        fp_to_canonical[s] = canonical

    return fp_to_canonical


# ------------------------
# Step 1: Load all tool results and construct truthset
# ------------------------

tool_results = {}
for tool in tools:
    tp_path = os.path.join(truvari_res_dir, tool, sample, tp_file)
    fp_path = os.path.join(truvari_res_dir, tool, sample, fp_file) 
    fn_path = os.path.join(truvari_res_dir, tool, sample, fn_file)
    
    tp_svs = parse_vcf(tp_path)
    fp_svs = parse_vcf(fp_path)
    fn_svs = parse_vcf(fn_path)
    
    tool_results[tool] = {
        'tp': tp_svs,
        'fp': fp_svs, 
        'fn': fn_svs,
        'calls': tp_svs.union(fp_svs)  # All calls made by this tool
    }
    
    print(f"{tool}: TP={len(tp_svs)}, FP={len(fp_svs)}, FN={len(fn_svs)}")

# Construct truthset from first tool (TP + FN)
first_tool = tools[0]
truthset = tool_results[first_tool]['tp'].union(tool_results[first_tool]['fn'])
print(f"\nTruthset size (from {first_tool} TP+FN): {len(truthset)}")


# ------------------------
# FP Clustering: merge nearby FP calls across tools
# ------------------------

all_fp_svs = set()
for tool in tools:
    all_fp_svs.update(tool_results[tool]['fp'])

print(f"  Unique FP loci pre-clustering:  {len(all_fp_svs)}")

fp_canonical_map = cluster_fps(all_fp_svs, window=50)

canonical_fps = set(fp_canonical_map.values())
print(f"  Unique FP loci post-clustering: {len(canonical_fps)}")
print(f"  FPs collapsed by clustering:    {len(all_fp_svs) - len(canonical_fps)}")

# Update each tool's calls to use clustered FP representatives
for tool in tools:
    fp_clustered = {fp_canonical_map[sv] for sv in tool_results[tool]['fp']}
    tool_results[tool]['calls'] = tool_results[tool]['tp'].union(fp_clustered)


# ------------------------
# Step 2: Evaluate individual tools
# ------------------------
print("\nEvaluating individual tools...")

combo_results = []

for tool in tools:
    calls = tool_results[tool]['calls']
    
    # Calculate metrics against truthset
    TP = len(calls.intersection(truthset))
    FP = len(calls - truthset)
    FN = len(truthset - calls)
    
    precision = TP / (TP + FP) if TP + FP > 0 else 0.0
    recall = TP / (TP + FN) if TP + FN > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall > 0 else 0.0
    
    combo_results.append({
        'combo': tool,
        'num_tools': 1,
        'TP': TP,
        'FP': FP,
        'FN': FN,
        'precision': precision,
        'recall': recall,
        'f1': f1,
        'total_calls': len(calls)
    })
    
    print(f"{tool}: P={precision:.3f}, R={recall:.3f}, F1={f1:.3f}")

# ------------------------
# Step 3: Evaluate tool combinations
# ------------------------
print(f"\nEvaluating tool combinations (min_support={min_support})...")

for r in range(2, max_combo_size + 1):
    print(f"Evaluating {r}-tool combinations...")
    
    for combo in combinations(tools, r):
        combo_name = '+'.join(combo)
        
        # Get union of calls from all tools in combination
        union_calls = set()
        for tool in combo:
            union_calls.update(tool_results[tool]['calls'])
        
        # Calculate metrics against truthset
        TP = len(union_calls.intersection(truthset))
        FP = len(union_calls - truthset)
        FN = len(truthset - union_calls)
        
        precision = TP / (TP + FP) if TP + FP > 0 else 0.0
        recall = TP / (TP + FN) if TP + FN > 0 else 0.0
        f1 = 2 * precision * recall / (precision + recall) if precision + recall > 0 else 0.0
        
        combo_results.append({
            'combo': combo_name,
            'num_tools': r,
            'TP': TP,
            'FP': FP,
            'FN': FN,
            'precision': precision,
            'recall': recall,
            'f1': f1,
            'total_calls': len(union_calls)
        })

# ------------------------ 
# Step 4: Alternative - Voting-based combinations
# ------------------------
print(f"\nEvaluating voting-based combinations (min_support={min_support})...")

# Create SV matrix for voting approach
all_svs = set()
for tool in tools:
    all_svs.update(tool_results[tool]['calls'])
all_svs.update(truthset)

sv_matrix = {}
for sv in all_svs:
    is_truth = 1 if sv in truthset else 0
    sv_matrix[sv] = {'truth': is_truth, **{t: 0 for t in tools}}

# Populate tool columns
for tool in tools:
    calls = tool_results[tool]['calls']
    for sv in calls:
        if sv in sv_matrix:
            sv_matrix[sv][tool] = 1

# Convert to DataFrame
df = pd.DataFrame([{'sv_key': sv, **vals} for sv, vals in sv_matrix.items()])


voting_results = []
for r in range(2, max_combo_size + 1):
    for combo in combinations(tools, r):
        combo_name = '+'.join(combo) + f'_vote{min_support}'
        vote = df[list(combo)].sum(axis=1)

        # TP: consensus required
        TP = ((vote >= min_support) & (df['truth'] == 1)).sum()
        # FN: in truth but didn't reach consensus
        FN = ((vote < min_support)  & (df['truth'] == 1)).sum()
        # FP: any tool in combo called it, regardless of support
        FP = ((vote >= 1)           & (df['truth'] == 0)).sum()

        precision = TP / (TP + FP) if TP + FP > 0 else 0.0
        recall    = TP / (TP + FN) if TP + FN > 0 else 0.0
        f1 = 2 * precision * recall / (precision + recall) if precision + recall > 0 else 0.0

        voting_results.append({
            'combo': combo_name,
            'num_tools': r,
            'TP': TP,
            'FP': FP,
            'FN': FN,
            'precision': precision,
            'recall': recall,
            'f1': f1,
            'total_calls': TP + FP
        })

# ------------------------
# Step 5: Save results separately
# ------------------------
# Save union results
union_df = pd.DataFrame(combo_results).sort_values(by='f1', ascending=False)
union_df.to_csv(os.path.join(results_dir, 'union_performance.csv'), index=False)

# Save voting results  
voting_df = pd.DataFrame(voting_results).sort_values(by='f1', ascending=False)
voting_df.to_csv(os.path.join(results_dir, 'voting_performance.csv'), index=False)

# Combine results for plotting
all_results = combo_results + voting_results
combo_df = pd.DataFrame(all_results).sort_values(by='f1', ascending=False)
combo_df.to_csv(os.path.join(results_dir, 'combo_performance.csv'), index=False)

# ------------------------
# Step 6: Analyze caller patterns
# ------------------------
print(f"\n=== CALLER PATTERN ANALYSIS ===")

# Correlation analysis between tools
print("\n Tool Correlation Analysis:")
# Calculate pairwise correlation between tools
tool_correlations = df[tools].corr()

print("\n True Positive Correlation Analysis:")
tp_matrix = df[df['truth'] == 1][tools]  # Only true SVs
tp_correlations = tp_matrix.corr()

tp_corr_pairs = []
for i in range(len(tools)):
    for j in range(i+1, len(tools)):
        tp_corr_val = tp_correlations.iloc[i, j]
        tp_corr_pairs.append((tools[i], tools[j], tp_corr_val))


# Save correlation matrix
tool_correlations.to_csv(os.path.join(results_dir, 'tool_correlations.csv'))
tp_correlations.to_csv(os.path.join(results_dir, 'tp_correlations.csv'))

# Save SV matrix
df.to_csv(os.path.join(results_dir, 'sv_matrix.csv'), index=False)


# ------------------------
# Visualization: Tool Correlation Heatmaps (Side by Side)
# ------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))

# All SVs correlation heatmap
mask = np.triu(np.ones_like(tool_correlations, dtype=bool))
im1 = ax1.imshow(tool_correlations.values, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')

# Add colorbar for all SVs
cbar1 = plt.colorbar(im1, ax=ax1, shrink=0.8)
cbar1.set_label('Correlation Coefficient', rotation=270, labelpad=20)

# Set ticks and labels for all SVs
ax1.set_xticks(range(len(tools)))
ax1.set_xticklabels(tools, rotation=45, ha='right')
ax1.set_yticks(range(len(tools)))
ax1.set_yticklabels(tools)

# Add correlation values as text for all SVs
for i in range(len(tools)):
    for j in range(len(tools)):
        if not mask[i, j]:  # Only show lower triangle
            text = ax1.text(j, i, f'{tool_correlations.iloc[i, j]:.2f}',
                          ha="center", va="center", 
                          color="white" if abs(tool_correlations.iloc[i, j]) > 0.5 else "black",
                          fontsize=8)

ax1.set_title('All SVs Correlation Matrix\n(Including TPs and FPs)', fontsize=12, pad=15)

# TP-only correlation heatmap
im2 = ax2.imshow(tp_correlations.values, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')

# Add colorbar for TP only
cbar2 = plt.colorbar(im2, ax=ax2, shrink=0.8)
cbar2.set_label('Correlation Coefficient', rotation=270, labelpad=20)

# Set ticks and labels for TP only
ax2.set_xticks(range(len(tools)))
ax2.set_xticklabels(tools, rotation=45, ha='right')
ax2.set_yticks(range(len(tools)))
ax2.set_yticklabels(tools)

# Add correlation values as text for TP only
for i in range(len(tools)):
    for j in range(len(tools)):
        if not mask[i, j]:  # Only show lower triangle
            text = ax2.text(j, i, f'{tp_correlations.iloc[i, j]:.2f}',
                          ha="center", va="center", 
                          color="white" if abs(tp_correlations.iloc[i, j]) > 0.5 else "black",
                          fontsize=8)

ax2.set_title('True Positives Only Correlation Matrix\n(Excluding False Positives)', fontsize=12, pad=15)

plt.tight_layout()
plt.savefig(os.path.join(results_dir, 'tool_correlation_comparison.pdf'), dpi=300, bbox_inches='tight')
plt.close()

