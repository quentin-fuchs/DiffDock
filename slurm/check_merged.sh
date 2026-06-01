#!/bin/bash
#SBATCH --job-name=check_merged
#SBATCH --account=MPHIL-DIS-SL2-CPU
#SBATCH --partition=icelake
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:10:00
#SBATCH --output=/home/qf226/MProject/DiffDock/logs/check_merged_%j.out
#SBATCH --error=/home/qf226/MProject/DiffDock/logs/check_merged_%j.err

source ~/.bashrc
conda activate diffdock

python3 - << 'PYEOF'
import numpy as np

merged_dir = "/home/qf226/rds/hpc-work/results/pb_evaluate_merged"
failed_dir = "/home/qf226/rds/hpc-work/results/pb_evaluate_out/chunk_failed"

# Load merged (chunks 0-4)
names_m = np.load(f"{merged_dir}/complex_names.npy", allow_pickle=True)
rmsds_m = np.load(f"{merged_dir}/rmsds.npy")
confs_m  = np.load(f"{merged_dir}/confidences.npy")
print(f"Merged (chunks 0-4): {len(names_m)} complexes, {rmsds_m.shape[1]} samples each")

# Load chunk_failed
names_f = np.load(f"{failed_dir}/complex_names.npy", allow_pickle=True)
rmsds_f = np.load(f"{failed_dir}/rmsds.npy")
confs_f  = np.load(f"{failed_dir}/confidences.npy")
print(f"Chunk_failed: {len(names_f)} complexes, {rmsds_f.shape[1]} samples each")

# Confidence top-1 for merged
ord_m = np.argsort(confs_m, axis=1)[:, ::-1]
top1_m = np.take_along_axis(rmsds_m, ord_m, axis=1)[:, 0]

# Confidence top-1 for chunk_failed
ord_f = np.argsort(confs_f, axis=1)[:, ::-1]
top1_f = np.take_along_axis(rmsds_f, ord_f, axis=1)[:, 0]

print(f"\nMerged top-1 (confidence): {(top1_m < 2).sum()}/{len(names_m)} = {100*(top1_m<2).mean():.1f}%")
print(f"Chunk_failed top-1 (confidence): {(top1_f < 2).sum()}/{len(names_f)} = {100*(top1_f<2).mean():.1f}%")

# Check overlap (are any chunk_failed names already in merged?)
overlap = set(names_m) & set(names_f)
print(f"\nOverlap between merged and chunk_failed: {len(overlap)} complexes")
if overlap:
    print(f"  Overlapping: {list(overlap)[:5]}")

# Combined (treating 3 unrecoverable as failures)
total_combined = len(names_m) + len(names_f)
successes_combined = (top1_m < 2).sum() + (top1_f < 2).sum()
print(f"\nCombined ({total_combined} complexes, 308-3={305} reachable):")
print(f"  Top-1 RMSD<2Å: {successes_combined}/{total_combined} = {100*successes_combined/total_combined:.1f}%")
print(f"  Top-1 RMSD<2Å (out of 308): {successes_combined}/308 = {100*successes_combined/308:.1f}%")

# Oracle for merged
oracle_m = (np.min(rmsds_m, axis=1) < 2).sum()
oracle_f = (np.min(rmsds_f, axis=1) < 2).sum()
print(f"\nOracle merged: {oracle_m}/{len(names_m)} = {100*oracle_m/len(names_m):.1f}%")
print(f"Oracle chunk_failed: {oracle_f}/{len(names_f)} = {100*oracle_f/len(names_f):.1f}%")
print(f"Oracle combined (out of 308): {oracle_m+oracle_f}/308 = {100*(oracle_m+oracle_f)/308:.1f}%")
PYEOF
