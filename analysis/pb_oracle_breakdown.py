"""
Break down PoseBusters oracle/top-1 by:
  1. Displaced vs non-displaced complexes (to see if centroid-offset hurts)
  2. n_copies in _ligands.sdf
  3. Complexes that are NEVER predicted correctly (oracle > 5Å) — who are they?
  4. Effect of excluding displaced on overall accuracy

Also: compute what oracle would need to be for 50% top-1 to be achievable,
assuming our confidence model ranking efficiency is ~74%.
"""

import os, numpy as np
from collections import defaultdict
from rdkit import Chem
from rdkit.Chem import RemoveAllHs

DATA_DIR   = "/home/qf226/MProject/DiffDock/data/posebusters_benchmark_set"
MERGED_DIR = "/home/qf226/rds/hpc-work/results/pb_evaluate_merged"

# --- load results ---
rmsds   = np.load(f"{MERGED_DIR}/rmsds.npy")
names   = np.load(f"{MERGED_DIR}/complex_names.npy", allow_pickle=True)
confs   = np.load(f"{MERGED_DIR}/confidences.npy")

conf_top1 = np.argmax(confs, axis=1)
top1_rmsd = rmsds[np.arange(len(rmsds)), conf_top1]
oracle    = np.min(rmsds, axis=1)

def ca_centroid(pdb_path):
    coords = []
    with open(pdb_path) as f:
        for line in f:
            if (line.startswith("ATOM") or line.startswith("HETATM")) and line[12:16].strip() == "CA":
                try:
                    coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                except ValueError:
                    pass
    return np.array(coords).mean(axis=0) if coords else None

def lig_centroid(sdf_path):
    supp = Chem.SDMolSupplier(sdf_path, sanitize=False, removeHs=False)
    mol = supp[0]
    if mol is None: return None
    try:
        Chem.SanitizeMol(mol)
        mol = RemoveAllHs(mol)
    except Exception: pass
    if mol.GetNumConformers() == 0: return None
    return mol.GetConformer().GetPositions().mean(axis=0)

def mol_count(path):
    try:
        return sum(1 for _ in Chem.SDMolSupplier(path, sanitize=False, removeHs=False))
    except Exception:
        return 0

print("Computing per-complex features...")
displacements = {}
n_copies_map  = {}
for name in names:
    d = os.path.join(DATA_DIR, name)
    prot = os.path.join(d, f"{name}_protein.pdb")
    ligs = os.path.join(d, f"{name}_ligands.sdf")
    lig  = os.path.join(d, f"{name}_ligand.sdf")
    ca = ca_centroid(prot)
    lc = lig_centroid(lig)
    if ca is not None and lc is not None:
        displacements[name] = float(np.linalg.norm(lc - ca))
    else:
        displacements[name] = np.nan
    n_copies_map[name] = mol_count(ligs)

dists  = np.array([displacements.get(n, np.nan) for n in names])
copies = np.array([n_copies_map.get(n, 0)       for n in names])

# --- split by displacement ---
CUTOFF = 30.0
displaced_mask = np.isfinite(dists) & (dists > CUTOFF)
normal_mask    = np.isfinite(dists) & (dists <= CUTOFF)

print(f"\n{'='*60}")
print(f"Overall: {len(names)} complexes")
print(f"  Top-1 < 2Å : {100*(top1_rmsd<2).mean():.2f}%  ({(top1_rmsd<2).sum()}/{len(top1_rmsd)})")
print(f"  Oracle < 2Å: {100*(oracle<2).mean():.2f}%  ({(oracle<2).sum()}/{len(oracle)})")

print(f"\nNormal (dist ≤ {CUTOFF}Å): n={normal_mask.sum()}")
print(f"  Top-1 < 2Å : {100*(top1_rmsd[normal_mask]<2).mean():.2f}%  ({(top1_rmsd[normal_mask]<2).sum()}/{normal_mask.sum()})")
print(f"  Oracle < 2Å: {100*(oracle[normal_mask]<2).mean():.2f}%  ({(oracle[normal_mask]<2).sum()}/{normal_mask.sum()})")

print(f"\nDisplaced (dist > {CUTOFF}Å): n={displaced_mask.sum()}")
print(f"  Top-1 < 2Å : {100*(top1_rmsd[displaced_mask]<2).mean():.2f}%  ({(top1_rmsd[displaced_mask]<2).sum()}/{displaced_mask.sum()})")
print(f"  Oracle < 2Å: {100*(oracle[displaced_mask]<2).mean():.2f}%  ({(oracle[displaced_mask]<2).sum()}/{displaced_mask.sum()})")

print(f"\nTop-1 efficiency (top1/oracle among oracle<2): ", end="")
or_pass = oracle < 2
if or_pass.sum() > 0:
    top1_when_oracle_passes = (top1_rmsd[or_pass] < 2).sum()
    print(f"{100*top1_when_oracle_passes/or_pass.sum():.1f}%  ({top1_when_oracle_passes}/{or_pass.sum()})")

# --- breakdown by n_copies ---
print(f"\n{'='*60}")
print(f"Breakdown by number of ligand copies in _ligands.sdf:")
for nc in sorted(set(copies)):
    mask = copies == nc
    n = mask.sum()
    if n == 0: continue
    t1 = (top1_rmsd[mask] < 2).mean()
    orc = (oracle[mask] < 2).mean()
    print(f"  n_copies={nc}: n={n}  top1={100*t1:.1f}%  oracle={100*orc:.1f}%")

# --- always wrong (oracle > 5Å) ---
print(f"\n{'='*60}")
always_bad = oracle > 5.0
print(f"Complexes always wrong (oracle > 5Å): {always_bad.sum()}/{len(oracle)}")
order = np.argsort(oracle[always_bad])[::-1]
bad_names = names[always_bad][order]
bad_oracle = oracle[always_bad][order]
bad_top1   = top1_rmsd[always_bad][order]
bad_dist   = dists[always_bad][order]
bad_copies = copies[always_bad][order]
print(f"{'Complex':15s} {'oracle':>8s} {'top1':>8s} {'dist':>7s} {'copies':>7s}")
for n, orc, t1, d, nc in zip(bad_names, bad_oracle, bad_top1, bad_dist, bad_copies):
    print(f"  {n:13s}  {orc:7.2f}  {t1:7.2f}  {d:7.1f}  {nc}")

# --- displacement distribution for oracle pass/fail ---
print(f"\n{'='*60}")
print(f"Displacement stats for oracle < 2Å vs oracle ≥ 2Å:")
pass_d = dists[oracle < 2]
fail_d = dists[oracle >= 2]
print(f"  Oracle pass: mean={np.nanmean(pass_d):.1f}Å  median={np.nanmedian(pass_d):.1f}Å  max={np.nanmax(pass_d):.1f}Å")
print(f"  Oracle fail: mean={np.nanmean(fail_d):.1f}Å  median={np.nanmedian(fail_d):.1f}Å  max={np.nanmax(fail_d):.1f}Å")

# --- what oracle would we need for 50% top-1 at 74% efficiency? ---
print(f"\n{'='*60}")
eff = (top1_rmsd[or_pass] < 2).mean()
print(f"Confidence model efficiency: {100*eff:.1f}%")
print(f"Oracle needed for 50% top-1 at {100*eff:.1f}% efficiency: {50/eff:.1f}%")
print(f"(Currently oracle = {100*(oracle<2).mean():.1f}%)")
print(f"Missing oracle complexes to reach 50% top-1: {round(0.50*len(oracle)/eff) - (oracle<2).sum()} more oracle-pass needed")
