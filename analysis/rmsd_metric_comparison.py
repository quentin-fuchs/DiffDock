"""
Experiment A: compare spyrmsd (unaligned, symmetry-aware) vs RDKit GetBestRMS
(aligned, symmetry-aware) on the merged PB inference set.

Reads existing rmsds.npy (spyrmsd, computed by evaluate.py) and recomputes
RMSD for each complex's confidence-top-1 pose using RDKit alignment.
"""
import argparse
import numpy as np
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdMolAlign

MERGED_DIR   = Path("/home/qf226/rds/hpc-work/results/pb_evaluate_merged")
CRYSTAL_DIR  = Path("/home/qf226/MProject/DiffDock/data/posebusters_benchmark_set")
OUTPUT_NPZ   = MERGED_DIR / "rmsd_metric_comparison.npz"


def rdkit_rmsd(crystal_sdf: Path, pred_sdf: Path) -> float:
    mol_ref = Chem.MolFromMolFile(str(crystal_sdf), sanitize=True, removeHs=True)
    mol_prb = Chem.MolFromMolFile(str(pred_sdf),   sanitize=True, removeHs=True)
    if mol_ref is None or mol_prb is None:
        return float("nan")
    if mol_ref.GetNumAtoms() != mol_prb.GetNumAtoms():
        return float("nan")
    try:
        return rdMolAlign.GetBestRMS(mol_ref, mol_prb)
    except Exception:
        return float("nan")


def main():
    names      = np.load(MERGED_DIR / "complex_names.npy",  allow_pickle=True)
    rmsds_all  = np.load(MERGED_DIR / "rmsds.npy")           # [N, S]  spyrmsd, sample order
    confs_all  = np.load(MERGED_DIR / "confidences.npy")     # [N, S]

    N = len(names)
    print(f"Complexes: {N}, samples each: {rmsds_all.shape[1]}")

    spyrmsd_top1 = np.full(N, float("nan"))
    rdkit_top1   = np.full(N, float("nan"))
    skipped      = []

    for i, name in enumerate(names):
        # confidence-ranked top-1 index (same as rank1.sdf written by evaluate.py)
        top1_idx = int(np.argmax(confs_all[i]))
        spyrmsd_top1[i] = rmsds_all[i, top1_idx]

        crystal_sdf = CRYSTAL_DIR / name / f"{name}_ligand.sdf"
        pred_sdf    = MERGED_DIR  / name / "rank1.sdf"

        if not crystal_sdf.exists() or not pred_sdf.exists():
            skipped.append(name)
            continue

        rdkit_top1[i] = rdkit_rmsd(crystal_sdf, pred_sdf)

    valid = ~np.isnan(rdkit_top1)
    n_valid = valid.sum()
    n_skip  = np.isnan(rdkit_top1).sum()

    print(f"\nSuccessful: {n_valid}/{N}  |  Skipped/failed: {n_skip}")
    if skipped:
        print(f"  Missing files for: {skipped[:10]}{'...' if len(skipped)>10 else ''}")

    for thresh in [2.0, 5.0]:
        spy = (spyrmsd_top1[valid] < thresh).mean() * 100
        rdk = (rdkit_top1[valid]   < thresh).mean() * 100
        print(f"\nRMSD < {thresh} Å  (on {n_valid} complexes):")
        print(f"  spyrmsd (unaligned): {spy:.1f}%")
        print(f"  RDKit   (aligned):   {rdk:.1f}%")
        print(f"  delta:               {rdk - spy:+.1f} pp")

    # per-complex comparison
    flips_up   = ((spyrmsd_top1[valid] >= 2.0) & (rdkit_top1[valid] < 2.0)).sum()
    flips_down = ((spyrmsd_top1[valid] <  2.0) & (rdkit_top1[valid] >= 2.0)).sum()
    print(f"\nPer-complex flips at 2 Å threshold:")
    print(f"  spyrmsd>=2 → RDKit<2 (newly correct): {flips_up}")
    print(f"  spyrmsd<2  → RDKit>=2 (newly wrong):  {flips_down}")

    diffs = rdkit_top1[valid] - spyrmsd_top1[valid]
    print(f"\nRDKit - spyrmsd RMSD difference (on {n_valid} complexes):")
    print(f"  mean: {np.mean(diffs):+.3f} Å")
    print(f"  median: {np.median(diffs):+.3f} Å")
    print(f"  p25/p75: {np.percentile(diffs,25):+.3f} / {np.percentile(diffs,75):+.3f} Å")

    np.savez(OUTPUT_NPZ, names=names, spyrmsd_top1=spyrmsd_top1, rdkit_top1=rdkit_top1)
    print(f"\nSaved to {OUTPUT_NPZ}")


if __name__ == "__main__":
    main()
