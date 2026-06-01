"""
Compare Exp 1 (fresh cache, new conformers) against the May-30 baseline run
to determine whether NaN failures are caused by conformer non-determinism or
are an intrinsic property of certain complexes.

Usage:
    python analysis/compare_exp1_vs_baseline.py \
        --exp1_dir /home/qf226/rds/hpc-work/results/pdbbind_eval_exp1_rebuild \
        --baseline_dir /home/qf226/rds/hpc-work/results/pdbbind_eval_<jobid>

Reports:
    - Overall accuracy for both runs
    - Set of NaN-sentinel complexes in each run
    - Overlap: same complexes failing in both → intrinsic
    - Different complexes failing → conformer-dependent
"""

import argparse
import os
import numpy as np


NAN_SENTINEL = 9000  # evaluate.py uses 10000 for all-failed complexes


def load_run(run_dir):
    """Load rmsds and complex names from an evaluate.py output directory."""
    rmsds = np.load(os.path.join(run_dir, "rmsds.npy"))
    names = np.load(os.path.join(run_dir, "complex_names.npy"))
    confidences = np.load(os.path.join(run_dir, "confidences.npy"))
    return rmsds, names, confidences


def top1_stats(rmsds, confidences, label, cutoff=2.0):
    N, S = rmsds.shape
    conf_rank1_idx = np.argmax(confidences, axis=1)
    top1 = rmsds[np.arange(N), conf_rank1_idx]
    nan_mask = top1 >= NAN_SENTINEL
    valid = top1[~nan_mask]
    print(f"\n=== {label} (n={N}) ===")
    print(f"  Top-1 RMSD < {cutoff} Å : {100*(top1 < cutoff).mean():.2f}%  ({(top1 < cutoff).sum()}/{N})")
    print(f"  NaN/failed complexes    : {nan_mask.sum()}")
    print(f"  Median RMSD (valid)     : {np.median(valid):.2f} Å" if len(valid) else "  Median RMSD: N/A")
    return top1, nan_mask


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--exp1_dir", required=True, help="Fresh-cache run output dir")
    parser.add_argument("--baseline_dir", required=True, help="May-30 baseline run output dir")
    parser.add_argument("--cutoff", type=float, default=2.0)
    args = parser.parse_args()

    exp1_rmsds, exp1_names, exp1_confs = load_run(args.exp1_dir)
    base_rmsds, base_names, base_confs = load_run(args.baseline_dir)

    exp1_top1, exp1_nan = top1_stats(exp1_rmsds, exp1_confs, "Exp 1 (fresh cache)", args.cutoff)
    base_top1, base_nan = top1_stats(base_rmsds, base_confs, "Baseline (May-30)", args.cutoff)

    exp1_nan_set = set(exp1_names[exp1_nan])
    base_nan_set = set(base_names[base_nan])

    both_fail = exp1_nan_set & base_nan_set
    only_exp1 = exp1_nan_set - base_nan_set
    only_base = base_nan_set - exp1_nan_set

    print(f"\n--- NaN overlap analysis ---")
    print(f"  Both fail (intrinsic)       : {len(both_fail)}  {sorted(both_fail)}")
    print(f"  Only Exp1 fails (new cache)  : {len(only_exp1)}  {sorted(only_exp1)}")
    print(f"  Only baseline fails (old cache): {len(only_base)}  {sorted(only_base)}")

    if len(both_fail) == len(exp1_nan_set) == len(base_nan_set):
        print("\nConclusion: NaN failures are DETERMINISTIC — same complexes fail regardless of conformer.")
    elif len(both_fail) == 0:
        print("\nConclusion: NaN failures are FULLY NON-DETERMINISTIC — different conformers cause different failures.")
    else:
        print(f"\nConclusion: MIXED — {len(both_fail)} intrinsic failures, some conformer-dependent variation.")

    # Per-complex RMSD comparison for all complexes in both runs
    exp1_dict = dict(zip(exp1_names, exp1_top1))
    base_dict = dict(zip(base_names, base_top1))
    common = sorted(set(exp1_names) & set(base_names))
    if common:
        print(f"\n--- Per-complex RMSD comparison (n={len(common)} common complexes) ---")
        diffs = [abs(exp1_dict[n] - base_dict[n]) for n in common if exp1_dict[n] < NAN_SENTINEL and base_dict[n] < NAN_SENTINEL]
        if diffs:
            print(f"  Mean |RMSD diff|  : {np.mean(diffs):.3f} Å")
            print(f"  Max  |RMSD diff|  : {np.max(diffs):.3f} Å")
            print(f"  Complexes where outcome flipped (one <2Å, other >=2Å):")
            flipped = [n for n in common
                       if exp1_dict[n] < NAN_SENTINEL and base_dict[n] < NAN_SENTINEL
                       and (exp1_dict[n] < args.cutoff) != (base_dict[n] < args.cutoff)]
            for n in flipped:
                print(f"    {n}  exp1={exp1_dict[n]:.2f}  base={base_dict[n]:.2f}")
            if not flipped:
                print("    (none — conformer variation does not change pass/fail)")


if __name__ == "__main__":
    main()
