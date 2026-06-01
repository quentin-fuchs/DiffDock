"""
Aggregate per-chunk evaluate.py outputs from a SLURM array job into a single
accuracy report.

Each chunk writes:
    {results_root}/chunk_{i}/rmsds.npy           shape (n_complexes, n_samples)
    {results_root}/chunk_{i}/centroid_distances.npy
    {results_root}/chunk_{i}/confidences.npy
    {results_root}/chunk_{i}/complex_names.npy   shape (n_complexes,)

Args:
    results_root: directory containing chunk_0/, chunk_1/, ... subdirectories
    n_chunks: number of chunks (default 6)

Returns (prints):
    - Total complexes evaluated
    - Top-1 filtered RMSD < 2 Å  (confidence-ranked)
    - Top-1 unfiltered RMSD < 2 Å (best-of-40 oracle)
    - Median and 75th percentile RMSD
    - Per-complex RMSD table for the worst 20
"""

import argparse
import os
import numpy as np


def load_chunk(chunk_dir):
    """Load numpy arrays for one chunk, return None if not yet written."""
    rmsds_path = os.path.join(chunk_dir, "rmsds.npy")
    names_path = os.path.join(chunk_dir, "complex_names.npy")
    conf_path = os.path.join(chunk_dir, "confidences.npy")
    centroid_path = os.path.join(chunk_dir, "centroid_distances.npy")
    if not os.path.exists(rmsds_path) or not os.path.exists(names_path):
        return None
    rmsds = np.load(rmsds_path)           # (n, samples)
    names = np.load(names_path)           # (n,)
    confidences = np.load(conf_path)      # (n, samples)
    centroids = np.load(centroid_path)    # (n, samples)
    return rmsds, names, confidences, centroids


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("results_root", help="Directory containing chunk_0/, chunk_1/, ...")
    parser.add_argument("--n_chunks", type=int, default=6)
    parser.add_argument("--rmsd_cutoff", type=float, default=2.0)
    parser.add_argument("--verbose", action="store_true", help="Print per-complex RMSDs")
    args = parser.parse_args()

    all_rmsds, all_names, all_confidences, all_centroids = [], [], [], []
    missing_chunks = []

    for i in range(args.n_chunks):
        chunk_dir = os.path.join(args.results_root, f"chunk_{i}")
        result = load_chunk(chunk_dir)
        if result is None:
            missing_chunks.append(i)
            print(f"[WARN] chunk {i} not yet complete — skipping")
            continue
        rmsds, names, confidences, centroids = result
        all_rmsds.append(rmsds)
        all_names.append(names)
        all_confidences.append(confidences)
        all_centroids.append(centroids)
        print(f"chunk {i}: {len(names)} complexes loaded")

    if not all_rmsds:
        print("No chunks complete yet.")
        return

    rmsds = np.concatenate(all_rmsds, axis=0)           # (N, samples)
    names = np.concatenate(all_names, axis=0)             # (N,)
    confidences = np.concatenate(all_confidences, axis=0) # (N, samples)
    centroids = np.concatenate(all_centroids, axis=0)     # (N, samples)

    N, S = rmsds.shape
    print(f"\nTotal complexes: {N} across {args.n_chunks - len(missing_chunks)} chunks")
    if missing_chunks:
        print(f"Missing chunks: {missing_chunks}")

    # Confidence-ranked top-1 (filtered): pick the sample with highest confidence
    conf_rank1_idx = np.argmax(confidences, axis=1)        # (N,)
    top1_rmsds = rmsds[np.arange(N), conf_rank1_idx]       # (N,)
    top1_centroids = centroids[np.arange(N), conf_rank1_idx]

    # Oracle top-1: best RMSD out of all samples
    oracle_top1 = np.min(rmsds, axis=1)                    # (N,)

    cutoff = args.rmsd_cutoff
    print(f"\n--- Top-1 (confidence-ranked) ---")
    print(f"  RMSD < {cutoff} Å : {100*(top1_rmsds < cutoff).mean():.2f}%  ({(top1_rmsds < cutoff).sum()}/{N})")
    print(f"  RMSD < 5 Å        : {100*(top1_rmsds < 5).mean():.2f}%")
    print(f"  Centroid < {cutoff} Å: {100*(top1_centroids < cutoff).mean():.2f}%")
    print(f"  Median RMSD       : {np.median(top1_rmsds):.2f} Å")
    print(f"  75th pct RMSD     : {np.percentile(top1_rmsds, 75):.2f} Å")

    print(f"\n--- Oracle top-1 (best of {S} samples) ---")
    print(f"  RMSD < {cutoff} Å : {100*(oracle_top1 < cutoff).mean():.2f}%  ({(oracle_top1 < cutoff).sum()}/{N})")
    print(f"  Median RMSD       : {np.median(oracle_top1):.2f} Å")

    # NaN detection: complexes where top-1 RMSD is very large (10000 sentinel from evaluate.py)
    nan_mask = top1_rmsds >= 9000
    if nan_mask.any():
        print(f"\n[WARN] {nan_mask.sum()} complexes had sentinel RMSD=10000 (all samples failed):")
        for name in names[nan_mask]:
            print(f"  {name}")

    if args.verbose:
        print(f"\n--- Per-complex top-1 RMSD (sorted) ---")
        order = np.argsort(top1_rmsds)[::-1]
        for idx in order:
            marker = "*" if top1_rmsds[idx] >= 9000 else (" " if top1_rmsds[idx] < cutoff else "X")
            print(f"  [{marker}] {names[idx]:12s}  rmsd={top1_rmsds[idx]:8.2f}  centroid={top1_centroids[idx]:8.2f}")


if __name__ == "__main__":
    main()
