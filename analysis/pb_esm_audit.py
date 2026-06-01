"""
Check ESM embedding coverage for PoseBusters benchmark.

Verifies:
  1. How many PB complexes have ESM embeddings in posebusters_esm2_embeddings.pt
  2. How many PB complexes in our split file are missing embeddings
  3. Chain ID coverage (protein.pdb chains vs embedding chains)
  4. Cross-check: which complexes in the split file were skipped during evaluation
     due to ESM mismatch (tensor size error in logs)
"""

import os, sys, torch, numpy as np
from collections import defaultdict

DATA_DIR    = "/home/qf226/MProject/DiffDock/data/posebusters_benchmark_set"
EMB_PATH    = "/home/qf226/MProject/DiffDock/data/posebusters_esm2_embeddings.pt"
SPLIT_PATH  = "/home/qf226/MProject/DiffDock/data/posebusters_pdb_set_correct.txt"
MERGED_NAMES = "/home/qf226/rds/hpc-work/results/pb_evaluate_merged/complex_names.npy"

# Load split list
with open(SPLIT_PATH) as f:
    split_ids = [l.strip() for l in f if l.strip()]
print(f"Split IDs: {len(split_ids)}")

# Load evaluated IDs
eval_names = set(np.load(MERGED_NAMES, allow_pickle=True))
print(f"Evaluated:  {len(eval_names)}")

# Missing from evaluation
missing_from_eval = [n for n in split_ids if n not in eval_names]
print(f"In split but NOT evaluated: {len(missing_from_eval)}")
for n in missing_from_eval:
    print(f"  {n}")

# Load ESM embeddings
print(f"\nLoading ESM embeddings from {EMB_PATH}...")
emb = torch.load(EMB_PATH, map_location='cpu')
print(f"  Keys in embedding dict: {len(emb)}")
print(f"  Example keys: {list(emb.keys())[:5]}")

# Check coverage
emb_keys = set(emb.keys())
split_set = set(split_ids)
has_emb   = split_set & emb_keys
no_emb    = split_set - emb_keys
print(f"\nESM embedding coverage for split set:")
print(f"  Have embeddings: {len(has_emb)} / {len(split_ids)}")
print(f"  Missing embeddings: {len(no_emb)}")
for n in sorted(no_emb):
    print(f"    {n}")

# Chain mismatch check: compare chains in protein.pdb vs embedding chains
print("\nChecking chain coverage for first 20 complexes...")
chain_issues = []
for name in split_ids[:20]:
    pdb_path = os.path.join(DATA_DIR, name, f"{name}_protein.pdb")
    if not os.path.exists(pdb_path):
        continue
    pdb_chains = set()
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and len(line) > 21:
                pdb_chains.add(line[21])

    if name in emb:
        emb_chains = list(emb[name].keys()) if isinstance(emb[name], dict) else ["(not_dict)"]
    else:
        emb_chains = ["MISSING"]

    print(f"  {name}: PDB chains={sorted(pdb_chains)}  ESM chains={emb_chains}")

# Distribution of embedding shapes
print("\nEmbedding shape distribution (first 10):")
for i, (k, v) in enumerate(emb.items()):
    if i >= 10:
        break
    if isinstance(v, dict):
        shapes = {ch: arr.shape for ch, arr in v.items()}
        print(f"  {k}: {shapes}")
    elif isinstance(v, torch.Tensor):
        print(f"  {k}: tensor {v.shape}")
    else:
        print(f"  {k}: {type(v)}")
