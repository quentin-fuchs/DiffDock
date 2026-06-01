#!/bin/bash
#SBATCH --job-name=check_esm
#SBATCH --account=MPHIL-DIS-SL2-CPU
#SBATCH --partition=icelake
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:10:00
#SBATCH --output=/home/qf226/MProject/DiffDock/logs/check_esm_%j.out
#SBATCH --error=/home/qf226/MProject/DiffDock/logs/check_esm_%j.err

source ~/.bashrc
conda activate diffdock
cd /home/qf226/MProject/DiffDock

python3 -c "
import torch, numpy as np

pb  = torch.load('data/posebusters_esm2_embeddings.pt', map_location='cpu')
pdb = torch.load('data/esm2_embeddings.pt', map_location='cpu')

pb_vals  = [v for v in pb.values()  if isinstance(v, torch.Tensor)]
pdb_vals = [v for v in pdb.values() if isinstance(v, torch.Tensor)]

print(f'PDBBind:     {len(pdb)} entries')
print(f'PoseBusters: {len(pb)} entries')
print()

# Shape and type
v0 = pdb_vals[0]; v1 = pb_vals[0]
print(f'PDBBind  shape example: {tuple(v0.shape)}, dtype={v0.dtype}')
print(f'PB       shape example: {tuple(v1.shape)}, dtype={v1.dtype}')
print()

# Check for empty tensors first
pb_empty  = [k for k, v in pb.items()  if isinstance(v, torch.Tensor) and v.numel() == 0]
pdb_empty = [k for k, v in pdb.items() if isinstance(v, torch.Tensor) and v.numel() == 0]
print(f'PDBBind empty tensors:  {len(pdb_empty)}')
print(f'PB      empty tensors:  {len(pb_empty)}')
if pb_empty:
    print(f'  PB empty keys (up to 10): {pb_empty[:10]}')
print()

# NaN / zero check for PB (skip empty tensors)
pb_nonempty  = [v for v in pb_vals  if v.numel() > 0]
pdb_nonempty = [v for v in pdb_vals if v.numel() > 0]
nan_count  = sum(1 for v in pb_nonempty if v.isnan().any())
zero_count = sum(1 for v in pb_nonempty if v.float().abs().max() < 1e-6)
print(f'PB non-empty:           {len(pb_nonempty)}/{len(pb_vals)}')
print(f'PB NaN embeddings:      {nan_count}/{len(pb_nonempty)}')
print(f'PB all-zero embeddings: {zero_count}/{len(pb_nonempty)}')
print()

# Norm distributions
def norm_stats(vals):
    norms = np.array([v.float().norm().item() for v in vals if v.numel() > 0])
    return f'mean={norms.mean():.1f}  std={norms.std():.1f}  min={norms.min():.1f}  max={norms.max():.1f}'

print(f'PDBBind norms: {norm_stats(pdb_vals)}')
print(f'PB      norms: {norm_stats(pb_nonempty)}')
print()

# Per-residue embedding dimension check (are they the same ESM layer?)
print(f'PDBBind embedding dim (per residue): {v0.shape[-1]}')
print(f'PB      embedding dim (per residue): {v1.shape[-1]}')
print()

# Check a few PB complexes that we KNOW are in pb_evaluate_merged
# (complexes where DiffDock ran successfully but scored poorly)
check_keys = ['5SAK_ZRY', '7OP9_06K', '7WL4_JFU', '7BJJ_TVW']
for k in check_keys:
    if k in pb:
        v = pb[k]
        print(f'{k}: shape={tuple(v.shape)}, norm={v.float().norm():.1f}, NaN={v.isnan().any().item()}')
    else:
        print(f'{k}: MISSING from PB embeddings!')
"
