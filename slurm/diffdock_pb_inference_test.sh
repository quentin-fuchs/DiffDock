#!/bin/bash
#SBATCH --job-name=pb_infer_test
#SBATCH --account=MPHIL-DIS-SL2-GPU
#SBATCH --partition=ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/home/qf226/MProject/DiffDock/logs/pb_infer_test_%j.out
#SBATCH --error=/home/qf226/MProject/DiffDock/logs/pb_infer_test_%j.err

# Test: run inference.py directly (PoseBusters paper protocol) on 20 "always wrong"
# complexes (oracle > 10Å in evaluate.py run) using _ligand_start_conf.sdf.
# Compare RMSD to evaluate.py results to check whether the evaluate.py
# preprocessing (matching step, channel handling) is causing accuracy loss.

DIFFDOCK_DIR=/home/qf226/MProject/DiffDock
DATA_DIR=$DIFFDOCK_DIR/data/posebusters_benchmark_set
OUT_DIR=/home/qf226/rds/hpc-work/results/pb_infer_test
CSV_PATH=$DIFFDOCK_DIR/data/pb_infer_test_20.csv

source ~/.bashrc
conda activate diffdock
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
export PYTHONPATH="$DIFFDOCK_DIR:$PYTHONPATH"

cd "$DIFFDOCK_DIR"
mkdir -p "$OUT_DIR"

# Build CSV from the 20 always-wrong complexes (oracle > 10Å from evaluate.py run)
python3 - <<'PYEOF'
always_wrong_10 = [
    "7MSR_DCA", "7DQL_4CL", "8AY3_OE3", "7Q5I_I0F", "7MGY_ZD1",
    "7A1P_QW2", "7NSW_HC4", "7JY3_VUD", "7SUC_COM", "7MWU_ZPM",
    "7CL8_TES", "7EBG_J0L", "7XJN_NSD", "7NPL_UKZ", "7N6F_0I1",
    "7NP6_UK8", "7VBU_6I4", "7D6O_MTE", "7UJ5_DGL", "7BCP_GCO",
]
import os
DATA_DIR = "/home/qf226/MProject/DiffDock/data/posebusters_benchmark_set"
rows = ["protein_path,ligand_description,complex_name"]
for name in always_wrong_10:
    prot = os.path.join(DATA_DIR, name, f"{name}_protein.pdb")
    lig  = os.path.join(DATA_DIR, name, f"{name}_ligand_start_conf.sdf")
    if os.path.exists(prot) and os.path.exists(lig):
        rows.append(f"{prot},{lig},{name}")
    else:
        print(f"  WARNING: missing files for {name}")
with open("/home/qf226/MProject/DiffDock/data/pb_infer_test_20.csv", "w") as f:
    f.write("\n".join(rows) + "\n")
print(f"Wrote {len(rows)-1} entries to CSV")
PYEOF

echo "=== inference.py run on 20 always-wrong PoseBusters complexes ==="
echo "Using _ligand_start_conf.sdf (PoseBusters paper protocol)"
echo ""

python inference.py \
    --config default_inference_args.yaml \
    --protein_ligand_csv "$CSV_PATH" \
    --out_dir "$OUT_DIR" \
    --esm_embeddings_path data/posebusters_esm2_embeddings.pt \
    --samples_per_complex 40 \
    --batch_size 40

echo ""
echo "=== Computing oracle RMSD for inference.py results ==="
python3 - <<'PYEOF'
import os, sys, warnings, numpy as np
warnings.filterwarnings('ignore')
from rdkit import Chem
from rdkit.Chem import RemoveAllHs

sys.path.insert(0, '/home/qf226/MProject/DiffDock')

try:
    from utils.molecules_utils import get_symmetry_rmsd
except ImportError:
    def get_symmetry_rmsd(mol, ref_pos, pred_positions):
        # Fallback: simple RMSD
        return [np.sqrt(((p - ref_pos)**2).sum(axis=1).mean()) for p in pred_positions]

DATA_DIR = "/home/qf226/MProject/DiffDock/data/posebusters_benchmark_set"
OUT_DIR  = "/home/qf226/rds/hpc-work/results/pb_infer_test"

results = {}
for name in sorted(os.listdir(OUT_DIR)):
    complex_dir = os.path.join(OUT_DIR, name)
    if not os.path.isdir(complex_dir):
        continue

    # Load crystal reference
    crystal_path = os.path.join(DATA_DIR, name, f"{name}_ligand.sdf")
    if not os.path.exists(crystal_path):
        print(f"  {name}: no crystal ligand"); continue
    ref_mol = Chem.SDMolSupplier(crystal_path, sanitize=True, removeHs=True)[0]
    if ref_mol is None:
        print(f"  {name}: failed to parse crystal"); continue
    ref_pos = ref_mol.GetConformer().GetPositions()

    # Load all ranked poses
    pose_rmsds = []
    for rank in range(1, 41):
        sdf = os.path.join(complex_dir, f"rank{rank}.sdf")
        if not os.path.exists(sdf):
            break
        pred_mol = Chem.SDMolSupplier(sdf, sanitize=True, removeHs=True)[0]
        if pred_mol is None:
            pose_rmsds.append(np.inf); continue
        pred_pos = pred_mol.GetConformer().GetPositions()
        if pred_pos.shape != ref_pos.shape:
            pose_rmsds.append(np.inf); continue
        try:
            rmsd = get_symmetry_rmsd(ref_mol, ref_pos, [pred_pos])[0]
        except Exception:
            rmsd = np.sqrt(((pred_pos - ref_pos)**2).sum(axis=1).mean())
        pose_rmsds.append(rmsd)

    if not pose_rmsds:
        print(f"  {name}: no poses found"); continue

    rmsds = np.array(pose_rmsds)
    oracle = np.min(rmsds[np.isfinite(rmsds)]) if np.any(np.isfinite(rmsds)) else np.inf
    top1   = rmsds[0]
    results[name] = (top1, oracle, len(rmsds))
    print(f"  {name}: top1={top1:.2f}  oracle={oracle:.2f}  n_poses={len(rmsds)}")

if results:
    top1s   = [v[0] for v in results.values()]
    oracles = [v[1] for v in results.values()]
    print(f"\nSummary ({len(results)} complexes):")
    print(f"  Top-1 < 2Å : {100*np.mean(np.array(top1s)<2):.1f}%")
    print(f"  Oracle < 2Å: {100*np.mean(np.array(oracles)<2):.1f}%")
    print(f"  Mean oracle : {np.nanmean(oracles):.2f}Å")
PYEOF
