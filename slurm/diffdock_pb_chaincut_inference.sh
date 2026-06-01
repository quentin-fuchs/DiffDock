#!/bin/bash
#SBATCH --job-name=pb_chaincut_inf
#SBATCH --account=MPHIL-DIS-SL2-GPU
#SBATCH --partition=ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=/home/qf226/MProject/DiffDock/logs/pb_chaincut_inf_%j.out
#SBATCH --error=/home/qf226/MProject/DiffDock/logs/pb_chaincut_inf_%j.err

# Inference on 100 randomly sampled PoseBusters complexes using
# proteins pre-filtered to chains within 10Å of the ligand.
# Tests whether chain_cutoff (applied manually) improves accuracy
# vs our current 29% top-1 with full-protein loading.

DIFFDOCK_DIR=/home/qf226/MProject/DiffDock
OUT_DIR=/home/qf226/rds/hpc-work/results/pb_chaincut_100

source ~/.bashrc
conda activate diffdock
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
export PYTHONPATH="$DIFFDOCK_DIR:$PYTHONPATH"

cd "$DIFFDOCK_DIR"
mkdir -p "$OUT_DIR"

echo "=== Inference on 100 chain-cut PoseBusters complexes ==="
python inference.py \
    --config default_inference_args.yaml \
    --protein_ligand_csv data/pb_chaincut_100.csv \
    --out_dir "$OUT_DIR" \
    --esm_embeddings_path data/posebusters_esm2_embeddings.pt \
    --samples_per_complex 40 \
    --batch_size 40

echo ""
echo "=== RMSD evaluation ==="
python analysis/pb_rmsd_from_inference.py \
    --results_dir "$OUT_DIR" \
    --data_dir    data/posebusters_benchmark_set \
    --label       "chain_cutoff=10 (100 complexes)" \
    --n_samples   40
