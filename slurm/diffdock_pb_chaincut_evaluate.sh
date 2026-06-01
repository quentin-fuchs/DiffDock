#!/bin/bash
#SBATCH --job-name=pb_chaincut_eval
#SBATCH --account=MPHIL-DIS-SL2-GPU
#SBATCH --partition=ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=/home/qf226/MProject/DiffDock/logs/pb_chaincut_eval_%j.out
#SBATCH --error=/home/qf226/MProject/DiffDock/logs/pb_chaincut_eval_%j.err

# Evaluate on 100 randomly sampled PoseBusters complexes using proteins
# pre-filtered to chains within 10Å of the ligand (protein_chaincut.pdb).
# Uses --dataset posebusters for proper multi-copy RMSD.
#
# Step 1: reindex ESM embeddings so chain indices match the trimmed proteins
#   (removing chains invalidates the original 0-based index ordering).
# Step 2: run evaluate.py with the corrected embeddings and a fresh cache.

DIFFDOCK_DIR=/home/qf226/MProject/DiffDock
OUT_DIR=/home/qf226/rds/hpc-work/results/pb_chaincut_eval

source ~/.bashrc
conda activate diffdock
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
export PYTHONPATH="$DIFFDOCK_DIR:$PYTHONPATH"

cd "$DIFFDOCK_DIR"
mkdir -p "$OUT_DIR"

# Generate split file from CSV (skip header, extract complex_name column)
awk -F',' 'NR>1 {print $3}' data/pb_chaincut_100.csv > data/pb_chaincut_100_split.txt
echo "Split file: $(wc -l < data/pb_chaincut_100_split.txt) complexes"

echo ""
echo "=== Step 1: Reindex ESM embeddings for chain-cut proteins ==="
python analysis/pb_reindex_chaincut_esm.py \
    --split_path data/pb_chaincut_100_split.txt \
    --data_dir   data/posebusters_benchmark_set \
    --esm_in     data/posebusters_esm2_embeddings.pt \
    --esm_out    data/posebusters_esm2_embeddings_chaincut.pt \
    --protein_orig protein \
    --protein_cut  protein_chaincut

echo ""
echo "=== Step 2: evaluate.py on 100 chain-cut PoseBusters complexes ==="
python evaluate.py \
    --config default_inference_args.yaml \
    --dataset posebusters \
    --data_dir data/posebusters_benchmark_set \
    --split_path data/pb_chaincut_100_split.txt \
    --cache_path data/cache_pb_chaincut \
    --out_dir "$OUT_DIR" \
    --esm_embeddings_path data/posebusters_esm2_embeddings_chaincut.pt \
    --protein_file protein_chaincut \
    --ligand_file ligands \
    --samples_per_complex 40 \
    --batch_size 40 \
    --save_predictions \
    --num_workers 4
