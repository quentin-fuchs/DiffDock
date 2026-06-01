#!/bin/bash
#SBATCH --job-name=pb_full_eval
#SBATCH --account=MPHIL-DIS-SL2-GPU
#SBATCH --partition=ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=/home/qf226/MProject/DiffDock/logs/pb_full_eval_%j.out
#SBATCH --error=/home/qf226/MProject/DiffDock/logs/pb_full_eval_%j.err

# evaluate.py on the 308 curated PoseBusters complexes.
# Usage: sbatch diffdock_pb_full_evaluate.sh [split_path] [cache_path] [out_dir]
# Defaults to posebusters_pdb_set_correct.txt (308 complexes), cache_pb_eval, and
# /home/qf226/rds/hpc-work/results/pb_full_eval.

DIFFDOCK_DIR=/home/qf226/MProject/DiffDock
SPLIT_PATH=${1:-$DIFFDOCK_DIR/data/posebusters_pdb_set_correct.txt}
CACHE_PATH=${2:-data/cache_pb_eval}
OUT_DIR=${3:-/home/qf226/rds/hpc-work/results/pb_full_eval}

source ~/.bashrc
conda activate diffdock
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
export PYTHONPATH="$DIFFDOCK_DIR:$PYTHONPATH"

cd "$DIFFDOCK_DIR"
mkdir -p "$OUT_DIR"

echo "=== PoseBusters evaluate.py ==="
echo "Split:  $SPLIT_PATH ($(wc -l < "$SPLIT_PATH") complexes)"
echo "Cache:  $CACHE_PATH"
echo "Output: $OUT_DIR"
echo ""

python evaluate.py \
    --config default_inference_args.yaml \
    --dataset posebusters \
    --data_dir data/posebusters_benchmark_set \
    --split_path "$SPLIT_PATH" \
    --cache_path "$CACHE_PATH" \
    --out_dir "$OUT_DIR" \
    --esm_embeddings_path data/posebusters_esm2_embeddings.pt \
    --protein_file protein \
    --ligand_file ligands \
    --samples_per_complex 40 \
    --batch_size 40 \
    --save_predictions \
    --num_workers 4
