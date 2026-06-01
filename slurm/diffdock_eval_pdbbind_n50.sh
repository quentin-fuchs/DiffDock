#!/bin/bash
#SBATCH --job-name=dd_pdb_n50
#SBATCH --account=MPHIL-DIS-SL2-GPU
#SBATCH --partition=ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/home/qf226/MProject/DiffDock/logs/dd_pdb_n50_%j.out
#SBATCH --error=/home/qf226/MProject/DiffDock/logs/dd_pdb_n50_%j.err

# 50 random PDBBind test complexes, 40 samples each.
# Split: pdbbind_test_50_random.txt (seed 42, sampled from timesplit_test)
# Cache rebuilt fresh after rdkit downgrade to 2022.03.3.

DIFFDOCK_DIR=/home/qf226/MProject/DiffDock
SPLIT_PATH=${1:-$DIFFDOCK_DIR/data/splits/pdbbind_test_50_random.txt}
OUT_DIR=${2:-/home/qf226/rds/hpc-work/results/pdbbind_eval_n50}

source ~/.bashrc
conda activate diffdock
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
export PYTHONPATH="$DIFFDOCK_DIR:$PYTHONPATH"

cd "$DIFFDOCK_DIR"
mkdir -p "$OUT_DIR"

python evaluate.py \
    --config default_inference_args.yaml \
    --dataset pdbbind \
    --data_dir data/PDBBind_processed \
    --split_path "$SPLIT_PATH" \
    --cache_path data/cache \
    --out_dir "$OUT_DIR" \
    --esm_embeddings_path data/esm2_embeddings.pt \
    --samples_per_complex 40 \
    --batch_size 40 \
    --save_predictions \
    --num_workers 4
