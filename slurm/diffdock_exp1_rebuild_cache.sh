#!/bin/bash
#SBATCH --job-name=dd_exp1_cache
#SBATCH --account=MPHIL-DIS-SL2-GPU
#SBATCH --partition=ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/home/qf226/MProject/DiffDock/logs/dd_exp1_cache_%j.out
#SBATCH --error=/home/qf226/MProject/DiffDock/logs/dd_exp1_cache_%j.err

# Exp 1: rebuild the random-100 cache at a fresh path (data/cache_exp1) to force
# new RDKit conformer generation, then run inference with those conformers.
# Goal: check whether NaN failure count is reproducible (i.e. deterministic) or
# varies run-to-run (non-deterministic conformer generation is the suspected root cause).

DIFFDOCK_DIR=/home/qf226/MProject/DiffDock
SPLIT_PATH=$DIFFDOCK_DIR/data/splits/pdbbind_test_100_random.txt
OUT_DIR=/home/qf226/rds/hpc-work/results/pdbbind_eval_exp1_rebuild

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
    --cache_path data/cache_exp1 \
    --out_dir "$OUT_DIR" \
    --esm_embeddings_path data/esm2_embeddings.pt \
    --samples_per_complex 40 \
    --batch_size 40 \
    --save_predictions \
    --num_workers 4
