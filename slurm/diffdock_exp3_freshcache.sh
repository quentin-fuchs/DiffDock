#!/bin/bash
#SBATCH --job-name=dd_exp3_fresh
#SBATCH --account=MPHIL-DIS-SL2-GPU
#SBATCH --partition=ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --array=0-5
#SBATCH --output=/home/qf226/MProject/DiffDock/logs/dd_exp3_fresh_%A_%a.out
#SBATCH --error=/home/qf226/MProject/DiffDock/logs/dd_exp3_fresh_%A_%a.err

# Exp 3: full timesplit_test (322 complexes, 6 chunks) with a brand-new cache
# directory (data/cache_freshtest), forcing a complete cache rebuild.
# Compare top-1 accuracy to Exp 2 (~38.8%) to test whether cache building
# is the source of any accuracy gap vs the paper's reported numbers.

CHUNK=${SLURM_ARRAY_TASK_ID}
DIFFDOCK_DIR=/home/qf226/MProject/DiffDock
SPLIT_PATH=$DIFFDOCK_DIR/data/splits/chunks/chunk_${CHUNK}.txt
OUT_DIR=/home/qf226/rds/hpc-work/results/pdbbind_exp3_freshcache/chunk_${CHUNK}

source ~/.bashrc
conda activate diffdock
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
export PYTHONPATH="$DIFFDOCK_DIR:$PYTHONPATH"

cd "$DIFFDOCK_DIR"
mkdir -p "$OUT_DIR"

echo "[exp3] Chunk ${CHUNK}: $(wc -l < "$SPLIT_PATH") complexes — fresh cache"

python evaluate.py \
    --config default_inference_args.yaml \
    --dataset pdbbind \
    --data_dir data/PDBBind_processed \
    --split_path "$SPLIT_PATH" \
    --cache_path data/cache_freshtest \
    --out_dir "$OUT_DIR" \
    --esm_embeddings_path data/esm2_embeddings.pt \
    --samples_per_complex 40 \
    --batch_size 40 \
    --save_predictions \
    --num_workers 4
