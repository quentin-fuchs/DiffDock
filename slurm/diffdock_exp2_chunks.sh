#!/bin/bash
#SBATCH --job-name=dd_exp2_chunk
#SBATCH --account=MPHIL-DIS-SL2-GPU
#SBATCH --partition=ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --array=0-5
#SBATCH --output=/home/qf226/MProject/DiffDock/logs/dd_exp2_chunk_%A_%a.out
#SBATCH --error=/home/qf226/MProject/DiffDock/logs/dd_exp2_chunk_%A_%a.err

# Exp 2: run inference on all 6 timesplit_test chunks in parallel, each using
# the May-14 chunk cache (known-good conformers) at data/cache.
# Result should reproduce ~40% top-1 and show zero NaN failures.
# Aggregate per-complex results with analysis/aggregate_chunk_results.py afterwards.

CHUNK=${SLURM_ARRAY_TASK_ID}
DIFFDOCK_DIR=/home/qf226/MProject/DiffDock
SPLIT_PATH=$DIFFDOCK_DIR/data/splits/chunks/chunk_${CHUNK}.txt
OUT_DIR=/home/qf226/rds/hpc-work/results/pdbbind_exp2_chunks/chunk_${CHUNK}

source ~/.bashrc
conda activate diffdock
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
export PYTHONPATH="$DIFFDOCK_DIR:$PYTHONPATH"

cd "$DIFFDOCK_DIR"
mkdir -p "$OUT_DIR"

echo "[exp2] Chunk ${CHUNK}: $(wc -l < "$SPLIT_PATH") complexes"

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
