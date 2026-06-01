#!/bin/bash
#SBATCH --job-name=pb_chaincut_pre
#SBATCH --account=MPHIL-DIS-SL2-CPU
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=/home/qf226/MProject/DiffDock/logs/pb_chaincut_pre_%j.out
#SBATCH --error=/home/qf226/MProject/DiffDock/logs/pb_chaincut_pre_%j.err

DIFFDOCK_DIR=/home/qf226/MProject/DiffDock

source ~/.bashrc
conda activate diffdock
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
export PYTHONPATH="$DIFFDOCK_DIR:$PYTHONPATH"

cd "$DIFFDOCK_DIR"

echo "=== PoseBusters chain cutoff preprocessing (cutoff=10Å) ==="
python analysis/pb_preprocess_chaincut.py \
    --data_dir data/posebusters_benchmark_set \
    --split    data/posebusters_pdb_set_correct.txt \
    --cutoff   10.0 \
    --csv_out  data/pb_chaincut_100.csv \
    --n_sample 100

echo "Done. CSV written to data/pb_chaincut_100.csv"
