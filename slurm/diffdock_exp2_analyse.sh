#!/bin/bash
#SBATCH --job-name=dd_exp2_analyse
#SBATCH --account=MPHIL-DIS-SL2-CPU
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=/home/qf226/MProject/DiffDock/logs/dd_exp2_analyse_%j.out
#SBATCH --error=/home/qf226/MProject/DiffDock/logs/dd_exp2_analyse_%j.err

DIFFDOCK_DIR=/home/qf226/MProject/DiffDock
RESULTS_ROOT=/home/qf226/rds/hpc-work/results/pdbbind_exp2_chunks

source ~/.bashrc
conda activate diffdock
export PYTHONPATH="$DIFFDOCK_DIR:$PYTHONPATH"

cd "$DIFFDOCK_DIR"

echo "=== Exp 2: Full timesplit_test (May-14 caches) ==="
python analysis/aggregate_chunk_results.py "$RESULTS_ROOT" --n_chunks 6 --verbose
