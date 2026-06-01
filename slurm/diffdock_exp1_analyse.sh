#!/bin/bash
#SBATCH --job-name=dd_exp1_analyse
#SBATCH --account=MPHIL-DIS-SL2-CPU
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=/home/qf226/MProject/DiffDock/logs/dd_exp1_analyse_%j.out
#SBATCH --error=/home/qf226/MProject/DiffDock/logs/dd_exp1_analyse_%j.err

DIFFDOCK_DIR=/home/qf226/MProject/DiffDock
EXP1_DIR=/home/qf226/rds/hpc-work/results/pdbbind_eval_exp1_rebuild
BASELINE_DIR=/home/qf226/rds/hpc-work/results/pdbbind_eval_n100

source ~/.bashrc
conda activate diffdock
export PYTHONPATH="$DIFFDOCK_DIR:$PYTHONPATH"

cd "$DIFFDOCK_DIR"

echo "=== Exp 1: Fresh conformers vs May-30 baseline ==="
python analysis/compare_exp1_vs_baseline.py \
    --exp1_dir "$EXP1_DIR" \
    --baseline_dir "$BASELINE_DIR"
