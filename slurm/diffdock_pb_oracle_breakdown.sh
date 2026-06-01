#!/bin/bash
#SBATCH --job-name=pb_oracle_bd
#SBATCH --account=MPHIL-DIS-SL2-CPU
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=/home/qf226/MProject/DiffDock/logs/pb_oracle_breakdown_%j.out
#SBATCH --error=/home/qf226/MProject/DiffDock/logs/pb_oracle_breakdown_%j.err

DIFFDOCK_DIR=/home/qf226/MProject/DiffDock

source ~/.bashrc
conda activate diffdock
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
export PYTHONPATH="$DIFFDOCK_DIR:$PYTHONPATH"

cd "$DIFFDOCK_DIR"

echo "=== PoseBusters oracle breakdown ==="
python analysis/pb_oracle_breakdown.py
