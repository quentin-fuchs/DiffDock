#!/bin/bash
#SBATCH --job-name=pb_file_audit
#SBATCH --account=MPHIL-DIS-SL2-CPU
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=/home/qf226/MProject/DiffDock/logs/pb_file_audit_%j.out
#SBATCH --error=/home/qf226/MProject/DiffDock/logs/pb_file_audit_%j.err

DIFFDOCK_DIR=/home/qf226/MProject/DiffDock

source ~/.bashrc
conda activate diffdock
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
export PYTHONPATH="$DIFFDOCK_DIR:$PYTHONPATH"

cd "$DIFFDOCK_DIR"

echo "=== PoseBusters file audit ==="
python analysis/pb_file_audit.py
