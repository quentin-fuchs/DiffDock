#!/bin/bash
#SBATCH --job-name=rmsd_metric_cmp
#SBATCH --account=FERGUSSON-SL3-CPU
#SBATCH --partition=sapphire
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/home/qf226/MProject/DiffDock/logs/rmsd_metric_cmp_%j.out
#SBATCH --error=/home/qf226/MProject/DiffDock/logs/rmsd_metric_cmp_%j.err

source ~/.bashrc
conda activate diffdock

cd /home/qf226/MProject/DiffDock
python3 analysis/rmsd_metric_comparison.py
