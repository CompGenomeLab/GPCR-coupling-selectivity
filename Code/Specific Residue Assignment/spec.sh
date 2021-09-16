#!/bin/bash
#
# CompecTA (c) 2018
#
# makeblastdb_trial01 job submission script
#
# TODO:
#   - Set name of the job below changing "NAMD" value.
#   - Set the requested number of nodes (servers) with --nodes parameter.
#   - Set the requested number of tasks (cpu cores) with --ntasks parameter. (Total accross all nodes)
#   - Select the partition (queue) you want to run the job in:
#     - short : For jobs that have maximum run time of 120 mins. Has higher priority.
#     - mid   : For jobs that have maximum run time of 1 days. Lower priority than short.
#     - long  : For jobs that have maximum run time of 7 days. Lower priority than long.
#     - longer: For testing purposes, queue has 15 days limit but only 2 nodes.
#     - cuda  : For CUDA jobs. Solver that can utilize CUDA acceleration can use this queue. 7 days limit.
#   - Set the required time limit for the job with --time parameter.
#     - Acceptable time formats include "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"
#   - Put this script and all the input file under the same directory.
#   - Set the required parameters, input/output file names below.
#   - If you do not want mail please remove the line that has --mail-type and --mail-user. If you do want to get notification emails, set your email address.
#   - Put this script and all the input file under the same directory.
#   - Submit this file using:
#      sbatch slurm_submit.sh
#
# -= Resources =-
#
#SBATCH --job-name=spec
#SBATCH --account=investor
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --qos=mid_investor
#SBATCH --partition=mid_investor
#SBATCH --time=1-0:0
#SBATCH --output=%j-GNA15_Inoue-and-Avet.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bselcuk@sabanciuniv.edu
#SBATCH --open-mode=append


module load py-pip-19.3-gcc-9.2.0-tfsp3uc
python3.7 spec.py > GNA13comparison_13-12-2020.txt
