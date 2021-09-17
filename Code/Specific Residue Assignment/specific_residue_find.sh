#!/bin/bash
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
#SBATCH --output=%j-specific_residue_find.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bselcuk@sabanciuniv.edu
#SBATCH --open-mode=append


module load py-pip-19.3-gcc-9.2.0-tfsp3uc
python3.7 specific_residue_find.py > output_file.txt
