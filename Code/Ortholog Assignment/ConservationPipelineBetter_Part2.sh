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
#SBATCH --job-name=GPCRpip
#SBATCH --account=investor
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --qos=mid_investor
#SBATCH --partition=mid_investor
#SBATCH --time=1-0:0
#SBATCH --output=%j-GPCR-Pipeline.out
#SBATCH --mail-type=ALL
# #SBATCH --mail-user=bselcuk@sabanciuniv.edu
#SBATCH --array=0
#SBATCH --open-mode=append

################################################################################
source /etc/profile.d/modules.sh
echo "source /etc/profile.d/modules.sh"
################################################################################

echo ===== Obtaining gene clade =====
#protein_array=(DRD2 DRD3 DRD4 ADA2A ADA2B ADA2C HRH1 HRH3 HRH4 ACM1 ACM2 ACM3 ACM4 ACM5 5HT1A 5HT1B 5HT1D 5HT1E 5HT1F 5HT7R 5HT2B 5HT2C 5HT2A)
#date="19-06-2020"
#protein_array=(ADRB1 ADRB2 ADRB3)
#date=07-06-2020
#protein_array=(DRD1 DRD5 5HT6R)
#date=08-06-2020
#protein_array=(ADA1A ADA1B ADA1D)
#date=14-06-2020
#protein_array=(AA1R AA2AR AA2BR AA3R OPRM OPRX OPRK OPRD APJ AGTR1 GRPR BKRB1 BKRB2 NMBR CNR1 CNR2 CCKAR GASR GP119 GP132 GPR17 GP183 GPR34 GPR35 GPR55 GPR84 GP174 MRGX2 P2Y10 EDNRA EDNRB FPR1 FPR2 FFAR1 FFAR2 FFAR3 FFAR4 GHSR GNRHR KISSR CLTR1 CLTR2 FPR2 LT4R1 LT4R2 LPAR1 LPAR2 LPAR3 LPAR4 LPAR5 LPAR6 S1PR1 S1PR2 S1PR3 S1PR5 MCHR1 MCHR2 MSHR MC3R MC4R MC5R MTR1A MTR1B MTLR NMUR1 NMUR2 NPFF1 NPFF2 NPBW1 OX1R OX2R OXGR1 P2RY1 P2Y11 P2Y12 P2Y13 P2Y14 P2RY2 P2RY4 P2RY6 PTH1R PTH2R PTAFR PRLHR PD2R PE2R1 PE2R2 PE2R3 PE2R4 PF2R PI2R TA2R PAR1 PAR2 PAR3 PAR4 SSR1 SSR2 SSR3 SSR4 SSR5 NK3R NK2R NK1R UR2R PACR V1AR V1BR V2R OXYR)
#date=06-07-2020
#protein_array=(GALR1 GALR2 GALR3)
#date=13-07-2020
protein_array=(5HT5A)
date=06-01-2021
protein=${protein_array[$SLURM_ARRAY_TASK_ID]}
cd $protein"_"$date
module load py-pip-19.3-gcc-9.2.0-tfsp3uc
python ../gene_clade_find.py --tree $protein"_blasttree_"$date".treefile" --blastout $protein"_result_tabular.txt" --out $protein"_"$date"_subtree.nwk" 
python ../subtree_fastaV2.py --tree $protein"_"$date"_subtree.nwk" --out $protein"_"$date"_subtree.fasta"
echo ===== Second MSA =====
module load mafft-7.407-gcc-9.2.0-fww6zrx
mafft --ep 0 --thread -1 --genafpair --maxiterate 1000 $protein"_"$date"_subtree.fasta" > $protein"_subtreeMSA"$date".fasta"  
echo ===== Alignment Trimming =====
module load trimal-1.4.1-gcc-9.2.0-m27o55z
trimal -in $protein"_subtreeMSA"$date".fasta" -out $protein"_subtreeMSA"$date"_trimmed.fasta" -automated1
echo ===== Second Tree =====
module load raxml-ng-0.9.0
raxml-ng --threads 4 --search --data-type AA --model JTT+I+G4+F --force msa_dups --force msa_names --seed 2 -msa $protein"_subtreeMSA"$date"_trimmed.fasta"  --prefix $protein"_subtree_"$date # --redo #--bs-metric tbe,fbp
echo ===== Paralog Triming =====
module load py-pip-19.3-gcc-9.2.0-tfsp3uc
module load py-biopython-1.73-gcc-9.2.0-gzrw2dv
python ../ete3trimmerV2.py --tree $protein"_subtree_"$date".raxml.bestTree" --blastout $protein"_result_tabular.txt" --protein $protein --fasta $protein"_blast.fasta" --out $protein"_subtree_"$date"_orthologs.nwk"
python ../ete3order.py --tree $protein"_subtree_"$date"_orthologs.nwk" --out $protein"_subtree_"$date"_orthologs_ordered.nwk"
python ../subtree_fastaV2.py --tree $protein"_subtree_"$date"_orthologs_ordered.nwk"  --out $protein"_"$date"_orthologs_ordered.fasta"
