#!/bin/bash
#
# CompecTA (c) 2018

# -= Resources =-
#
#SBATCH --job-name=GPCRpip
#SBATCH --account=investor
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
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

# Module File
echo "Performing BLAST!"
protein_array=(5HT5A)
date=06-01-2021
protein=${protein_array[$SLURM_ARRAY_TASK_ID]}
mkdir $protein"_"$date 
cd $protein"_"$date
blastp="/cta/users/bselcuk/GPCRA/conservation_pipelineV2/ncbi-blast-2.6.0+/bin/blastp"
$blastp -query "../"$protein"_HUMAN.fasta" -max_target_seqs 1000 -num_threads 10 -db "/cta/users/bselcuk/GPCRA/blast_documents/databases/eukaryote_all" -outfmt "6 qseqid sseqid pident length mismatch gapopen send evalue bitscore stitle" -out $protein"_result_tabular.txt" 
echo "BLAST is done!"
echo ===========================
echo "Obtaining FASTA files to form tree..."
module load py-pip-19.3-gcc-9.2.0-tfsp3uc
python ../FastaObtainer.py --blastout $protein"_result_tabular.txt" --num 3 --taxid 9606 --out $protein"_blast.fasta" >>$protein"_"$date"_pipeline.log"
echo "FASTA file for first tree is produced."
echo "Making multiple sequence alignment with mafft."
module load mafft-7.407-gcc-9.2.0-fww6zrx
mafft --ep 0 --genafpair --thread -1 --maxiterate 1000  $protein"_blast.fasta" > $protein"_MSA"$date".fasta" 
echo "Alignment file is produced in FASTA format."
echo ===========================
echo "Starting to construct an evolutionary tree with IQtree version 1.6.12 with 1000 UFboot..."
iqtree2="/cta/users/bselcuk/GPCRA/conservation_pipelineV2/iqtree-2.0.6-Linux/bin/iqtree2"
$iqtree2 -s $protein"_MSA"$date".fasta" -m JTT+I+G4+F -B 1000 -bnni --prefix $protein"_blasttree_"$date -T 10 -wbt -seed 1 >>$protein"_"$date"_pipeline.log"


