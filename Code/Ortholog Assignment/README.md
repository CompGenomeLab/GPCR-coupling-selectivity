# Ortholog Assignment  
Ortholog assignment pipline contain two main sh files: ConservationPipeline_Part1.sh and ConservationPipeline_Part2.sh. In part 1, first maximum-likelihood tree is constructed. 

**ConservationPipeline_Part1.sh with explanations of codes:**


```bash
#SBATCH --array=0     #This number should change based on number of proteins we analyze
protein_array=(5HT5A)   #Protein of interest. It can be a list also.
date=06-01-2021         #Date of the analysis
protein=${protein_array[$SLURM_ARRAY_TASK_ID]}  
mkdir $protein"_"$date  #Creating a folder for each protein
cd $protein"_"$date
blastp="/cta/users/bselcuk/GPCRA/conservation_pipelineV2/ncbi-blast-2.6.0+/bin/blastp" #blast version 2.6.0 used
$blastp -query "../"$protein"_HUMAN.fasta" -max_target_seqs 1000 -num_threads 10 -db "/cta/users/bselcuk/GPCRA/blast_documents/databases/eukaryote_all" -outfmt "6 qseqid sseqid pident length mismatch gapopen send evalue bitscore stitle" -out $protein"_result_tabular.txt"  
#blastp output in tabular format. Protein of interest was given in "$protein"_HUMAN.fasta" format
python ../FastaObtainer.py --blastout $protein"_result_tabular.txt" --num 3 --taxid 9606 --out $protein"_blast.fasta" >>$protein"_"$date"_pipeline.log" 
#FastaObtainer obtains fasta sequences of blast results. 
#--num indicates the number of human (--taxid 9606) sequences we obtain for tree construction
#--blastout is the input blastp result file produced in previous step
module load mafft-7.407-gcc-9.2.0-fww6zrx 
#MAFFT version 7.407 
mafft --ep 0 --genafpair --thread -1 --maxiterate 1000  $protein"_blast.fasta" > $protein"_MSA"$date".fasta" 
#IQtree2 version 2.0.6 was used to for tree construction with ultra-fast bootstrapping.
iqtree2="/cta/users/bselcuk/GPCRA/conservation_pipelineV2/iqtree-2.0.6-Linux/bin/iqtree2"
$iqtree2 -s $protein"_MSA"$date".fasta" -m JTT+I+G4+F -B 1000 -bnni --prefix $protein"_blasttree_"$date -T 10 -wbt -seed 1 >>$protein"_"$date"_pipeline.log"

```
**ConservationPipeline_Part2.sh with explanations of codes:**

```bash
#!/bin/bash
#SBATCH --array=0
protein_array=(5HT5A)
date=06-01-2021
protein=${protein_array[$SLURM_ARRAY_TASK_ID]}
cd $protein"_"$date
module load py-pip-19.3-gcc-9.2.0-tfsp3uc
python ../gene_clade_find.py --tree $protein"_blasttree_"$date".treefile" --blastout $protein"_result_tabular.txt" --out $protein"_"$date"_subtree.nwk"
#gene_clade_find uses ete3 to obtain the subtree for our protein of interest 
python ../subtree_fasta.py --tree $protein"_"$date"_subtree.nwk" --out $protein"_"$date"_subtree.fasta"
#subtree_fasta obtains fasta sequences for the subtree we obtained
module load mafft-7.407-gcc-9.2.0-fww6zrx
mafft --ep 0 --thread -1 --genafpair --maxiterate 1000 $protein"_"$date"_subtree.fasta" > $protein"_subtreeMSA"$date".fasta"  
echo ===== Alignment Trimming =====
module load trimal-1.4.1-gcc-9.2.0-m27o55z
#trimal version 1.4.1
trimal -in $protein"_subtreeMSA"$date".fasta" -out $protein"_subtreeMSA"$date"_trimmed.fasta" -automated1
echo ===== Second Tree =====
module load raxml-ng-0.9.0
#raxml-ng version 0.9.0
raxml-ng --threads 4 --search --data-type AA --model JTT+I+G4+F --force msa_dups --force msa_names --seed 2 -msa $protein"_subtreeMSA"$date"_trimmed.fasta"  --prefix $protein"_subtree_"$date # --redo #--bs-metric tbe,fbp
echo ===== Paralog Triming =====
module load py-pip-19.3-gcc-9.2.0-tfsp3uc
module load py-biopython-1.73-gcc-9.2.0-gzrw2dv
python ../ete3trimmer.py --tree $protein"_subtree_"$date".raxml.bestTree" --blastout $protein"_result_tabular.txt" --protein $protein --fasta $protein"_blast.fasta" --out $protein"_subtree_"$date"_orthologs.nwk"
#--tree input tree in newick format
#--blastout blast output in tabular format
#--protein our protein of interest
#--fasta fasta sequence source produced from blast results
#--out output trimmed tree in newick format

python ../ete3order.py --tree $protein"_subtree_"$date"_orthologs.nwk" --out $protein"_subtree_"$date"_orthologs_ordered.nwk"
#--tree input tree directory
#--out output tree directory
#Human leaf node is at the top of the tree.
python ../subtree_fasta.py --tree $protein"_subtree_"$date"_orthologs_ordered.nwk"  --out $protein"_"$date"_orthologs_ordered.fasta"
#Obtaining orthologous sequences in the same order with the ordered tree.
```