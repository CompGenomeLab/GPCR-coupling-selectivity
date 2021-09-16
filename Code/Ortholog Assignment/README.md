# Ortholog Assignment  
Ortholog assignment pipline contain two main sh files: ConservationPipeline_Part1.sh and ConservationPipeline_Part2.sh. In part 1, first maximum-likelihood tree is constructed. 

**ConservationPipeline_Part1.sh:**
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
    
    
    