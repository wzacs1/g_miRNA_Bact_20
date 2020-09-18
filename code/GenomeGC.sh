#!/bin/bash

#SBATCH --account=round
#SBATCH --partition=notchpeak-shared
#SBATCH -J par_GC
#SBATCH -n 24
#SBATCH -t 48:00:00
#SBATCH -D /uufs/chpc.utah.edu/common/home/u0210816/Projects/g_miRNA_Bact_20
#SBATCH -o /uufs/chpc.utah.edu/common/home/u0210816/Projects/g_miRNA_Bact_20/code/GenomeGC_parallel.outerror

# Purpose:
# This script is written to allow functions to go into the original "miRNA_Scan_Count" scripts, which it should have already. Thus I maintain the initial setup and genome copy commands until new commands noted.

######### User-defined variables ######################
# Working directory to store results.  (Created if not present)
WrkDir=/uufs/chpc.utah.edu/common/home/u0210816/Projects/g_miRNA_Bact_20
GenomeSetName=all_pangenomes_parallel
# Scratch directory. (Created if not present). Potentially large amount of space required for genomes copy.
SCRATCH=/scratch/general/nfs1/u0210816/miRNA
# Target seq to search for. In future, iterate through single full sequence instead.
TargetSeqA=AGTTCTCA
TargetSeqB=GTTCTCA
TargetSeqC=TTCTCA
TargetSeqD=AGTTCT
TargetSeqE=AGTTCTC
TargetSeqF=GTTCTC
# Pangenomes files to search. Text file with one genome (gzipped) per line.
Genomes=/uufs/chpc.utah.edu/common/home/u0210816/Projects/g_miRNA_Bact_20/data/genome_lists/full_chocophlan.txt
# Location of reference genomes to search
GenomesDir=/uufs/chpc.utah.edu/common/home/round-group1/reference_seq_dbs/humann2/chocophlan
# Number of parallel jobs to start (2X processes should be okay)
NumProc=48
#######################################################

########### Setup #####################################
mkdir -p $SCRATCH
mkdir -p $WrkDir
# Module loading
module load fastx_toolkit
#######################################################

mkdir -p ${SCRATCH}/${GenomeSetName}/ref_genomes

while read line
    do
    cp ${GenomesDir}/${line} ${SCRATCH}/${GenomeSetName}/ref_genomes/
    gunzip -c ${SCRATCH}/${GenomeSetName}/ref_genomes/${line} | fasta_formatter -t -o ${SCRATCH}/${GenomeSetName}/ref_genomes/${line%.gz}.tab
done < $Genomes

cd ${SCRATCH}/${GenomeSetName}/ref_genomes

### Begin new commands to ID GC content and size of search space (not exactly genome b/c only coding regions included in .fna files, but should be nearly identical to genome)

# First, collate the sequences together to form a temp genome file

GenomeStats () {
    GenomeName=${1%.centroids.v0.1.1.ffn.tab}
    echo "GenomeName,A_content,T_content,G_content,C_content,Size" > ${GenomeName}_GenomeStats.csv
    cut -f 2 $1 | tr -d '\n' > tmp_${GenomeName}.txt
    Size=`cat tmp_${GenomeName}.txt | wc -m`
    Acontent=`cat tmp_${GenomeName}.txt | tr -d "GCT" | wc -m`
    Gcontent=`cat tmp_${GenomeName}.txt | tr -d "TCA" | wc -m`
    Ccontent=`cat tmp_${GenomeName}.txt | tr -d "GAT" | wc -m`
    Tcontent=`cat tmp_${GenomeName}.txt | tr -d "GCA" | wc -m`
    echo "$GenomeName,$Acontent,$Gcontent,$Ccontent,$Tcontent,$Size" >> ${GenomeName}_GenomeStats.csv
    rm tmp_${GenomeName}.txt
}

export -f GenomeStats

# Run Command in parallel

ls -1 --color=never *.tab | parallel -j $NumProc "GenomeStats {}"

# Collate results to one file:
echo "GenomeName,A_content,T_content,G_content,C_content,Size" > all_GenomeStats.tmp
for f in *_GenomeStats.csv
    do
    tail -n +2 $f >> all_GenomeStats.tmp
done
mv all_GenomeStats.tmp all_GenomeStats.csv

# Cleanup, copy results to work dir
cp all_GenomeStats.csv ${WrkDir}/${GenomeSetName}_GenomeStats.csv
rm -R ${SCRATCH}/${GenomeSetName}/ref_genomes/
