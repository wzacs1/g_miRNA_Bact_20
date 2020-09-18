#!/bin/bash

#SBATCH --account=round-np
#SBATCH --partition=round-shared-np
#SBATCH -J ScanGC
#SBATCH -n 12
#SBATCH -t 24:00:00
#SBATCH -D /uufs/chpc.utah.edu/common/home/u0210816/Projects/g_miRNA_Bact_20
#SBATCH -o /uufs/chpc.utah.edu/common/home/u0210816/Projects/g_miRNA_Bact_20/code/Genome_miRNA_ScanAndGC.outerror

# Purpose:
# This script is written to allow functions to go into the original "miRNA_Scan_Count" scripts, which it should have already. Thus I maintain the initial setup and genome copy commands until new commands noted.

######### User-defined variables ######################
# Working directory to store results.  (Created if not present)
WrkDir=/uufs/chpc.utah.edu/common/home/u0210816/Projects/g_miRNA_Bact_20
GenomeSetName=all_pangenomes_wGenomeStats
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

miRNASearchGenome () {
    GenomeName=${1%.centroids.v0.1.1.ffn.tab}
    echo "GenomeName,A_content,T_content,G_content,C_content,Size,${2}_counts,${3}_counts,${4}_counts,${5}_counts,${6}_counts,${7}_counts" > ${GenomeName}_miRNA_hits.csv
    cut -f 2 $1 | tr -d '\n' > tmp_${GenomeName}.txt
    Size=`cat tmp_${GenomeName}.txt | wc -m`
    Acontent=`cat tmp_${GenomeName}.txt | tr -d "GCT" | wc -m`
    Gcontent=`cat tmp_${GenomeName}.txt | tr -d "TCA" | wc -m`
    Ccontent=`cat tmp_${GenomeName}.txt | tr -d "GAT" | wc -m`
    Tcontent=`cat tmp_${GenomeName}.txt | tr -d "GCA" | wc -m`
    TargetACounts=`cat tmp_${GenomeName}.txt | grep -o "$2" | wc -l`
    TargetBCounts=`cat tmp_${GenomeName}.txt | grep -o "$3" | wc -l`
    TargetCCounts=`cat tmp_${GenomeName}.txt | grep -o "$4" | wc -l`
    TargetDCounts=`cat tmp_${GenomeName}.txt | grep -o "$5" | wc -l`
    TargetECounts=`cat tmp_${GenomeName}.txt | grep -o "$6" | wc -l`
    TargetFCounts=`cat tmp_${GenomeName}.txt | grep -o "$7" | wc -l`
    echo "$GenomeName,$Acontent,$Gcontent,$Ccontent,$Tcontent,$Size,$TargetACounts,$TargetBCounts,$TargetCCounts,$TargetDCounts,$TargetECounts,$TargetFCounts" >> ${GenomeName}_miRNA_hits.csv
    rm tmp_${GenomeName}.txt
}

export -f miRNASearchGenome

ls -1 --color=never *.tab | parallel "miRNASearchGenome {} $TargetSeqA $TargetSeqB $TargetSeqC $TargetSeqD $TargetSeqE $TargetSeqF"

# Collate results to one file:
echo "GenomeName,A_content,T_content,G_content,C_content,Size,${TargetSeqA}_counts,${TargetSeqB}_counts,${TargetSeqC}_counts,${TargetSeqD}_counts,${TargetSeqE}_counts,${TargetSeqF}_counts" > all_miRNA_hits.tmp
for f in *_miRNA_hits.csv
    do
    tail -n +2 $f >> all_miRNA_hits.tmp
done
mv all_miRNA_hits.tmp all_miRNA_hits.csv

# Cleanup, copy results to work dir
cp all_miRNA_hits.csv ${WrkDir}/${GenomeSetName}_miRNA_hits.csv
rm -R ${SCRATCH}/${GenomeSetName}/

