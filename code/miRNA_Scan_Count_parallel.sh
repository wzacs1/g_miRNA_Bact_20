#!/bin/bash

#SBATCH --account=round
#SBATCH --partition=notchpeak-shared
#SBATCH -J par_miRNA
#SBATCH -n 24
#SBATCH -t 48:00:00
#SBATCH -D /uufs/chpc.utah.edu/common/home/u0210816/Projects/g_miRNA_Bact_20
#SBATCH -o /uufs/chpc.utah.edu/common/home/u0210816/Projects/g_miRNA_Bact_20/code/miRNA_Scan_Count_parallel_all_pangenomes.outerror

# Purpose:
# Find potential miRNA targets in sets of microbial genomes and return counts per genome and gene names.

# Future improvements to make:
# 1. User should be able to just provide full seed sequence, then extract subseqs (down to 6), put into array and loop through elements of array in search
# 2. Instead of seaprate script for each genome set (though this is faster), user should be able to provide list of genome sets to loop over. Could then incorporate analysis in single script as well.

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
NumProc=256
#######################################################

########### Setup #####################################
mkdir -p $SCRATCH
mkdir -p $WrkDir
# Module loading
module load fastx_toolkit
#######################################################


# Copy and extract genomes to search. (run through fastx_toolkit to convert fixed width fasta to unlimited line width. Change to tab separated instead of fasta makes search easier.)
# Currently for chocophlan pangenomes format.
mkdir -p ${SCRATCH}/${GenomeSetName}/ref_genomes

while read line
    do
    cp ${GenomesDir}/${line} ${SCRATCH}/${GenomeSetName}/ref_genomes/
    gunzip -c ${SCRATCH}/${GenomeSetName}/ref_genomes/${line} | fasta_formatter -t -o ${SCRATCH}/${GenomeSetName}/ref_genomes/${line%.gz}.tab
done < $Genomes

# Search genomes for target. Based on pangenomes format currently.
# Chocophlan pangenomes mantain the genome name in each seq identifier, but this is likely not the case for other ref datasets, thus genome name field.
cd ${SCRATCH}/${GenomeSetName}/ref_genomes

# Serial version for reference and clarity:
# for ingenome in *.tab
#     do
#     GenomeName=${ingenome%.centroids.v0.1.1.ffn.tab}
#     while read line
#         do
#         GeneID=`echo "$line" | cut -f 1`
#         GeneSeq=`echo "$line" | cut -f 2`
#         TargetACounts=`echo "$GeneSeq" | grep -o "$TargetSeqA" | wc -l`
#         TargetBCounts=`echo "$GeneSeq" | grep -o "$TargetSeqB" | wc -l`
#         TargetCCounts=`echo "$GeneSeq" | grep -o "$TargetSeqC" | wc -l`
#         TargetDCounts=`echo "$GeneSeq" | grep -o "$TargetSeqD" | wc -l`
#         TargetECounts=`echo "$GeneSeq" | grep -o "$TargetSeqE" | wc -l`
#         TargetFCounts=`echo "$GeneSeq" | grep -o "$TargetSeqF" | wc -l`
#         echo "$GenomeName,$GeneID,$TargetACounts,$TargetBCounts,$TargetCCounts,$TargetDCounts,$TargetECounts,$TargetFCounts" >> miRNA_hits.csv
#     done < $ingenome
# done

# Define the function. Write each genome to sep file since not sure if doing in parallel will mess up write order.
miRNASearch () {
    GenomeName=${1%.centroids.v0.1.1.ffn.tab}
    rm -f ${GenomeName}_miRNA_hits.csv
    echo "GenomeName,GeneID,${2}_counts,${3}_counts,${4}_counts,${5}_counts,${6}_counts,${7}_counts" > ${GenomeName}_miRNA_hits.csv
    while read line
        do
        GeneID=`echo "$line" | cut -f 1`
        GeneSeq=`echo "$line" | cut -f 2`
        TargetACounts=`echo "$GeneSeq" | grep -o "$2" | wc -l`
        TargetBCounts=`echo "$GeneSeq" | grep -o "$3" | wc -l`
        TargetCCounts=`echo "$GeneSeq" | grep -o "$4" | wc -l`
        TargetDCounts=`echo "$GeneSeq" | grep -o "$5" | wc -l`
        TargetECounts=`echo "$GeneSeq" | grep -o "$6" | wc -l`
        TargetFCounts=`echo "$GeneSeq" | grep -o "$7" | wc -l`
        echo "$GenomeName,$GeneID,$TargetACounts,$TargetBCounts,$TargetCCounts,$TargetDCounts,$TargetECounts,$TargetFCounts" >> ${GenomeName}_miRNA_hits.csv
    done < $1
}

export -f miRNASearch

# The base command for a single file would look like this (and does appear to work):
# miRNASearch <GenomeFile> "$TargetSeqA" "$TargetSeqB" "$TargetSeqC" "$TargetSeqD" "$TargetSeqE" "$TargetSeqF"

# Pass each

ls -1 --color=never *.tab | parallel -j $NumProc "miRNASearch {} $TargetSeqA $TargetSeqB $TargetSeqC $TargetSeqD $TargetSeqE $TargetSeqF"

# Collate results to one file:
echo "GenomeName,GeneID,${TargetSeqA}_counts,${TargetSeqB}_counts,${TargetSeqC}_counts,${TargetSeqD}_counts,${TargetSeqE}_counts,${TargetSeqF}_counts" > all_miRNA_hits.tmp
for f in *_miRNA_hits.csv
    do
    tail -n +2 $f >> all_miRNA_hits.tmp
done
mv all_miRNA_hits.tmp all_miRNA_hits.csv

# Cleanup, copy results to work dir
cp all_miRNA_hits.csv ${WrkDir}/${GenomeSetName}_miRNA_hits.csv
rm -R ${SCRATCH}/${GenomeSetName}/
