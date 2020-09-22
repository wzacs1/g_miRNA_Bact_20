#!/bin/bash

#SBATCH --account=round
#SBATCH --partition=notchpeak-shared
#SBATCH -J PullScanGC
#SBATCH -n 24
#SBATCH -t 24:00:00
#SBATCH -D /uufs/chpc.utah.edu/common/home/u0210816/Projects/g_miRNA_Bact_20
#SBATCH -o /uufs/chpc.utah.edu/common/home/u0210816/Projects/g_miRNA_Bact_20/code/Pull_PATRIC_genomes_non-vertebrate-host.outerror

# Purpose:
# This script is written to allow functions to go into the original "miRNA_Scan_Count" scripts, which it should have already. Thus I maintain the initial setup and genome copy commands until new commands noted.

######### User-defined variables ######################
# Working directory to store results.  (Created if not present)
WrkDir=/uufs/chpc.utah.edu/common/home/u0210816/Projects/g_miRNA_Bact_20
GenomeSetName=patric_NonVertHost
# Scratch directory. (Created if not present). Potentially large amount of space required for genomes copy.
SCRATCH=/scratch/general/nfs1/u0210816/miRNA
# Target seq to search for. In future, iterate through single full sequence instead.
TargetSeqA=AGTTCTCA
TargetSeqB=GTTCTCA
TargetSeqC=TTCTCA
TargetSeqD=AGTTCT
TargetSeqE=AGTTCTC
TargetSeqF=GTTCTC
# Genome metadata file from PATRIC search with genome ID in first field, tsv with double quotes for each field
PatricGenomeMetadata=/uufs/chpc.utah.edu/common/home/u0210816/Projects/g_miRNA_Bact_20/data/genome_lists/PATRIC_non-vertebrate-host_RepresentativeGenomes.txt

# BOTH CREATED ON FLY when geiven Genome Metadata input. LEAVE EMPTY. Pangenomes files to search. Text file with one genome (gzipped) per line.
# Genomes=
# Location of reference genomes to search
# GenomesDir=
#######################################################

########### Setup #####################################
mkdir -p $SCRATCH
mkdir -p $WrkDir
# Module loading
module load fastx_toolkit
#######################################################

# Get List of Genomes To Retrieve and Set Names for Input into Scan Command
mkdir -p ${SCRATCH}/${GenomeSetName}/tmp_lists
mkdir -p ${SCRATCH}/${GenomeSetName}/ref_genomes
cd ${SCRATCH}/${GenomeSetName}/ref_genomes

cat $PatricGenomeMetadata | cut -f 1 | sed 's/"//g' | tail -n +2 > ${SCRATCH}/${GenomeSetName}/tmp_lists/gl.txt

while read line
    do
    echo "${line}.gz" >> ${SCRATCH}/${GenomeSetName}/tmp_lists/Genomes.txt
done < ${SCRATCH}/${GenomeSetName}/tmp_lists/gl.txt

Genomes=${SCRATCH}/${GenomeSetName}/tmp_lists/Genomes.txt
GenomesDir=${SCRATCH}/${GenomeSetName}/ref_genomes

# Pull from PATRIC
while read line
    do
    wget -O ${line}.ffn ftp://ftp.patricbrc.org/genomes/${line}/${line}.PATRIC.ffn
    gzip ${line}.ffn
done < ${SCRATCH}/${GenomeSetName}/tmp_lists/gl.txt

# Run miRNA Scan and GC stats

for f in *.gz
    do gunzip -c ${f} | fasta_formatter -t -o ${f%.ffn.gz}.tab
done
    
# First, collate the sequences together to form a temp genome file

miRNASearchGenome () {
    GenomeName=${1%.tab}
    echo "GenomeName,A_content,T_content,G_content,C_content,Size,${2}_counts,${3}_counts,${4}_counts,${5}_counts,${6}_counts,${7}_counts" > ${GenomeName}_miRNA_hits.csv
    cut -f 2 $1 | tr -d '\n' | tr "acgtuU" "ACGTTT" > tmp_${GenomeName}.txt
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
# rm -R ${SCRATCH}/${GenomeSetName}/


# Example PATRIC directory sturcure
# curl -l ftp://ftp.patricbrc.org/genomes/1000561.3/
# -rw-r--r--   1 p3       p3        2243968 Aug  4  2018 1000561.3.PATRIC.faa
# -rw-r--r--   1 p3       p3        1613039 Aug  4  2018 1000561.3.PATRIC.features.tab
# -rw-r--r--   1 p3       p3        5435958 Aug  4  2018 1000561.3.PATRIC.ffn
# -rw-r--r--   1 p3       p3          12033 Aug  4  2018 1000561.3.PATRIC.frn
# -rw-r--r--   1 p3       p3         839710 Aug  4  2018 1000561.3.PATRIC.gff
# -rw-r--r--   1 p3       p3         614155 Aug  4  2018 1000561.3.PATRIC.pathway.tab
# -rw-r--r--   1 p3       p3         204195 Aug  4  2018 1000561.3.PATRIC.spgene.tab
# -rw-r--r--   1 p3       p3         596882 Sep 20  2019 1000561.3.PATRIC.subsystem.tab
# -rw-r--r--   1 p3       p3         165350 Oct  4  2017 1000561.3.RefSeq.features.tab
# -rw-r--r--   1 p3       p3       11590473 Apr 21  2015 1000561.3.RefSeq.gbf
# -rw-r--r--   1 p3       p3          76418 Oct  4  2017 1000561.3.RefSeq.gff
# -rw-r--r--   1 p3       p3        6617253 Aug  4  2018 1000561.3.fna
