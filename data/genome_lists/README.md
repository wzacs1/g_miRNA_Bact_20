# Genome Subsets Acquisition Methods/Notes

## PATRIC genomes
- All genome lists from PATRIC were subsetted with web interface. Bacteria only. [https://www.patricbrc.org/view/Taxonomy/2](https://www.patricbrc.org/view/Taxonomy/2)
- Version 3.6.6.
- Subsets only include "Representative" genomes in "Reference Genome" field, then are subsetted as noted by "Host Name" field.

1. PATRIC_human-associated_RepresentativeGenomes.txt
2. PATRIC_non-mammalian-host_RepresentativeGenomes.txt
3. PATRIC_non-vertebrate-host_RepresentativeGenomes.txt

## Pangenomes / Chocophlan genomes
- Version 0.1.1 Chocophlan database pangenomes (part of Humann2)


## Disease subsets / Disbiome
- Subets from Disbiome database [https://disbiome.ugent.be/home](https://disbiome.ugent.be/home)
- Version unclear. Date of access = 2020-09-16
- Subsetted on Query term "Colorectal cancer" to get all "Organism" listed as changed. This includes both up and down in CRC and includes "Invasive colorectal cancer" as term.
- Resulting .json file is only export available. Poor data export limits ease of use of this database, despite tabular results.
 - disbiome_CRC_all.json
- Also just retrieved list of organisms associated with CRC (no values or up/down) as TSV from webpage display. This could be used akin to a "Gene Set".
 - disbiome_CRC_all_min.tsv
- Further filtering:
 - 

