## Aim/Questions
1. Scan bacterial genomes for miRNA146 potential targets. 
2. Count how prevalent targets are in each genome and to what genes they are connected.
3. Are the targets more abundant than would be expected from random observation of sequences.
4. Do host-associated bacterial genomes have more targets than non-host-associated.

## Method (# corresponds to above Aim/Questions)
1. Genome Scan
	a. Scan pangenomes.
	b. Use simple regular expression pattern matching. Start by breaking each potential target up (random distribution more complicated otherwise?)
2. Counting
	a. Count with grep built-in or line count. Counts by gene ID and by genome; 2 separate fields
3. Normalization
	a. Determine GC content of each pangenome using `tr`
4. Enrichment
	a. Use neg binomial or hypergeometric distribution to test if observed is different from expecrted.

## Potential issues
1. Genome Scan
	a. Pangenomes can aggregate non-host-associated and host-associated to a single genome. Refseq would keep strains separate, but has redundancy among species genomes (eg. hundreds of E. coli genomes), and don't have a way to partition into host-associated (human wgs not run on RefSeq)
		- Probably need to stick with pangenomes for now. Parsing to host and non-host associated to slow/manual. 
		- Choosing a strictly non-host-associated group of bacteria to compare (acid mine drainage, taxa never/rarely observed with humans, etc.) 
2. Counting
3. Normalization
4. Enrichment
