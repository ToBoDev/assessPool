## What is assessPool?
Pooling individual samples before sequencing (pool-seq) can greatly reduce cost while producing allele frequencies of single nucleotide polymorphisms (SNPs), but can be difficult to analyze. assessPool is a user-friendly R pipeline designed to analyze population structure from pool-seq data. Our software accepts a variant-call format (VCF) file and a FASTA-formatted reference, providing a straightforward transition from commonly used pipelines such as Stacks or dDocent. AssessPool can handle a variable number of pools and uses Popoolation2 to generate locus-by-locus pairwise FST values and associated Fisher T-test values as measures of population structure, and generates visuals to aid in assessing relationships between pools as well as coverage-based trends.

## What does this software do?
* Filters SNPs based on adjustable criteria using vcftools
* Determines pool number and prepares proper data structure for analysis
* Creates a customizable run script for Popoolation2 for all pairwise comparisons
* Runs Popoolation2 for each pair of pooled samples
* Generates population genetic statistics and plots for data visualization
