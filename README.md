# assessPool
Pool-Seq analysis software (in dev)

* Filters SNPs based on adjustable criterion with suggestions for pooled data 
* Determines pool number and prepares proper data structure for analysis
* Creates a customizable run script for Popoolation2 for all pairwise comparisons
* Runs Popoolation2
* Imports Popoolation2 output
* Generates population genetic statistics and plots for data visualization.

### Quick Start
1. To run assessPool, you will need: [samtools](http://www.htslib.org/), [cpan](https://www.cpan.org/), [vcftools](http://vcftools.sourceforge.net/) v0.1.12 or higher, [vcflib](https://github.com/vcflib/vcflib), and the R library [packrat](https://rstudio.github.io/packrat)
2. Navigate to desired working directory on your computer
3. Clone repo into directory: `git clone https://github.com/ToBoDev/assessPool.git`
4. Open assessPool Rproject in RStudio
5. Open assessPool_main.Rmd file in RStudio
6. Make sure you have your SNP file (.vcf) and reference file (.fasta) in the working directory
7. Note - the first time you run the script, library installation (via packrat) will take some time. 
8. Follow the steps in the main file to complete your analysis - run each code chunk and adjust parameters as needed 
  
 See the [wiki](https://github.com/ToBoDev/assessPool/wiki) for more information!
