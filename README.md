![assesspool_logo_sm](https://github.com/user-attachments/assets/cfa839d5-086b-434b-9880-9e69cb5cc265)

assessPool v2.0.0
Pool-Seq analysis software

* Filters SNPs based on adjustable criterion with suggestions for pooled data 
* Determines pool number and prepares proper data structure for analysis
* Creates a customizable run script for Popoolation2 for all pairwise comparisons
* Runs Popoolation2 and poolfstat
* Imports Popoolation2 and poolfstat output
* Generates population genetic statistics and plots for data visualization.

### Quick Start
1. To run assessPool, you will need: [R](https://www.r-project.org/) v3.6.0 or higher, [RStudio](https://www.rstudio.com/), [samtools](http://www.htslib.org/), [CPAN](https://www.cpan.org/), [vcftools](http://vcftools.sourceforge.net/) v0.1.12 or higher, [vcflib](https://github.com/vcflib/vcflib), and (optionally) [GNU parallel](https://www.gnu.org/software/parallel/)
2. Navigate to desired working directory on your computer
3. Clone repo into directory: `git clone https://github.com/ToBoDev/assessPool.git`
4. Open assessPool Rproject in RStudio
5. Open assessPool_main.Rmd file in RStudio
6. Make sure you have your SNP file (.vcf) and reference file (.fasta) in the working directory
7. Note - the first time you run the script, library installation (via {renv}) will take some time. 
8. Follow the steps in the main file to complete your analysis - run each code chunk and adjust parameters as needed 
  
 See the [wiki](https://github.com/ToBoDev/assessPool/wiki) for more information!
