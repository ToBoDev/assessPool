#sets up new filtering file and outputs initial data
vcf_init <- function(vcf_file, wd, output){
  setwd(working_dir)
  message(paste(snp_count(vcf_file, working_dir), "SNPs before filtering.", sep=" "))
  
  if (output){
    system(paste(paste(working_dir,"scripts/filtering/ErrorCount.sh",sep="/"), paste(working_dir, vcf_file, sep="/"), sep=" "))
  }
  
  #copies vcf file
  filtered_filename <- paste("filtered_",vcf_file,sep="")
  sink("/dev/null"); file.copy(from=vcf_file, to=filtered_filename, overwrite = TRUE); sink()
  
  filtering_steps <- data.frame(Filter=NA, SNPs=NA)
  return(filtering_steps)
}

#counts SNPs in given vcf file
snp_count <- function(vcf, wd) {
  system(paste("mawk '!/#/'", paste(wd, as.character(vcf), sep="/"), "| wc -l", sep=" "), intern=T)
}

#calculates differences in # SNPs between two vcf files
snp_diff <- function(vcf1, wd, vcf2) {
  message(
    snp_count(vcf2, wd),
    " SNPs kept. ",
    as.numeric(snp_count(vcf1, wd)) - as.numeric(snp_count(vcf2, wd)),
    " (",
    round(100 - (as.numeric(snp_count(vcf2, wd)) / as.numeric(snp_count(vcf1, wd))) * 100, 2),
    "%) filtered. Cumulative filtered: ",
    as.numeric(snp_count(vcf_file, wd)) - as.numeric(snp_count(vcf2, wd)),
    " (",
    round(100 - (as.numeric(snp_count(vcf2, wd)) / as.numeric(snp_count(vcf_file, wd))) * 100, 2),
    "%)"
  )
}

### Filter by number of pools
filter_numpools <- function(working_dir, project_name, vcf_file, x, fdf, out){
  setwd(working_dir)
  system(paste(paste(vcflib_PATH, 'vcffilter', sep=""), '-f "NS >', x, '| NS =', x, '"', paste("filtered_",vcf_file,sep=""), ">", paste("temp_",vcf_file,sep="")))
  
  message("After filtering by minimum number of pools (",x,"):")
  snp_diff(paste("filtered_",vcf_file,sep=""), working_dir, paste("temp_",vcf_file,sep=""))
  if (out){
    system(paste(paste(working_dir,"/scripts/filtering/ErrorCount.sh", sep=""), paste("temp_",vcf_file,sep=""), sep=" "))
  }
    
  #rewrite temp folder
  sink("/dev/null"); file.copy(from=paste("temp_",vcf_file,sep=""), to=paste("filtered_",vcf_file,sep=""), overwrite = TRUE); file.remove(paste("temp_",vcf_file,sep="")); sink()
  message("Wrote filtered SNPs vcf to filtered_",vcf_file)
  
  #Save # SNPs to filtering log
  numSNPs <- system(paste("mawk '!/#/'", paste("filtered_",vcf_file,sep=""), "| wc -l", sep=" "), intern=T)
  fdf <- na.omit(rbind(fdf, c(paste("MinPools",x,sep="_"), as.numeric(numSNPs))))
  return(fdf)
}

### Filter by quality 
filter_quality <- function(working_dir, project_name, vcf_file, x, fdf, out){
  setwd(working_dir)
  system(paste(vcftools_path_export,"vcftools --vcf", paste("filtered_",vcf_file,sep=""), "--minQ", x, "--recode --recode-INFO-all --out temp"))
  
  message("")
  message("After filtering by minimum quality (",x,"):")
  snp_diff(paste("filtered_",vcf_file,sep=""), working_dir, paste("temp.recode.vcf",sep=""))
  if (out){
    system(paste(paste(working_dir,"/scripts/filtering/ErrorCount.sh", sep=""), paste("temp.recode.vcf"), sep=" "))
  }
  
  #rewrite temp folder
  sink("/dev/null"); file.copy(from=paste("temp.recode.vcf",sep=""), to=paste("filtered_",vcf_file,sep=""), overwrite = TRUE); file.remove(paste("temp.recode.vcf",sep="")); sink()
  message("Wrote filtered SNPs vcf to filtered_",vcf_file)
  
  #Save # SNPs to filtering log
  numSNPs <- system(paste("mawk '!/#/'", paste("filtered_",vcf_file,sep=""), "| wc -l", sep=" "), intern=T)
  fdf <- na.omit(rbind(fdf, c(paste("MinQuality",x,sep="_"), as.numeric(numSNPs))))
  return(fdf)
}

### Filter by minimum depth
filter_mindepth <- function(working_dir, project_name, vcf_file, x, y, fdf, out){
  setwd(working_dir)
  system(paste(vcftools_path_export,"vcftools --vcf", paste("filtered_",vcf_file,sep=""), "--minDP", x, "--remove-filtered-geno-all --max-missing-count", y, "--recode --recode-INFO-all --out temp"))
  
  message("")
  message("After filtering by minimum depth (",x,"):")
  snp_diff(paste("filtered_",vcf_file,sep=""), working_dir, paste("temp.recode.vcf",sep=""))
  if (out){
    system(paste(paste(working_dir,"/scripts/filtering/ErrorCount.sh", sep=""), "temp.recode.vcf", sep=" "))
  }
  
  #rewrite temp folder
  sink("/dev/null"); file.copy(from=paste("temp.recode.vcf",sep=""), to=paste("filtered_",vcf_file,sep=""), overwrite = TRUE); file.remove(paste("temp.recode.vcf",sep="")); sink()
  message("Wrote filtered SNPs vcf to filtered_",vcf_file)
  
  #Save # SNPs to filtering log
  numSNPs <- system(paste("mawk '!/#/'", paste("filtered_",vcf_file,sep=""), "| wc -l", sep=" "), intern=T)
  fdf <- na.omit(rbind(fdf, c(paste("MinDepth",x,sep="_"), as.numeric(numSNPs))))
  return(fdf)
}

### Filter by max allele length
filter_maxallelelength <- function(working_dir, project_name, vcf_file, x, fdf, out){
  setwd(working_dir)
  system(paste(paste(vcflib_PATH, 'vcffilter', sep=""), '-f "LEN <', x, '| LEN =', x, '"', paste("filtered_",vcf_file,sep=""), ">", paste("temp_",vcf_file,sep="")))
  
  message("After filtering by maximum allele length (",x,"):")
  snp_diff(paste("filtered_",vcf_file,sep=""), working_dir, paste("temp_",vcf_file,sep=""))
  if (out){
    system(paste(paste(working_dir,"/scripts/filtering/ErrorCount.sh", sep=""), paste("temp_",vcf_file,sep=""), sep=" "))
  }
  
  #rewrite temp folder
  sink("/dev/null"); file.copy(from=paste("temp_",vcf_file,sep=""), to=paste("filtered_",vcf_file,sep=""), overwrite = TRUE); file.remove(paste("temp_",vcf_file,sep="")); sink()
  message("Wrote filtered SNPs vcf to filtered_",vcf_file)
  
  #Save # SNPs to filtering log
  numSNPs <- system(paste("mawk '!/#/'", paste("filtered_",vcf_file,sep=""), "| wc -l", sep=" "), intern=T)
  fdf <- na.omit(rbind(fdf, c(paste("MaxAlleleLength",x,sep="_"), as.numeric(numSNPs))))
  return(fdf)
}

### Filter by mispaired reads
filter_mispaired <- function(working_dir, project_name, vcf_file, fdf, out){
  setwd(working_dir)
  system(paste(paste(vcflib_PATH, 'vcffilter', sep=""), '-f "PAIRED > 0.05 & PAIREDR > 0.05 & PAIREDR / PAIRED < 1.75 & PAIREDR / PAIRED > 0.25 | PAIRED < 0.05 & PAIREDR < 0.05" -s', paste("filtered_",vcf_file,sep=""), ">", paste("temp_",vcf_file,sep="")))
  
  message("After filtering for mispaired reads:")
  snp_diff(paste("filtered_",vcf_file,sep=""), working_dir, paste("temp_",vcf_file,sep=""))
  if (out){
    system(paste(paste(working_dir,"/scripts/filtering/ErrorCount.sh", sep=""), paste("temp_",vcf_file,sep=""), sep=" "))
  }
  
  #rewrite temp folder
  sink("/dev/null"); file.copy(from=paste("temp_",vcf_file,sep=""), to=paste("filtered_",vcf_file,sep=""), overwrite = TRUE); file.remove(paste("temp_",vcf_file,sep="")); sink()
  message("Wrote filtered SNPs vcf to filtered_",vcf_file)
  
  #Save # SNPs to filtering log
  numSNPs <- system(paste("mawk '!/#/'", paste("filtered_",vcf_file,sep=""), "| wc -l", sep=" "), intern=T)
  fdf <- na.omit(rbind(fdf, c("MispairedReads", as.numeric(numSNPs))))
  return(fdf)
}

### Filter by quality/depth ratio (helps remove low-quality, high-depth reads)
filter_qualdepth <- function(working_dir, project_name, vcf_file, x, fdf, out){
  setwd(working_dir)
  system(paste(paste(vcflib_PATH, 'vcffilter', sep=""), '-f "QUAL / DP >', x, '| QUAL / DP =', x, '"', paste("filtered_",vcf_file,sep=""), ">", paste("temp_",vcf_file,sep="")))
  
  message("After filtering by minimum quality/depth ratio: (",x,"):")
  snp_diff(paste("filtered_",vcf_file,sep=""), working_dir, paste("temp_",vcf_file,sep=""))
  if (out){
    system(paste(paste(working_dir,"/scripts/filtering/ErrorCount.sh", sep=""), paste("temp_",vcf_file,sep=""), sep=" "))
  }
  
  #rewrite temp folder
  sink("/dev/null"); file.copy(from=paste("temp_",vcf_file,sep=""), to=paste("filtered_",vcf_file,sep=""), overwrite = TRUE); file.remove(paste("temp_",vcf_file,sep="")); sink()
  message("Wrote filtered SNPs vcf to filtered_",vcf_file)
  
  #Save # SNPs to filtering log
  numSNPs <- system(paste("mawk '!/#/'", paste("filtered_",vcf_file,sep=""), "| wc -l", sep=" "), intern=T)
  fdf <- na.omit(rbind(fdf, c(paste("QualDepthRatio",x,sep="_"), as.numeric(numSNPs))))
  return(fdf)
}

### dDocent filtering based on mean depth per site vs. quality score
#helps filter out true variants vs false variants
filter_ddocent <- function(working_dir, project_name, vcf_file, fdf, out){
  setwd(working_dir)
  
  system(paste("cut -f8", paste("filtered_",vcf_file,sep=""), "| grep -oe ","DP=[0-9]*"," | sed -s 's/DP=//g' > F5.DEPTH")) #pulls out read depth per SNP
  system(paste("mawk '!/#/'", paste("filtered_",vcf_file,sep=""), "| cut -f1,2,6 > F5.loci.qual")) #pulls out snpid (CHROM/POS) and quality score per SNP
  dp1 <- system("mawk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' F5.DEPTH", intern=T)
  dp2 <- system(paste('python -c "print int(', dp1,'+3*(', dp1,'**0.5))"', sep=""), intern=T)
  system(paste("paste F5.loci.qual F5.DEPTH | mawk -v x=",dp2," '$4 > x' | mawk '$3 < 2 * $4' > F5.lowQDloci", sep=""))
  system(paste(vcftools_path_export, ' vcftools --vcf', paste("filtered_",vcf_file,sep=""), '--recode-INFO-all --out temp --exclude-positions F5.lowQDloci --recode'))
  
  message("")
  message("After filtering by quality score to mean depth to depth ratio:")
  snp_diff(paste("filtered_",vcf_file,sep=""), working_dir, paste("temp.recode.vcf",sep=""))
  if (out){
    system(paste(paste(working_dir,"/scripts/filtering/ErrorCount.sh", sep=""), paste("temp.recode.vcf"), sep=" "))
  }
  
  #rewrite temp folder
  sink("/dev/null"); file.copy(from=paste("temp.recode.vcf",sep=""), to=paste("filtered_",vcf_file,sep=""), overwrite = TRUE); 
  file.remove(c("F5.loci.qual", "F5.lowQDloci", "temp.recode.vcf")); sink()
  message("Wrote filtered SNPs vcf to filtered_",vcf_file)
  
  #Save # SNPs to filtering log
  numSNPs <- system(paste("mawk '!/#/'", paste("filtered_",vcf_file,sep=""), "| wc -l", sep=" "), intern=T)
  fdf <- na.omit(rbind(fdf, c("QualMeanDepthRatio", as.numeric(numSNPs))))
  return(fdf)
}

### Filter by maximum mean depth
filter_maxmeandepth <- function(working_dir, project_name, vcf_file, x, fdf, out){
  setwd(working_dir)
  system(paste(vcftools_path_export,"vcftools --vcf", paste("filtered_",vcf_file,sep=""), "--max-meanDP", x, "--recode --recode-INFO-all --out temp"))
  
  message("")
  message("After filtering by maximum mean depth (",x,"):")
  snp_diff(paste("filtered_",vcf_file,sep=""), working_dir, paste("temp.recode.vcf",sep=""))
  if (out){
    system(paste(paste(working_dir,"/scripts/filtering/ErrorCount.sh", sep=""), paste("temp.recode.vcf"), sep=" "))
  }
  
  #rewrite temp folder
  sink("/dev/null"); file.copy(from=paste("temp.recode.vcf",sep=""), to=paste("filtered_",vcf_file,sep=""), overwrite = TRUE); file.remove(paste("temp.recode.vcf",sep="")); sink()
  message("Wrote filtered SNPs vcf to filtered_",vcf_file)
  system(paste("cut -f8", paste("filtered_",vcf_file,sep=""), "| grep -oe ","DP=[0-9]*"," | sed -s 's/DP=//g' > F5.DEPTH"))
  
  #Save # SNPs to filtering log
  numSNPs <- system(paste("mawk '!/#/'", paste("filtered_",vcf_file,sep=""), "| wc -l", sep=" "), intern=T)
  fdf <- na.omit(rbind(fdf, c(paste("MaxMeanDepth",x,sep="_"), as.numeric(numSNPs))))
  return(fdf)
}

