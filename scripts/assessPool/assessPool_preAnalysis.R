#cessPOOL_preAnalysis.R

#VCF INFO
#NS=Number samples with data
#DP=Total read depth (all pops)
#NUMALT=Number of alternate alleles 
#AN=Allele Number (total, inc. ref)
#AF=Allele Freq
#LEN=Allele Length (1=SNP, 2+=MNP)
#PAIRED=Prop.of AltOBs that cam from properly paired reads
#PAIREDR=Prop. of REF.OBS that came from properly paired reads
#FORMAT
#DP=Read Depth per sample
#RO=REF obs. count
#AO=ALT obs. count
#AD=Number of each allele

#function to write given block of text to log file
write.log <- function(textin, fileout) {
  #if log file doesn't exist, create it
  if (file.exists(fileout)){
      write(textin, file=fileout, append=TRUE)
  } else{
      write(textin, file=fileout, append=FALSE)
  }    
}

#given a variable, pulls a subset out of saved vcf and returns it
subset.per.pool <- function(depth_var, as.vcf, POPS) {
  depth.tmp <- character()
  for(i in 1:length(POPS)) depth.tmp[i] <- c(paste(paste(depth_var,".",sep=""),POPS[i], sep=""))
  depth.df <- as.data.frame(geno(as.vcf)[depth_var])
  #cleanup dataframe
  depth.df <- depth.df[, -c(1,2)]
  depth.df$snpid <- row.names(depth.df)
  colnames(depth.df) <- c(depth.tmp,"snpid"); rownames(depth.df) <- NULL
  depth.df #returns subsetted dataframe
}

#given a character column from a dataframe, tidy up strings and return
clean.up.list <- function(data_in) {
  data_in <- gsub("c\\(", "", data_in)
  data_in <- gsub("\"", "", data_in)
  data_in <- gsub("\\)", "", data_in)
  data_in <- gsub(":", ",", data_in)
  data_in <- as.character(gsub(" ", "", data_in))
  data_in <- gsub("NA", "0", data_in)
  data_in #return cleaned data
}

#creates master dataframe (wide-form) and stacked dataframe (long-form) for SYNC file creation
preAnalysis <-function(working_dir, project_name, POPS, min.pool.number, min.depth.threshold, max.indel.length, include.multiallelic, include.indels, vcf_file, ref_file) {
  
    setwd(working_dir)
  
    #test to see if filtering occured - if so, use that vcf file
    if (file.exists(paste("filtered_", vcf_file, sep=""))){
      vcf_file <- (paste("filtered_", vcf_file, sep=""))
    }
  
    #creates directories for output and analysis
    #if directories already exist, will throw a warning message
    # To retrieve a subset of data from a VCF file, create a ScanVcfParam object. This object can specify genomic coordinates (ranges) or individual VCF elements to be extracted. When ranges are extracted, a tabix index file must exist for the VCF. See ?indexTabix for details.
    param <- ScanVcfParam(fixed = c("ALT","QUAL"), info = c("NS","DP","TYPE","LEN","NUMALT"), geno = c("DP","RO","AO","GT","GQ"))
    suppressWarnings(as.vcf <- readVcf(vcf_file,ref_file, param)); rm(param) ## Read-in data using VariantAnnotation
    project_name <- gsub(" ", "", project_name) 
  
    tryCatch({
      dir.create(paste(working_dir, project_name, sep="/"))
    }, warning=function(w){
      message("\nWarning: Project directory already exists - may overwrite files.")
    })
    setwd(paste(working_dir, project_name, sep="/"))
    suppressWarnings(dir.create(paste(working_dir, project_name, "output", sep="/")))
    suppressWarnings(dir.create(paste(working_dir, project_name, "popoolation", sep="/")))
    suppressWarnings(dir.create(paste(working_dir, project_name, "logs", sep="/")))
  
    #if filtering log exists, move it to new logs directory
    if (file.exists("filter.log")){
      file.copy(from="filter.log", to=paste(working_dir, project_name, "logs", "filter.log", sep="/"))
      file.remove("filter.log")
    }

    #if user did not define population names, get them from the VCF
    if (is.null(POPS)) { POPS <- samples(header(as.vcf)) }
    POPS <- gsub("_", "", POPS)
    
    #create population pairwise combination matrix (each column, represents a pair). If triplets desired, use `comb(pops, 3)`, quad `comb(pops, 4)`, etc
    popcomb <- combn(POPS,2) 
    ## create INFO df
    info.as <- as.data.frame(info(as.vcf))
    info.as$snpid <- row.names(info.as)
    colnames(info.as) <- c("NS","TDP","TYPE","LEN","AN","snpid")
    info.as$REF <- as.character(rowRanges(as.vcf)$REF)
    
    ## Coerce ALT (DNAstringsetlist) into comma separated list
    clist <- CharacterList(rowRanges(as.vcf)$ALT)
    mult <- elementNROWS(clist) > 1L #replaced elementLengths
    clist[mult] <- lapply(clist[mult], paste0, collapse=",")
    info.as$ALT <- unlist(clist)
    
    ## Add additional vectors
    info.as$Rlen <- with(info.as, nchar(as.character(REF)))
    info.as$Alen <- with(info.as, (nchar(as.character(info.as$ALT))-str_count(as.character(info.as$ALT), ","))/(1+str_count(as.character(info.as$ALT), ",")))
    info.as$INS.len <- info.as$Alen-info.as$Rlen 
    info.as$DEL.len <- info.as$Rlen-info.as$Alen 
    info.as$TYPE <- as.character(info.as$TYPE)
    info.as$LEN <- as.character(info.as$LEN)
    rm(clist)
    
    #Subset read depth (DP), reference depth (RO), alternate depth (AO), and genotype (GT) for each pool, combine into one dataframe
    #read depth, reference depth, alt depth, genotype subsets
    DP.as <- subset.per.pool("DP", as.vcf, POPS)
    RO.as <- subset.per.pool("RO", as.vcf, POPS)
    AO.as <- subset.per.pool("AO", as.vcf, POPS)

    # Merge summary dataframes by shared column (snpid) (will contain NAs)
    as <- Reduce(function(x, y) merge(x, y, by="snpid", all=TRUE), list(info.as,DP.as,RO.as,AO.as))
    rm(DP.as,RO.as,AO.as,info.as) #remove old temp dataframes
    
    #get chromosome id (CHROM) and position (POS) from snpid
    tst <- read.table(text = as$snpid, sep = ":", colClasses = "character")
    as$CHROM <- as.factor(tst$V1)
    tst2 <- read.table(text = tst$V2, sep = "_", colClasses = "character")
    as$POS <- as.integer(tst2$V1)
    rm(tst,tst2)
    
    #reorder columns and sort by chromosome and position
    as <- as[,c(1,which(names(as)=="CHROM"),which(names(as)=="POS"),2:(length(as)-2))]

    #if any population has missing data (e.g. all zeroes) remove them from analysis
    dp.cols <- grep("DP\\.", names(as))
    sums <- apply(as[,dp.cols], 2, function(x) sum(x, na.rm=T))
    if (any(as.numeric(sums) == 0)){
      invalid_pops <- data.frame(strsplit(names(which(sums == 0)), "\\."))[2,]
      for (i in invalid_pops){
        temp_inval <- as.character(i)
        temp.cols <- grep(temp_inval, names(as))
        
        #remove columns containing this population
        as <- as[,-temp.cols]
        POPS <- POPS[which(POPS!=temp_inval)]
        popcomb <- combn(POPS,2) 
        message("Removed ",temp_inval," pool due to lack of data.")
      }
    }     
                 
    #clean up list variables 
    #get ref depth and alternate depth columns, coerce to character
    ao.cols <- grep("AO\\.", names(as)); ro.cols <- grep("RO\\.", names(as))
    for(i in 1:length(ro.cols)) as[,ro.cols[i]] <- as.character(as[,ro.cols[i]]) 
    #clean up type and alternate allele counts
    for(i in 1:length(ao.cols)){
      as[,ao.cols[i]] <- as.character(as[,ao.cols[i]]) 
      as[,ao.cols[i]] <- clean.up.list(as[,ao.cols[i]])
    }
    as$TYPE <- clean.up.list(as.character(as$TYPE))
    
    #clean up snp labels
    as$TYPE <- gsub("snp", "bi-snp", as$TYPE)
    
    if (length(which(as$TYPE=="bi-snp,bi-snp" & as$AN==2)) != 0){
      as[which(as$TYPE=="bi-snp,bi-snp" & as$AN==2),]$TYPE <- "tri-snp"}
    if (length(which(as$TYPE=="bi-snp,bi-snp,bi-snp" & as$AN==3)) != 0){
      as[which(as$TYPE=="bi-snp,bi-snp,bi-snp" & as$AN==3),]$TYPE <- "quad-snp"}
    
    #keep only columns needed for SYNC file
    RO.tmp <- names(as)[ro.cols]; AO.tmp <- names(as)[ao.cols]
    as.tmp.cols <- c("snpid","CHROM","POS","REF","ALT","TYPE","AN","INS.len","DEL.len","Rlen","Alen",RO.tmp,AO.tmp)
    as.tmp <- as[,which(names(as) %in% as.tmp.cols)] 
    
    #Stack reference/alt allele counts by population 
    as.stack.ro <- melt(as.tmp, measure.vars=c(RO.tmp), variable_name="Pop") #convert from wide to long
    colnames(as.stack.ro)[which(names(as.stack.ro) == "value")] <- "RO"
    levels(as.stack.ro$Pop) <- POPS
    as.stack.ro <- as.stack.ro[-which(names(as.stack.ro) %in% AO.tmp)]
    
    as.stack.ao <- melt(as.tmp, measure.vars=AO.tmp, variable_name="Pop")
    colnames(as.stack.ao)[which(names(as.stack.ao) == "value")] <- "AO"
    levels(as.stack.ao$Pop) <- POPS
    as.stack.ao <- as.stack.ao[,which(names(as.stack.ao) %in% c("snpid","Pop","AO"))]
    
    #Merge together
    as.st <- merge(as.stack.ro, as.stack.ao, by=c("snpid","Pop"))
    as.st <- arrange(as.st, CHROM, POS) #sort by chromosome and position
    rm(as.stack.ro, as.stack.ao, as.tmp)
    #as.st <- as.st[,c(1,4,5,2,6,7,10,11,12,13,8,14,9,3)] #reorder columns
    as.st <- as.st[,c(1,3:6,9:12,7,13,8,14,2)]
    
    #coerce to data types
    as.st$RO <- as.integer(as.character(as.st$RO))
    as.st$AO <- as.character(as.st$AO)
    as.st$Rlen <- as.integer(as.character(as.st$Rlen))
    as.st$Alen <- as.integer(as.character(as.st$Alen))
    as.st$TYPE <- as.factor(as.character(as.st$TYPE)) #coerce to factor
    as.st[is.na(as.st)] <- 0 #remove NA
    
    return(list("as"=as, "as.st"=as.st, "POPS"=POPS, "popcomb"=popcomb, "project_name"=project_name))

}
