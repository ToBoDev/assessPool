getFasta <- function(df, ref){
  df <- unique(data.frame(CHROM=df$CHROM))
  return(merge(df, ref, by="CHROM", all.x=T, all.y=F))
}

getPopNames <- function(popsin){
  popsin <- lapply(unique(unlist(lapply(popsin, function(x) strsplit(x=gsub(paste(project_name,"_", sep=""), "", gsub(".fst", "", x)), split="_")), use.names = FALSE)), function(x) paste("DP", x, sep="."))
  return(popsin)
}

popNACheck <- function(row.in, tmp.cov){
  tmp.pops <- unlist(row.in["popIncl"], use.names=FALSE)
  tmp.inds <- row.in[names(row.in) %in% tmp.pops]
  bool.out <- all(na.omit(tmp.inds) >= tmp.cov)
  return(bool.out)
}

#define new sumlog function to deal with NAs and zeros 
postPopoolation <- function(filetype, project_name, as, popcomb, strong_diff, p_cutoff, fasta_generation, ref_file, first_mincov, last_mincov, cov_step){

    file.copy(from=vcf_file, to=paste(working_dir,project_name,vcf_file, sep="/"), overwrite = TRUE)
    file.copy(from=ref_file, to=paste(working_dir,project_name,ref_file, sep="/"), overwrite = TRUE)  
    setwd(paste(working_dir,"/",project_name,"/popoolation", sep=""))
  
    #add parameters to logfile
    write.log("Analysis Parameters:", paste(working_dir, project_name, "logs/analysis.log", sep="/"))
    write.log(paste("Strongly differentiated FST >=",strong_diff), paste(working_dir, project_name, "logs/analysis.log", sep="/"))
    write.log(paste("Significant FET <=", p_cutoff), paste(working_dir, project_name, "logs/analysis.log", sep="/"))
  
    C <- ncol(popcomb) #number of combinations
    combs <- list() #empty list for pairwise tags
    for(i in 1:C) combs[i] <- paste(project_name, popcomb[1,i], popcomb[2,i], sep='_') #pairwise tags
    df_names <- list(); df_index <- 1
    
    #import PoPoolation2 output, save as dataframes
    for(ft in filetype){
      for(i in 1:C){
        f=paste(combs[i], ft,sep=''); assign(f, read.table(f, sep='\t')) #read in dataframe
        d=get(f) #get dataframe object
        colnames(d)[c(1,2,6)]=c("CHROM","value"); d$pair <- combs[i]; d$analysis <- ft #reorganize columns
        d$value <- as.numeric(gsub(".*=","",as.character(d$value))) #remove unwanted characters, convert to numeric
        d$snpid <- paste(d[,1], d[,2], sep = '_') #add snpid column
        d <- d[,c(1:2,6:9)] #drop unneeded columns
        if (ft==".fet"){ d$value <- 10^(-(d$value)) } #convert log(p) to p-value
        assign(f,d)
        df_names[df_index] <- paste(combs[i],ft,sep=""); df_index <- df_index+1
      }
    }
    
    #save dataframe in long format for future data manipulation
    postpop.master.long <- bind_rows(mget(unlist(df_names)))
    for(ft in filetype) rm(list=ls(pattern=paste("*.",ft,sep=''))); rm(d)  
  
    #convert to wide format for export
    postpop.master.tmp <- postpop.master.long
    postpop.master.tmp$tag <- paste(postpop.master.tmp$pair, postpop.master.tmp$analysis,sep="")
    postpop.master.tmp$pair <- NULL; postpop.master.tmp$analysis <- NULL
    postpop.master.wide <- spread(postpop.master.tmp, tag, value)
    rm(postpop.master.tmp)
    
    #add columns from master dataframe
    as$snpid <- paste(as$CHROM, as$POS, sep="_")
    postpop.master.wide <- merge(as[,-which(names(as) %in% c("LEN","INS.len","DEL.len"))], postpop.master.wide[,-which(names(postpop.master.wide) %in% c("CHROM"))], by="snpid", all.y=TRUE, all.x=FALSE) 
    
    #filter out monomorphic SNPs with all pairwise Fst=0
    if (".fst" %in% filetype){
      postpop.master.wide$fst0s <- apply(postpop.master.wide[,grep(".fst", names(postpop.master.wide))], 1, function(x) length(which(x==0)))  #counts instances of Fst=0
      postpop.master.wide <- postpop.master.wide[which(postpop.master.wide$fst0s < C),]; postpop.master.wide$fst0s <- NULL
    }
                                         
    #keep track of which populations have a SNP at this position (slow...)
    postpop.master.wide$popIncl <- apply(postpop.master.wide[,grep(filetype[1],names(postpop.master.wide))], 1, function(x) names(postpop.master.wide[,grep(filetype[1],names(postpop.master.wide))])[which(!is.na(x))])
    postpop.master.wide$popIncl <- apply(postpop.master.wide[,"popIncl",drop=F], 1, getPopNames)
    
    #count number of populations in which a SNP is not called; should not be any NAs but to double check
    postpop.master.wide$fstNAs <- apply(postpop.master.wide[,grep(filetype[1],names(postpop.master.wide))], 1, function(x) sum(is.na(x)))
    postpop.master.wide <- postpop.master.wide[!(postpop.master.wide$fstNAs >= C),] # remove rows without SNPS called in any pair (i.e. those that did not pass through filters)
    m <- paste(nrow(postpop.master.wide), "informative variable sites retained after PoPoolation2.")
    message(m); write.log(m, paste(working_dir, project_name, "logs/analysis.log", sep="/"))
    
    #apply same filters to data stored in long format
    snpids_filter <- data.frame("snpid"=postpop.master.wide$snpid)
    postpop.master.long <- merge(postpop.master.long, snpids_filter, all=FALSE)
    
    #clean up environment and organize data
    postpop.master.wide$TYPE <- as.factor(postpop.master.wide$TYPE); postpop.master.long$CHROM <- as.factor(postpop.master.long$CHROM)
    postpop.master.wide <- arrange(postpop.master.wide, CHROM, POS) #sort by CHROM/POS
    
    ####################### File export
    
    #export all variable sites
    popl.mast.export <- data.frame(lapply(postpop.master.wide, as.character), stringsAsFactors=FALSE)
    write.csv(popl.mast.export, paste(working_dir,"/", project_name, "/output/", project_name,"_allvar_postPopoolation2.csv", sep=""), row.names=FALSE)
    m <- paste("Exported all variable sites to ","output/", project_name,"_allvar_postPopoolation2.csv", sep="")
    message(m); write.log(m, paste(working_dir, project_name, "logs/analysis.log", sep="/"))
    
    #Test how many SNPs are called in all populations
    called_allpops <- postpop.master.wide[postpop.master.wide$fstNAs==0,] 
    m <- paste(nrow(called_allpops), "SNPs called in all pools.")
    message(m); write.log(m, paste(working_dir, project_name, "logs/analysis.log", sep="/"))
    
    #export sites called in all populations
    popl.allpops.exports <- data.frame(lapply(called_allpops, as.character), stringsAsFactors=FALSE)
    write.csv(popl.mast.export, paste(working_dir,"/", project_name,"/output/", project_name, "_calledAllPools_postPopoolation2.csv", sep=""), row.names=FALSE)
    m <- paste("Exported variable sites called in all pools to ","output/", project_name,"_calledAllPools_postPopoolation2.csv", sep="")
    message(m); write.log(m, paste(working_dir, project_name, "logs/analysis.log", sep="/"))
    
    if (".fst" %in% filetype){
      
      #extract sites that show strong differentiation
      postpop.master.wide$DiffFST <- apply(postpop.master.wide[,grep(".fst", names(postpop.master.wide))], 1, function(x) any(x >= strong_diff))
      popl.appfx <- postpop.master.wide[which(postpop.master.wide$DiffFST==TRUE),]
      #popl.appfx <- na.omit(postpop.master.wide[apply(postpop.master.wide[,grep(".fst", names(postpop.master.wide))], 1, function(x) any(x >= strong_diff)),])
      m <- paste(nrow(popl.appfx), " SNPs are strongly differentiated (FST>=", strong_diff, ") in at least one comparison.", sep="")
      message(m); write.log(m, paste(working_dir, project_name, "logs/analysis.log", sep="/"))
      
      write.csv(popl.appfx, paste(working_dir,"/", project_name,"/output/", project_name, "_strongly_differentiated_sites.csv", sep=""), row.names=FALSE)
      m <- paste("Exported strongly differentiated sites to ","output/", project_name,"_strongly_differentiated_sites.csv", sep="")
      message(m); write.log(m, paste(working_dir, project_name, "logs/analysis.log", sep="/"))
      
      #extract alternatively fixed sites 
      postpop.master.wide$FixedFST <- apply(postpop.master.wide[,grep(".fst", names(postpop.master.wide))], 1, function(x) any(x == 1.0))
      popl.fixed <- postpop.master.wide[which(postpop.master.wide$FixedFST==TRUE),]
      #popl.fixed <- na.omit(postpop.master.wide[apply(postpop.master.wide[,grep(".fst", names(postpop.master.wide))], 1, function(x) any(x == 1.0)),])
      m <- paste(nrow(popl.fixed), "sites are alternatively fixed (FST=1) in an least one comparison.")
      message(m); write.log(m, paste(working_dir, project_name, "logs/analysis.log", sep="/"))
      
      write.csv(popl.fixed, paste(working_dir,"/", project_name,"/output/", project_name, "_fixed_snps.csv", sep=""), row.names=FALSE)
      m <- paste("Exported alternatively fixed sites to ","output/", project_name,"_fixed_snps.csv", sep="")
      message(m); write.log(m, paste(working_dir, project_name, "logs/analysis.log", sep="/"))
      
      if (fasta_generation==TRUE){
        ref_file_df <- read.fasta(paste(working_dir, ref_file, sep="/"))
        ref_file_df <- data.frame(CHROM=names(ref_file_df), Seq=toupper(unlist(getSequence(ref_file_df, as.string=T))))
        
        #generate FASTA from strongly differentiated sites
        sd_fasta <- getFasta(popl.appfx, ref_file_df)
        write.fasta(sequences=as.list(sd_fasta$Seq), names=sd_fasta$CHROM, file.out=paste(working_dir,"/", project_name,"/output/", project_name, "_strongly_differentiated_sites_fasta.fa", sep=""), nbchar=1000)
        m <- paste("Wrote ", nrow(sd_fasta), " contigs to /output/", project_name, "_strongly_differentiated_sites_fasta.fa", sep="")
        message(m); write.log(m, paste(working_dir, project_name, "logs/analysis.log", sep="/"))
        
        #generate FASTA from alternatively fixed sites
        alt_fasta <- getFasta(popl.fixed, ref_file_df)
        write.fasta(sequences=as.list(alt_fasta$Seq), names=alt_fasta$CHROM, file.out=paste(working_dir,"/", project_name,"/output/", project_name, "_fixed_sites_fasta.fa", sep=""), nbchar=1000)
        m <- paste("Wrote ", nrow(alt_fasta), " contigs to /output/", project_name, "_fixed_sites_fasta.fa", sep="")
        message(m); write.log(m, paste(working_dir, project_name, "logs/analysis.log", sep="/"))
      }
    }
    
    if (".fet" %in% filetype){
      postpop.master.wide$LowP <- apply(postpop.master.wide[,grep(".fet", names(postpop.master.wide))], 1, function(x) any(x <= p_cutoff))
      popl.lowP <- postpop.master.wide[which(postpop.master.wide$LowP==TRUE),]
      #popl.fixed <- na.omit(postpop.master.wide[apply(postpop.master.wide[,grep(".fst", names(postpop.master.wide))], 1, function(x) any(x == 1.0)),])
      m <- paste(nrow(popl.lowP), " SNPs have a low p-value (p<=", p_cutoff, ") in at least one comparison.", sep="")
      message(m); write.log(m, paste(working_dir, project_name, "logs/analysis.log", sep="/"))
      
      write.csv(popl.lowP, paste(working_dir,"/", project_name,"/output/", project_name, "_lowP_snps.csv", sep=""), row.names=FALSE)
      m <- paste("Exported low p-value sites to ","output/", project_name,"_lowP_snps.csv", sep="")
      message(m); write.log(m, paste(working_dir, project_name, "logs/analysis.log", sep="/"))
    }
    
    ###############add min coverage level to dataframe
    
    covs <- rev(c(seq(first_mincov,last_mincov,cov_step))) 
    #convert to numeric and replace NAs with 0
    #postpop.master.wide[,grep("DP\\.", names(postpop.master.wide))] <- apply(postpop.master.wide[,grep("DP\\.", names(postpop.master.wide))], 1, as.character)
    #postpop.master.wide[,grep("DP\\.", names(postpop.master.wide))] <- apply(postpop.master.wide[,grep("DP\\.", names(postpop.master.wide))], 1, as.numeric)
    postpop.master.wide[,grep("DP\\.", names(postpop.master.wide))][is.na(postpop.master.wide[,grep("DP\\.", names(postpop.master.wide))])] <- 0
    df_names_total <- list(); df_index <- 1
    
    message("")
    #assign min in-population coverage levels to SNPs

    tmp.idx <- c() #initialize list to keep track of rows already analyzed
    d <- data.frame()
    tmp.postpop <- postpop.master.wide[,-grep("RO\\.|AO\\.", names(postpop.master.wide))]
    tmp.postpop <- tmp.postpop[,-which(names(tmp.postpop) %in% c("CHROM","POS","NS","TDP","TYPE","AN","REF","ALT","Rlen","Alen"))]
    #remove uneeded columns
    rm(postpop.master.wide) #no longer needed
    gc()

    #get SNPs at each minimum coverage level
    for(i in 1:length(covs)){
      df_name <- paste("popl.master.", covs[i],"x", sep="")
      if (length(tmp.idx>0)){
        tmp.postpop <- tmp.postpop[-tmp.idx,]
      }
      bool.out <- apply(tmp.postpop[,c("popIncl", names(tmp.postpop)[grep("DP\\.", names(tmp.postpop))])], 1, function(x) popNACheck(x, covs[i]))
      tmp.idx <- as.numeric(which(bool.out))
      tmp.d <- tmp.postpop[tmp.idx,]
      tmp.d <- tmp.d[,-grep("DP\\.", names(tmp.d))]

      if (nrow(tmp.d)>0){
        if (nrow(d)>0){
          #new SNPs and old SNPs
          #add from previous step
          d$MinCoverage <- covs[i]
          tmp.d$MinCoverage <- covs[i]
          d <- rbind(d, tmp.d); rm(tmp.d)
        } else{
          #new SNPs but no old SNPs
          tmp.d$MinCoverage <- covs[i]
          d <- tmp.d; rm(tmp.d)
        }

        #save new dataframe
        assign(df_name, d)
        df_names_total[df_index] <- paste(df_name); df_index <- df_index+1
        message(paste(nrow(d)," SNPs at ",covs[i],"x coverage.", sep=""))

      } else{
        if (nrow(d)>0){
          #old SNPs but no new SNPs
          d$MinCoverage <- covs[i]
          assign(df_name, d)
          df_names_total[df_index] <- paste(df_name); df_index <- df_index+1
          message(paste(nrow(d)," SNPs at ",covs[i],"x coverage.", sep=""))

        } else{
          assign(df_name, NA)
          df_names_total[df_index] <- paste(df_name); df_index <- df_index+1
          message(paste("0 SNPs at ",covs[i],"x coverage.", sep=""))
        }
      }
    }
    
    #append dataframes to build master dataframes (wide and long) by coverage
    postpop.master.wide.bycov <- unique(bind_rows(mget(unlist(df_names_total))))
    min_covs <- data.frame(snpid=postpop.master.wide.bycov$snpid, MinCoverage=postpop.master.wide.bycov$MinCoverage, fstNAs=postpop.master.wide.bycov$fstNAs)
    postpop.master.long.bycov <- data.frame(merge(postpop.master.long, min_covs, all=FALSE)); rm(min_covs)
    postpop.master.long.bycov$pair <- as.character(postpop.master.long.bycov$pair)
    postpop.master.long.bycov$CHROM <- as.character(postpop.master.long.bycov$CHROM)
    
    summarizePopoolation <- function(data_in){
        #aggregate to get summary values (total sites)
        per.pop.meanFST <- aggregate(value~MinCoverage+pair, data=data_in[which(data_in$analysis==".fst"),], function(x) mean(x, na.rm = T))
        colnames(per.pop.meanFST)[3] <- "MeanFST"
        per.pop.sdFST <- aggregate(value~MinCoverage+pair, data=data_in[which(data_in$analysis==".fst"),], function(x) sd(x, na.rm = T))
        colnames(per.pop.sdFST)[3] <- "SdFST"
        per.pop.numSNPs <- aggregate(snpid~MinCoverage+pair, data=data_in[which(data_in$analysis==".fst"),], function(x) length(unique(x)))
        colnames(per.pop.numSNPs)[3] <- "NumSNPs"
        per.pop.numContigs <- aggregate(CHROM~MinCoverage+pair, data=data_in[which(data_in$analysis==".fst"),], function(x) length(unique(x)))
        colnames(per.pop.numContigs)[3] <- "NumContigs"
        
        #combine per-population summaries
        cov.perpair.table <- Reduce(function(x, y) merge(x, y, by=c("MinCoverage","pair")), list(per.pop.meanFST, per.pop.sdFST, per.pop.numSNPs, per.pop.numContigs))
        cov.perpair.table$MeanSNPsPerContig <- round((cov.perpair.table$NumSNPs / cov.perpair.table$NumContigs), 1)
        
        #across population summary 
        all.meanFST <- aggregate(MeanFST~MinCoverage, per.pop.meanFST, FUN=mean); colnames(all.meanFST)[2] <- "MeanFST"
        all.sdFST <- aggregate(MeanFST~MinCoverage, per.pop.meanFST, FUN=sd); colnames(all.sdFST)[2] <- "SdFST"
        all.numSNPs <- aggregate(snpid~MinCoverage, data=data_in[which(data_in$analysis==".fst"),], function(x) length(unique(x)))
        colnames(all.numSNPs)[2] <- "NumSNPs"
        all.numContigs <- aggregate(CHROM~MinCoverage, data=data_in[which(data_in$analysis==".fst"),], function(x) length(unique(x)))
        colnames(all.numContigs)[2] <- "NumContigs"
        
        #combine across-population summaries
        cov.allpairs.table <- Reduce(function(x, y) merge(x, y, by="MinCoverage", all=TRUE), list(all.meanFST, all.sdFST, all.numSNPs, all.numContigs))
        cov.allpairs.table$MeanSNPsPerContig <- round((cov.allpairs.table$NumSNPs / cov.allpairs.table$NumContigs), 1)
        
        #return per-population and across-population summaries
        return(list("cov.perpair.table"=cov.perpair.table, "cov.allpairs.table"=cov.allpairs.table))
    }
    
    #get summaries for all SNPs
    total.summary.out <- summarizePopoolation(postpop.master.long.bycov)
    cov.allpairs.table.total <- total.summary.out$cov.allpairs.table
    cov.perpair.table.total <- total.summary.out$cov.perpair.table
    cov.perpair.table.total$pair <- gsub(paste(project_name,"_",sep=""), "", cov.perpair.table.total$pair)
    
    #get summaries for SNPs called in all pools
    allpools <- data.frame("snpid"=as.character(called_allpops$snpid))
    postpop.master.long.bycov.allpools <- merge(postpop.master.long.bycov, allpools, all.x=FALSE, all.y=TRUE)
    allpools.summary.out <- summarizePopoolation(postpop.master.long.bycov.allpools)
    cov.allpairs.table.allpools <- allpools.summary.out$cov.allpairs.table
    cov.perpair.table.allpools <- allpools.summary.out$cov.perpair.table
    cov.perpair.table.allpools$pair <- gsub(paste(project_name,"_",sep=""), "", cov.perpair.table.allpools$pair)
        
    return(list("cov.allpairs.table.total"=cov.allpairs.table.total, 
                "cov.perpair.table.total"=cov.perpair.table.total, 
                "cov.allpairs.table.allpools"=cov.allpairs.table.allpools, 
                "cov.perpair.table.allpools"=cov.perpair.table.allpools,
                "postpop.master.long"=postpop.master.long.bycov,
                "postpop.master.long.allpools"=postpop.master.long.bycov.allpools))
    
}
