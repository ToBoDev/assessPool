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

postFST <- function(filetype, project_name, fst.calc, as, popcomb, strong_diff, p_cutoff, fasta_generation, ref_file, first_mincov, last_mincov, cov_step, lowP_removal, include_called_allpops){
  
  #file.copy(from=vcf_file, to=paste(working_dir,project_name,vcf_file, sep="/"), overwrite = TRUE)
  #file.copy(from=ref_file, to=paste(working_dir,project_name,ref_file, sep="/"), overwrite = TRUE) 
  

    
    setwd(paste(working_dir,"/",project_name,"/", fst.calc, sep=""))
  
    #add parameters to logfile
    write.log("Analysis Parameters:", paste(working_dir, project_name, "logs/analysis.log", sep="/"))
    write.log(paste("Strongly differentiated FST >=",strong_diff), paste(working_dir, project_name, "logs/analysis.log", sep="/"))
    write.log(paste("Significant FET <=", p_cutoff), paste(working_dir, project_name, "logs/analysis.log", sep="/"))
    
    C <- ncol(popcomb) #number of combinations
    combs <- list() #empty list for pairwise tags
    for(i in 1:C) combs[i] <- paste(project_name, popcomb[1,i], popcomb[2,i], sep='_') #pairwise tags
    df_names <- list(); df_index <- 1
    
    for(ft in filetype){
      for(i in 1:C){
        f=paste(combs[i], ft,sep='')
        if ("popoolation" %in% fst.calc){
          assign(f, read.table(f, sep='\t'))} #read in dataframe
        if ("poolfstat" %in% fst.calc){
          assign(f, read.table(f, sep=' '))}
        d=get(f) #get dataframe object
        colnames(d)[c(1,2,6)]=c("CHROM","POS","value"); d$pair <- combs[i]; d$analysis <- ft #reorganize columns
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
    
    if(fst.calc=="popoolation"){
      if(lowP_removal){
      postpop.master.long <- postpop.master.long %>% mutate(pair=unlist(pair)) %>% 
        pivot_wider(values_from = value, names_from = analysis, values_fn={mean}) %>% filter(.fet < p_cutoff) %>% 
        pivot_longer(cols = c(.fet,.fst), names_to = "analysis")
      }
    }
    
    #convert to wide format for export
    postpop.master.tmp <- postpop.master.long
    postpop.master.tmp$tag <- paste(postpop.master.tmp$pair, postpop.master.tmp$analysis,sep="")
    postpop.master.tmp$pair <- NULL; postpop.master.tmp$analysis <- NULL
    postpop.master.wide <- spread(postpop.master.tmp, tag, value)
    rm(postpop.master.tmp)
    
    #add columns from master dataframe
    as$snpid <- paste(as$CHROM, as$POS, sep="_")#; as$POS <- NULL
    #postpop.master.wide <- merge(as[,-which(names(as) %in% c("LEN","INS.len","DEL.len"))], postpop.master.wide[,-which(names(postpop.master.wide) %in% c("CHROM"))], by="snpid", all.y=TRUE, all.x=FALSE) 
    
    #convert above line to dplyr
    postpop.master.wide <- left_join(as[,-which(names(as) %in% c("LEN","INS.len","DEL.len"))], postpop.master.wide[,-which(names(postpop.master.wide) %in% c("CHROM"))])
    
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
    m <- paste0(nrow(postpop.master.wide), " informative variable sites retained after ", fst.calc)
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
    mkdirs(paste(working_dir,"/", project_name, "/output/", fst.calc, "/", sep=""))
    write.csv(popl.mast.export, paste0(working_dir,"/", project_name, "/output/", fst.calc, "/", project_name,"_allvar_post",fst.calc,".csv"), row.names=FALSE)
    m <- paste0("Exported all variable sites to ","output/", fst.calc, "/", project_name,"_allvar_post", fst.calc,".csv")
    message(m); write.log(m, paste(working_dir, project_name, "logs/analysis.log", sep="/"))
    
    #Test how many SNPs are called in all populations
    called_allpops <- postpop.master.wide[postpop.master.wide$fstNAs==0,] 
    m <- paste(nrow(called_allpops), "SNPs called in all pools.")
    if(nrow(called_allpops)==0){
      message("No SNPs shared accross all populations, we recommend using \"include_called_allpops=F\"")
    }
    message(m); write.log(m, paste(working_dir, project_name, "logs/analysis.log", sep="/"))
    
    
    
    if(include_called_allpops){

    #export sites called in all populations
    popl.allpops.exports <- data.frame(lapply(called_allpops, as.character), stringsAsFactors=FALSE)
    write.csv(popl.mast.export, paste0(working_dir,"/", project_name,"/output/", fst.calc, "/", project_name, "_calledAllPools_post",fst.calc,".csv"), row.names=FALSE)
    m <- paste0("Exported variable sites called in all pools to ","output/", fst.calc, "/", project_name,"_calledAllPools_post",fst.calc,".csv")
    message(m); write.log(m, paste(working_dir, project_name, "logs/analysis.log", sep="/"))
    
    }
    if (".fst" %in% filetype){
      
      #extract sites that show strong differentiation
      postpop.master.wide$DiffFST <- apply(postpop.master.wide[,grep(".fst", names(postpop.master.wide))], 1, function(x) any(x >= strong_diff))
      popl.appfx <- postpop.master.wide[which(postpop.master.wide$DiffFST==TRUE),]
      #popl.appfx <- na.omit(postpop.master.wide[apply(postpop.master.wide[,grep(".fst", names(postpop.master.wide))], 1, function(x) any(x >= strong_diff)),])
      m <- paste0(nrow(popl.appfx), " SNPs are strongly differentiated (FST>=", strong_diff, ") in at least one comparison.")
      message(m); write.log(m, paste(working_dir, project_name, "logs/analysis.log", sep="/"))
      
      popl.appfx.export <- data.frame(lapply(popl.appfx, as.character), stringsAsFactors=FALSE)
      write.csv(popl.appfx.export, paste0(working_dir,"/", project_name,"/output/", fst.calc, "/", project_name, "_strongly_differentiated_sites.csv"), row.names=FALSE)
      m <- paste0("Exported strongly differentiated sites to ","output/", fst.calc, "/", project_name,"_strongly_differentiated_sites.csv")
      message(m); write.log(m, paste(working_dir, project_name, "logs/analysis.log", sep="/"))
      
      #extract alternatively fixed sites 
      postpop.master.wide$FixedFST <- apply(postpop.master.wide[,grep(".fst", names(postpop.master.wide))], 1, function(x) any(x == 1.0))
      popl.fixed <- postpop.master.wide[which(postpop.master.wide$FixedFST==TRUE),]
      #popl.fixed <- na.omit(postpop.master.wide[apply(postpop.master.wide[,grep(".fst", names(postpop.master.wide))], 1, function(x) any(x == 1.0)),])
      m <- paste(nrow(popl.fixed), "sites are alternatively fixed (FST=1) in an least one comparison.")
      message(m); write.log(m, paste(working_dir, project_name, "logs/analysis.log", sep="/"))
      
      popl.fixed.export <- data.frame(lapply(popl.fixed, as.character), stringsAsFactors=FALSE)
      write.csv(popl.fixed.export, paste0(working_dir,"/", project_name,"/output/", fst.calc, "/", project_name, "_fixed_snps.csv"), row.names=FALSE)
      m <- paste0("Exported alternatively fixed sites to ","output/", fst.calc, "/", project_name,"_fixed_snps.csv")
      message(m); write.log(m, paste(working_dir, project_name, "logs/analysis.log", sep="/"))
      
      if (fasta_generation==TRUE){
        ref_file_df <- read.fasta(paste(working_dir, ref_file, sep="/"))
        ref_file_df <- data.frame(CHROM=names(ref_file_df), Seq=toupper(unlist(getSequence(ref_file_df, as.string=T))))
        
        #generate FASTA from strongly differentiated sites
        sd_fasta <- getFasta(popl.appfx, ref_file_df)
        write.fasta(sequences=as.list(sd_fasta$Seq), names=sd_fasta$CHROM, file.out=paste0(working_dir,"/", project_name,"/output/", fst.calc, "/", project_name, "_strongly_differentiated_sites_fasta.fa"), nbchar=1000)
        m <- paste0("Wrote ", nrow(sd_fasta), " contigs to /output/",fst.calc, "/", project_name, "_strongly_differentiated_sites_fasta.fa")
        message(m); write.log(m, paste(working_dir, project_name, "logs/analysis.log", sep="/"))
        
        #generate FASTA from alternatively fixed sites
        alt_fasta <- getFasta(popl.fixed, ref_file_df)
        write.fasta(sequences=as.list(alt_fasta$Seq), names=alt_fasta$CHROM, file.out=paste0(working_dir,"/", project_name,"/output/", fst.calc, "/", project_name, "_fixed_sites_fasta.fa"), nbchar=1000)
        m <- paste0("Wrote ", nrow(alt_fasta), " contigs to /output/", fst.calc, "/", project_name, "_fixed_sites_fasta.fa")
        message(m); write.log(m, paste(working_dir, project_name, "logs/analysis.log", sep="/"))
      }
    }
    
    if (".fet" %in% filetype){
      postpop.master.wide$LowP <- apply(postpop.master.wide[,grep(".fet", names(postpop.master.wide))], 1, function(x) any(x <= p_cutoff))
      popl.lowP <- postpop.master.wide[which(postpop.master.wide$LowP==TRUE),]
      #popl.fixed <- na.omit(postpop.master.wide[apply(postpop.master.wide[,grep(".fst", names(postpop.master.wide))], 1, function(x) any(x == 1.0)),])
      m <- paste(nrow(popl.lowP), " SNPs have a low p-value (p<=", p_cutoff, ") in at least one comparison.", sep="")
      message(m); write.log(m, paste(working_dir, project_name, "logs/analysis.log", sep="/"))
      
      popl.lowP.export <- data.frame(lapply(popl.lowP, as.character), stringsAsFactors=FALSE)
      write.csv(popl.lowP.export, paste(working_dir,"/", project_name,"/output/", fst.calc, "/", project_name, "_lowP_snps.csv", sep=""), row.names=FALSE)
      m <- paste("Exported low p-value sites to ","output/", fst.calc, "/", project_name,"_lowP_snps.csv", sep="")
      message(m); write.log(m, paste(working_dir, project_name, "logs/analysis.log", sep="/"))
    }
    
    
    ###############add min coverage level to dataframe
    
    covs <- rev(c(seq(first_mincov,last_mincov,cov_step))) 
    #convert to numeric and replace NAs with 0
    #postpop.master.wide[,grep("DP\\.", names(postpop.master.wide))] <- apply(postpop.master.wide[,grep("DP\\.", names(postpop.master.wide))], 1, as.character)
    #postpop.master.wide[,grep("DP\\.", names(postpop.master.wide))] <- apply(postpop.master.wide[,grep("DP\\.", names(postpop.master.wide))], 1, as.numeric)
    postpop.master.wide[,grep("DP\\.", names(postpop.master.wide))][is.na(postpop.master.wide[,grep("DP\\.", names(postpop.master.wide))])] <- 0
    df_names_total <- list(); df_index <- 1
    
    #message("")
    #assign min in-population coverage levels to SNPs
    
    tmp.idx <- c() #initialize list to keep track of rows already analyzed
    d <- data.frame()
    tmp.postpop <- postpop.master.wide[,-grep("RO\\.|AO\\.", names(postpop.master.wide))]
    tmp.postpop <- tmp.postpop[,-which(names(tmp.postpop) %in% c("CHROM","POS","NS","TDP","TYPE","AN","REF","ALT","Rlen","Alen"))]
    #remove uneeded columns
    #rm(postpop.master.wide) #no longer needed
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
    
    #append dataframes to build master dataframes (wide  and long) by coverage
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
    if(include_called_allpops){
      
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
    } else {
      
    allpools <- data.frame("snpid"=as.character(called_allpops$snpid))
    postpop.master.long.bycov.allpools <- merge(postpop.master.long.bycov, allpools, all.x=FALSE, all.y=TRUE)
    allpools.summary.out <- postpop.master.long.bycov.allpools
    cov.allpairs.table.allpools <- allpools.summary.out$cov.allpairs.table
    cov.perpair.table.allpools <- allpools.summary.out$cov.perpair.table
    
    return(list("cov.allpairs.table.total"=cov.allpairs.table.total, 
                "cov.perpair.table.total"=cov.perpair.table.total, 
                "cov.allpairs.table.allpools"=cov.allpairs.table.allpools, 
                "cov.perpair.table.allpools"=cov.perpair.table.allpools,
                "postpop.master.long"=postpop.master.long.bycov,
                "postpop.master.long.allpools"=postpop.master.long.bycov.allpools))
    }
   
}

#wrapper functions to run postFST()
get_popoolation_out <- function(working_dir, project_name, fst.calc){
  setwd(paste(working_dir, project_name, fst.calc, sep="/"))
  
  if (length(list.files(pattern="*.fst"))==0 & length(list.files(pattern="*.fet"))==0){
    message("\n\nERROR: No .fst/.fet PoPoolation2 output files found. Please run PoPoolation2.")
  } else {
    if (length(list.files(pattern="*.fst"))>0 & length(list.files(pattern="*.fet"))>0){ filetype = c(".fst",".fet")
    }else if(length(list.files(pattern="*.fst"))>0 & length(list.files(pattern="*.fet"))==0){ filetype = c(".fst")
    }else if (length(list.files(pattern="*.fst"))==0 & length(list.files(pattern="*.fet"))>0){ filetype = c(".fet")}
    
    pa_list_out_pp <- postFST(fst.calc="popoolation",
                              filetype=filetype,
                              project_name=project_name,
                              as=as, 
                              popcomb=popcomb, 
                              strong_diff=strong_diff,
                              p_cutoff=p_cutoff,
                              fasta_generation=fasta_generation,
                              ref_file=ref_file,
                              first_mincov=first_mincov,
                              last_mincov=last_mincov,
                              cov_step=cov_step,
                              lowP_removal=lowP_removal,
                              include_called_allpops=include_called_allpops)
    
    #extract needed data from return values
    assign(paste0("cov.perpair.table.total_", fst.calc), pa_list_out_pp$cov.perpair.table.total)
    assign(paste0("postpop.master.long_", fst.calc), pa_list_out_pp$postpop.master.long)
    assign(paste0("cov.allpairs.table.total_", fst.calc), pa_list_out_pp$cov.allpairs.table.total)
    
    if(include_called_allpops){
    assign(paste0("cov.allpairs.table.allpools_", fst.calc), pa_list_out_pp$cov.allpairs.table.allpools)
    assign(paste0("cov.perpair.table.allpools_", fst.calc), pa_list_out_pp$cov.perpair.table.allpools) 
    assign(paste0("postpop.master.long.allpools_", fst.calc), pa_list_out_pp$postpop.master.long.allpools) 
    }
    
    #save all output tables
    write.csv(get(paste0("cov.allpairs.table.total_", fst.calc)), 
              file = paste(working_dir, project_name,"output", fst.calc, paste0("cov.allpairs.table.total_", fst.calc,".csv"), sep="/"),
              row.names=F)
    
    write.csv(get(paste0("cov.perpair.table.total_", fst.calc)),
              file = paste(working_dir, project_name,"output", fst.calc, paste0("cov.perpair.table.total_", fst.calc,".csv"), sep="/"),
              row.names=F)
    
    write.csv(get(paste0("postpop.master.long_", fst.calc)),
              file = paste(working_dir, project_name,"output", fst.calc, paste0("postpop.master.long_", fst.calc,".csv"), sep="/"),
              row.names=F)
    
    if(include_called_allpops){
     write.csv(get(paste0("cov.allpairs.table.allpools_", fst.calc)), 
              file = paste(working_dir, project_name,"output", fst.calc, paste0("cov.allpairs.table.allpools_", fst.calc,".csv"), sep="/"),
              row.names=F)
    
    write.csv(get(paste0("cov.perpair.table.allpools_", fst.calc)), 
              file = paste(working_dir, project_name,"output", fst.calc, paste0("cov.perpair.table.allpools_", fst.calc,".csv"), sep="/"),
              row.names=F)
    write.csv(get(paste0("postpop.master.long.allpools_", fst.calc)), 
              file = paste(working_dir, project_name,"output", fst.calc, paste0("postpop.master.long.allpools_", fst.calc,".csv"), sep="/"),
              row.names=F)
    }
  }
}

get_poolfstat_out <- function(working_dir, project_name, fst.calc){
  setwd(paste(working_dir, project_name, fst.calc, sep="/"))
  
  if (length(list.files(pattern=".fst"))==0 & length(list.files(pattern=".fet"))==0){
    message("\n\nERROR: No poolfstat output files found. Please run poolfstat.")
  } else {
    if (length(list.files(pattern=".fst"))>0 & length(list.files(pattern=".fet"))>0){ filetype = c(".fst",".fet")
    } else if(length(list.files(pattern=".fst"))>0 & length(list.files(pattern=".fet"))==0){ filetype = c(".fst")
    } else if (length(list.files(pattern=".fst"))==0 & length(list.files(pattern=".fet"))>0){ filetype = c(".fet")}
    
    pa_list_out_pf <- postFST(fst.calc="poolfstat",
                              filetype=filetype,
                              project_name=project_name,
                              as=as, 
                              popcomb=popcomb, 
                              strong_diff=strong_diff,
                              p_cutoff=p_cutoff,
                              fasta_generation=fasta_generation,
                              ref_file=ref_file,
                              first_mincov=first_mincov,
                              last_mincov=last_mincov,
                              cov_step=cov_step,
                              lowP_removal=lowP_removal,
                              include_called_allpops=include_called_allpops)
    
    #extract needed data from return values
    assign(paste0("cov.perpair.table.total_", fst.calc), pa_list_out_pf$cov.perpair.table.total)
    assign(paste0("postpop.master.long_", fst.calc), pa_list_out_pf$postpop.master.long)
    assign(paste0("cov.allpairs.table.total_", fst.calc), pa_list_out_pf$cov.allpairs.table.total)
    
    if(include_called_allpops){
      assign(paste0("cov.allpairs.table.allpools_", fst.calc), pa_list_out_pf$cov.allpairs.table.allpools)
      assign(paste0("cov.perpair.table.allpools_", fst.calc), pa_list_out_pf$cov.perpair.table.allpools) 
      assign(paste0("postpop.master.long.allpools_", fst.calc), pa_list_out_pf$postpop.master.long.allpools) 
    }
    
    #save all output tables
    write.csv(get(paste0("cov.allpairs.table.total_", fst.calc)), 
              file = paste(working_dir, project_name,"output", fst.calc, paste0("cov.allpairs.table.total_", fst.calc,".csv"), sep="/"),
              row.names=F)
    
    write.csv(get(paste0("cov.perpair.table.total_", fst.calc)),
              file = paste(working_dir, project_name,"output", fst.calc, paste0("cov.perpair.table.total_", fst.calc,".csv"), sep="/"),
              row.names=F)
    
    write.csv(get(paste0("postpop.master.long_", fst.calc)),
              file = paste(working_dir, project_name,"output", fst.calc, paste0("postpop.master.long_", fst.calc,".csv"), sep="/"),
              row.names=F)
    
    if(include_called_allpops){
      write.csv(get(paste0("cov.allpairs.table.allpools_", fst.calc)), 
                file = paste(working_dir, project_name,"output", fst.calc, paste0("cov.allpairs.table.allpools_", fst.calc,".csv"), sep="/"),
                row.names=F)
      
      write.csv(get(paste0("cov.perpair.table.allpools_", fst.calc)), 
                file = paste(working_dir, project_name,"output", fst.calc, paste0("cov.perpair.table.allpools_", fst.calc,".csv"), sep="/"),
                row.names=F)
      write.csv(get(paste0("postpop.master.long.allpools_", fst.calc)), 
                file = paste(working_dir, project_name,"output", fst.calc, paste0("postpop.master.long.allpools_", fst.calc,".csv"), sep="/"),
                row.names=F)
    }
  }
}  
wrangle_out <- function(working_dir, project_name, fst.calc){
  
  setwd(paste(working_dir, project_name, fst.calc, sep="/"))
  
  if("popoolation" %in% fst.calc){
    
    fst.c <- "popoolation"
    cov.allpairs.table.total_popoolation <<- read.csv(file = paste(working_dir, project_name,"output", fst.c,
                                                                   paste0("cov.allpairs.table.total_", fst.c,".csv"), sep="/"))
    
    cov.perpair.table.total_popoolation <<- read.csv(file = paste(working_dir, project_name,"output", fst.c,
                                                                  paste0("cov.perpair.table.total_", fst.c,".csv"), sep="/")) %>% mutate(pair=as.character(pair))
    
    postpop.master.long_popoolation <<- read.csv(file = paste(working_dir, project_name,"output", fst.c,
                                                              paste0("postpop.master.long_", fst.c,".csv"), sep="/")) %>% mutate(pair=as.character(pair))
    
    if(include_called_allpops){
    postpop.master.long.allpools_popoolation <<- read.csv(file = paste(working_dir, project_name,"output", fst.c,
                                                                       paste0("postpop.master.long.allpools_", fst.c,".csv"), sep="/")) %>% mutate(pair=as.character(pair))
    
    cov.allpairs.table.allpools_popoolation <<- read.csv(file = paste(working_dir, project_name,"output", fst.c,
                                                                      paste0("cov.allpairs.table.allpools_", fst.c,".csv"), sep="/")) 
    
    cov.perpair.table.allpools_popoolation <<- read.csv(file = paste(working_dir, project_name,"output", fst.c,
                                                                     paste0("cov.perpair.table.allpools_", fst.c,".csv"), sep="/")) %>% mutate(pair=as.character(pair))
    } else {
      cov.allpairs.table.allpools_popoolation <<- NA
      cov.perpair.table.allpools_popoolation <<- NA
      postpop.master.long.allpools_popoolation <<- NA
    }
    
    postpop.master.wide_pp <<- spread(postpop.master.long_popoolation, analysis, value)
    postpop.master.wide_pp$pair <<- gsub(paste0(project_name,   "_"), "", postpop.master.wide_pp$pair) 
    postpop.master.wide_pp <<- separate(postpop.master.wide_pp, pair, c("popA","popB"), sep="_", remove=F)
    
    cov.perpair.table.total_pp <<- cov.perpair.table.total_popoolation[order(cov.perpair.table.total_popoolation$MinCoverage),]
    cov.perpair.table.total_pp <<- separate(cov.perpair.table.total_pp, pair, c("popA","popB"), sep="_", remove=F)
    cov.perpair.table.total_pp$NumSNPs <<- as.numeric(cov.perpair.table.total_pp$NumSNPs)
    cov.perpair.table.total_pp$NumContigs <<- as.numeric(cov.perpair.table.total_pp$NumContigs)
    cov.perpair.table.total_pp$pair <<- as.character(cov.perpair.table.total_pp$pair)
    
    }
  
  
  
  if("poolfstat" %in% fst.calc){
    
    fst.c <- "poolfstat"
    cov.allpairs.table.total_poolfstat <<- read.csv(file = paste(working_dir, project_name,"output", fst.c,
                                                                 paste0("cov.allpairs.table.total_", fst.c,".csv"), sep="/"))
    
    cov.perpair.table.total_poolfstat <<- read.csv(file = paste(working_dir, project_name,"output", fst.c,
                                                                paste0("cov.perpair.table.total_", fst.c,".csv"), sep="/")) %>% mutate(pair=as.character(pair))
    
    postpop.master.long_poolfstat <<- read.csv(file = paste(working_dir, project_name,"output", fst.c,
                                                            paste0("postpop.master.long_", fst.c,".csv"), sep="/")) %>% mutate(pair=as.character(pair))
    
    if(include_called_allpops){
    cov.allpairs.table.allpools_poolfstat <<- read.csv(file = paste(working_dir, project_name,"output", fst.c,
                                                                    paste0("cov.allpairs.table.allpools_", fst.c,".csv"), sep="/"))
    
    cov.perpair.table.allpools_poolfstat <<- read.csv(file = paste(working_dir, project_name,"output", fst.c,
                                                                   paste0("cov.perpair.table.allpools_", fst.c,".csv"), sep="/")) %>% mutate(pair=as.character(pair))
    
    postpop.master.long.allpools_poolfstat <<- read.csv(file = paste(working_dir, project_name,"output", fst.c,
                                                                     paste0("postpop.master.long.allpools_", fst.c,".csv"), sep="/")) %>% mutate(pair=as.character(pair))
    } else {
      cov.allpairs.table.allpools_poolfstat <<- NA
      cov.perpair.table.allpools_poolfstat <<- NA
      postpop.master.long.allpools_poolfstat <<- NA
    }
    
    postpop.master.wide_pf <<- spread(postpop.master.long_poolfstat, analysis, value)
    postpop.master.wide_pf$pair <<- gsub(paste0(project_name,"_"), "", postpop.master.wide_pf$pair) 
    postpop.master.wide_pf <<- separate(postpop.master.wide_pf, pair, c("popA","popB"), sep="_", remove=F)
    
    cov.perpair.table.total_pf <<- cov.perpair.table.total_poolfstat[order(cov.perpair.table.total_poolfstat$MinCoverage),]
    cov.perpair.table.total_pf <<- separate(cov.perpair.table.total_pf, pair, c("popA","popB"), sep="_", remove=F)
    cov.perpair.table.total_pf$NumSNPs <<- as.numeric(cov.perpair.table.total_pf$NumSNPs)
    cov.perpair.table.total_pf$NumContigs <<- as.numeric(cov.perpair.table.total_pf$NumContigs)
    cov.perpair.table.total_pf$pair <<- as.character(cov.perpair.table.total_pf$pair)

        }
}
