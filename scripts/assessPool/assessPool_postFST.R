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

summarizePopoolation <- function(data_in){
  
  #aggregate to get summary values (total sites)
  
 # data_in %>% filter(!is.na(.fst), MinCoverage >= first_mincov) %>% mutate(cov_bin = case_when(MinCoverage > covs[1] ~ "A"))
  
  cov.perpair.table <- data_in %>% 
    filter(!is.na(.fst), MinCoverage >= first_mincov) %>% dplyr::group_by(MinCoverage,pair) %>% 
    dplyr::summarize("MeanFST"=mean(.fst), 
              "SdFST"=sd(.fst),
              "SeFST"=sd(.fst)/sqrt(n()),
              "MeanDP"=mean(fst.dp), 
              "NumSNPs"=n_distinct(snpid), 
              "NumContigs"=n_distinct(contig), 
              "MeanSNPsPerContig"=(NumSNPs/NumContigs)) %>% 
    dplyr::mutate(ci_lower= MeanFST - qt(0.975, n() -1) * SeFST,
           ci_upper= MeanFST + qt(0.975, n() -1) * SeFST)
  
  cov.allpairs.table <- data_in %>% 
    filter(!is.na(.fst), MinCoverage >= first_mincov)  %>% group_by(MinCoverage) %>% 
    dplyr::summarize("MeanFST"=mean(.fst), 
              "SdFST"=sd(.fst),
              "SeFST"=sd(.fst)/sqrt(n()),
              "MeanDP"=mean(fst.dp), 
              "NumSNPs"=n_distinct(snpid), 
              "NumContigs"=n_distinct(contig), 
              "MeanSNPsPerContig"=(NumSNPs/NumContigs)) %>% 
    dplyr::mutate(ci_lower= MeanFST - qt(0.975, n() -1) * SeFST,
           ci_upper= MeanFST + qt(0.975, n() -1) * SeFST)
  
  cov.perpair.table.allpools <- data_in %>% 
    filter(!is.na(.fst), MinCoverage >= first_mincov) %>% filter(snpid %in% called_allpops$snpid) %>% group_by(MinCoverage,pair) %>% 
    dplyr::summarize("MeanFST"=mean(.fst), 
              "SdFST"=sd(.fst),
              "SeFST"=sd(.fst)/sqrt(n()),
              "MeanDP"=mean(fst.dp), 
              "NumSNPs"=n_distinct(snpid), 
              "NumContigs"=n_distinct(contig), 
              "MeanSNPsPerContig"=(NumSNPs/NumContigs)) %>% 
    dplyr::mutate(ci_lower= MeanFST - qt(0.975, n() -1) * SeFST,
           ci_upper= MeanFST + qt(0.975, n() -1) * SeFST)
  
  cov.allpairs.table.allpools <- data_in %>% 
    filter(!is.na(.fst), MinCoverage >= first_mincov) %>% filter(snpid %in% called_allpops$snpid) %>% group_by(MinCoverage) %>% 
    dplyr::summarize("MeanFST"=mean(.fst), 
              "SdFST"=sd(.fst),
              "SeFST"=sd(.fst)/sqrt(n()),
              "MeanDP"=mean(fst.dp), 
              "NumSNPs"=n_distinct(snpid), 
              "NumContigs"=n_distinct(contig), 
              "MeanSNPsPerContig"=(NumSNPs/NumContigs)) %>% 
    dplyr::mutate(ci_lower= MeanFST - qt(0.975, n() -1) * SeFST,
           ci_upper= MeanFST + qt(0.975, n() -1) * SeFST)
  

  #return per-population and across-population summaries
  return(list("cov.perpair.table"=cov.perpair.table, 
              "cov.allpairs.table"=cov.allpairs.table,
              "cov.perpair.table.allpools"=cov.perpair.table.allpools, 
              "cov.allpairs.table.allpools"=cov.allpairs.table.allpools))
}


postFST <- function(filetype, 
                    project_name, 
                    fst.calc, 
                    as, 
                    popcomb, 
                    strong_diff, 
                    p_cutoff, 
                    fasta_generation, 
                    ref_file, 
                    first_mincov, 
                    last_mincov, 
                    cov_step, 
                    lowP_removal, 
                    lowP_to_zeros, 
                    lowDP_to_zeros, 
                    include_called_allpops, 
                    min.pairwise.prop,
                    assessPool_thinning){

    
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
        colnames(d)[c(1,2,5,6)]=c("contig","pos","fst.dp","value"); d$pair <- combs[i]; d$analysis <- ft #reorganize columns
        d$value <- as.numeric(gsub(".*=","",as.character(d$value))) #remove unwanted characters, convert to numeric
        d$snpid <- paste(d[,1], d[,2], sep = '_') #add snpid column
        d <- d[,c("contig","pos","fst.dp","value","pair","analysis","snpid")] #drop unneeded columns
        if (ft==".fet"){ d$value <- 10^(-(d$value)) } #convert log(p) to p-value
        assign(f,d)
        df_names[df_index] <- paste(combs[i],ft,sep=""); df_index <- df_index+1
      }
    }
    
    

    #save dataframe in long format for future data manipulation
    postpop.master.long.backup <- bind_rows(mget(unlist(df_names))) %>% mutate(pair=unlist(pair))
    for(ft in filetype){
    rm(list=ls(pattern=paste("*.",ft,sep='')))}
    
    rm(d) 
    
    # read in fasta file as character vector
    fasta_lines <- readLines(paste0(working_dir, "/",ref_file))
    
    # initialize empty lists to store sequence names and lengths
    seq_names <- list()
    seq_lengths <- list()
    
    fa_headers <- which(grepl("^>", fasta_lines))
    # loop through each line in fasta file
    for (i in 1:length(fa_headers)) {
      
      line <- fasta_lines[fa_headers[i]]
      
      # check if line starts with ">"
      if (substr(line, 1, 1) == ">") {
        # get sequence name
        seq_name <- gsub("^>", "", line)
        
        # get length of sequence
        if(i < length(fa_headers)){
        seq_length <- nchar(paste(fasta_lines[(fa_headers[i]+1):(fa_headers[i+1]-1)], collapse=""))
        } else { 
          seq_length <- nchar(paste(fasta_lines[(fa_headers[i]+1):length(fasta_lines)], collapse=""))
        }
        
        # append sequence name and length to lists
        seq_names[[length(seq_names)+1]] <- seq_name
        seq_lengths[[length(seq_lengths)+1]] <- seq_length
      }
    }
    
    # combine sequence names and lengths into data frame
    ref_names_lengths <- data.frame("contig"=as.factor(unlist(seq_names)), "contig_length"=as.numeric(unlist(seq_lengths)), stringsAsFactors = FALSE)
    
    postpop.master.long <- left_join(postpop.master.long.backup, ref_names_lengths, by="contig") %>% 
      select(snpid, contig, pos, contig_length, fst.dp, pair, analysis, value) %>% mutate(contig=as.character(contig))
    
    # postpop.master.long.backup <- postpop.master.long.backup %>% 
    #   mutate(contig_length = as.numeric(str_extract(CHROM_full,  "(?<=length_)\\d+(?=\\_cov)"))) %>% 
    #   mutate(contig = str_extract(CHROM_full,  "NODE_\\d+")) %>%
    #   mutate(contig_cov = as.numeric(str_extract(CHROM_full,  "\\d+\\.\\d+$"))) %>% 
    #   select(snpid, contig, POS, contig_length, contig_cov, fst.dp, pair, analysis, value, -CHROM_full)
    
    #postpop.master.long <- postpop.master.long.backup
    
    #pivot to columns for .fst and .fet (if PoPoolation2)
    postpop.master.long <- postpop.master.long %>% 
      mutate(pair=gsub(paste0(project_name,"_"),"",pair)) %>% 
      pivot_wider(names_from=analysis, values_from=value)
    
    #convert nonsignificant Fst calls to 0
    if(fst.calc=="popoolation"){
      if(lowP_to_zeros){
    postpop.master.long <- postpop.master.long %>% 
      mutate(.fst = if_else(.fet > p_cutoff, 0, .fst))
      }
    }
    
    #convert low depth pairwise Fst calls to 0
    if(lowDP_to_zeros){
      postpop.master.long <- postpop.master.long %>% 
        mutate(.fst = if_else(fst.dp < first_mincov, 0, .fst))
    }
    
    #add MinCoverage column
    covs <<- c(seq(first_mincov,last_mincov,cov_step))
    postpop.master.long <- postpop.master.long %>% 
        mutate(CovStepCutoff = cut(fst.dp, breaks = c(0,(covs-1), max(fst.dp)))) %>% 
      mutate(CovStepCutoff = as.numeric(str_extract(CovStepCutoff, "(?<=\\()[^,]+"))+1) %>% 
      mutate(CovStepCutoff = case_when(CovStepCutoff==1 ~ 0, CovStepCutoff==CovStepCutoff ~ CovStepCutoff))
    

    
    #if removing all nonsignificant Fst calls
    if(fst.calc=="popoolation"){
        
      fail.fet <- unique((postpop.master.long %>% filter(.fet > p_cutoff) %>% filter(fst.dp > first_mincov))$snpid)
      total.loci <- unique(postpop.master.long$snpid)
      m <- message(paste0(length(fail.fet), " of ", length(total.loci), " had low significance (p>", p_cutoff,") in PoPoolation's Fisher's Exact Test in at least one comparison"))
      message(m); write.log(m, paste(working_dir, project_name, "logs/analysis.log", sep="/"))
      
      if(lowP_removal){
        postpop.master.long <- postpop.master.long %>% filter(.fet < p_cutoff) %>% filter(fst.dp > first_mincov)
      }
    }
    
    #add a column to get relative position of SNP along contig
    postpop.master.long <- postpop.master.long %>% 
      mutate(dist_contig_center=(1-(abs((contig_length/2)-pos)/contig_length)*2))
    
    #convert to wide format for export
    if(fst.calc=="popoolation"){
      
    postpop.master.wide.fst <- postpop.master.long %>%
      mutate(tag=paste0(pair,".fst")) %>% 
      select(-fst.dp, -.fet, -CovStepCutoff, -pair) %>%
      pivot_wider(values_from=.fst, names_from=tag) %>% 
      mutate(fstNAs = rowSums(is.na(.))) %>% 
      filter(fstNAs < C*min.pairwise.prop)
    
    #generate summary statistics to rank SNP postion along contig in order to select highest ranking SNP (thins data to one SNP per contig)
    wide_to_add <- postpop.master.long %>% 
      group_by(snpid) %>% 
      dplyr::summarise(mean_fst.dp=mean(fst.dp),
                mean_MinCov=mean(CovStepCutoff),
                mean_fst=mean(.fst))
    
    #Add summary statistics
    postpop.master.wide.fst.sum <- left_join(postpop.master.wide.fst, wide_to_add, by = c("snpid")) %>% 
      select(snpid, contig, pos, contig_length, dist_contig_center, mean_fst, mean_fst.dp, mean_MinCov, fstNAs, everything())
      
    #perform ranking
    # ranking order: dist_contig_center mean_fst.dp fstNAs mean_fst
    weights <- c(-1, -1, 1, -1)
    df_score <- postpop.master.wide.fst.sum %>%
      mutate(score = weights[1]*dist_contig_center*0 + weights[2]*(mean_fst.dp/first_mincov) + weights[3]*(1-(fstNAs/C)) + weights[4]*(100*mean_fst))
    
    # rank positions based on their score within each contig
    df_ranked <- df_score %>%
      group_by(contig) %>%
      mutate(rank = dense_rank(score)) %>%
      arrange(contig, rank)
    
    postpop.master.wide.fst <- df_ranked
    
    if(assessPool_thinning){
    # select top position for each contig
    postpop.master.wide.fst <- df_ranked %>%
      filter(rank == 1)
    }
    
    postpop.master.wide.fet <- postpop.master.long %>% 
      mutate(tag=paste0(pair,".fet")) %>% 
      select(-pair, -.fst, -fst.dp, -CovStepCutoff) %>% 
      pivot_wider(values_from=.fet, names_from=tag) %>% 
      filter(snpid %in% postpop.master.wide.fst$snpid)
    }
    
    if(fst.calc=="poolfstat"){
      
      postpop.master.wide.fst <- postpop.master.long %>%
        mutate(tag=paste0(pair,".fst")) %>% 
        select(-fst.dp , -CovStepCutoff, -pair) %>%
        pivot_wider(values_from=.fst, names_from=tag) %>% 
        mutate(fstNAs = rowSums(is.na(.))) %>% 
        filter(fstNAs < C*min.pairwise.prop)
      
      wide_to_add <- postpop.master.long %>% 
        group_by(snpid) %>% 
        dplyr::summarise(mean_fst.dp=mean(fst.dp),
                  mean_MinCov=mean(CovStepCutoff),
                  mean_fst=mean(.fst))
      
      #Add summary statistics
      postpop.master.wide.fst.sum <- left_join(postpop.master.wide.fst, wide_to_add, by = c("snpid")) %>% 
        select(snpid, contig, pos, contig_length, dist_contig_center, mean_fst, mean_fst.dp, mean_MinCov, fstNAs, everything())
      
      #perform ranking
      # ranking order: dist_contig_center mean_fst.dp fstNAs mean_fst
      weights <- c(-1, -1, 1, -1)
      df_score <- postpop.master.wide.fst.sum %>%
        mutate(score = weights[1]*dist_contig_center + 
                 weights[2]*(mean_fst.dp/first_mincov) + 
                 weights[3]*(1-(fstNAs/C)) + 
                 weights[4]*(100*mean_fst))
      
      # rank positions based on their score within each contig
      df_ranked <- df_score %>%
        group_by(contig) %>%
        mutate(rank = dense_rank(score)) %>%
        arrange(contig, rank)
      
      postpop.master.wide.fst <- df_ranked
      
      if(assessPool_thinning){
      # select top position for each contig
      postpop.master.wide.fst <- df_ranked %>%
        filter(rank == 1)
      }
      # postpop.master.wide.fst <- postpop.master.long %>% 
      #   mutate(tag=paste0(pair,".fst")) %>% select(-pair, -fst.dp, -MinCoverage) %>% pivot_wider(values_from=.fst, names_from=tag)
    }
    
        #add columns from master dataframe
    as$snpid <- paste(as$CHROM, as$POS, sep="_")#; as$POS <- NULL
    #postpop.master.wide <- merge(as[,-which(names(as) %in% c("LEN","INS.len","DEL.len"))], postpop.master.wide[,-which(names(postpop.master.wide) %in% c("CHROM"))], by="snpid", all.y=TRUE, all.x=FALSE) 
    
    ##MAYBE 
    postpop.master.wide <- right_join((as %>% select(-"LEN",-"INS.len",-"DEL.len")), 
                                     (postpop.master.wide.fst %>% select(-rank,-score)),
                                     by="snpid")
    
    # #filter out monomorphic SNPs with all pairwise Fst=0
    # if (".fst" %in% filetype){
    #   postpop.master.wide$fst0s <- apply(postpop.master.wide[,grep(".fst", names(postpop.master.wide))], 1, function(x) length(which(x==0)))  #counts instances of Fst=0
    #   postpop.master.wide <- postpop.master.wide[which(postpop.master.wide$fst0s < C),]; postpop.master.wide$fst0s <- NULL
    # }
    
    # #keep track of which populations have a SNP at this position (slow...)
    # postpop.master.wide$popIncl <- apply(postpop.master.wide[,grep(filetype[1],names(postpop.master.wide))], 1, function(x) names(postpop.master.wide[,grep(filetype[1],names(postpop.master.wide))])[which(!is.na(x))])
    # 
    # postpop.master.wide$popIncl <- apply(postpop.master.wide[,"popIncl",drop=F], 1, getPopNames)
    
    #count number of populations in which a SNP is not called; should not be any NAs but to double check
    # postpop.master.wide$fstNAs <- apply(postpop.master.wide[,grep(filetype[1],names(postpop.master.wide))], 1, function(x) sum(is.na(x)))
    # postpop.master.wide <- postpop.master.wide[!(postpop.master.wide$fstNAs >= C),] # remove rows without SNPS called in any pair (i.e. those that did not pass through filters)
    m <- paste0(nrow(postpop.master.wide), " informative variable sites retained after ", fst.calc)
    message(m); write.log(m, paste(working_dir, project_name, "logs/analysis.log", sep="/"))
    
    #apply same filters to data stored in long format
    postpop.master.long <- postpop.master.long %>% filter(snpid %in% postpop.master.wide$snpid)
    # merge(postpop.master.long, snpids_filter, all=FALSE)
    
    #clean up environment and organize data
    postpop.master.wide$TYPE <- as.factor(postpop.master.wide$TYPE); postpop.master.long$snpid <- as.factor(postpop.master.long$snpid)
    postpop.master.wide <- postpop.master.wide %>% arrange(snpid) #sort by CHROM/POS
    
    ####################### File export
    
    #export all variable sites
    #popl.mast.export <- data.frame(lapply(postpop.master.wide, as.character), stringsAsFactors=FALSE)
    mkdirs(paste(working_dir,"/", project_name, "/output/", fst.calc, "/", sep=""))
    write.csv(postpop.master.wide, paste0(working_dir,"/", project_name, "/output/", fst.calc, "/", project_name,"_allvar_post",fst.calc,".csv"), row.names=FALSE)
    m <- paste0("Exported all variable sites to ","output/", fst.calc, "/", project_name,"_allvar_post", fst.calc,".csv")
    message(m); write.log(m, paste(working_dir, project_name, "logs/analysis.log", sep="/"))
    
    #Test how many SNPs are called in all populations
    called_allpops <<- postpop.master.wide[postpop.master.wide$fstNAs==0,] 
    m <- paste(nrow(called_allpops), "SNPs called in all pools.")
    if(nrow(called_allpops)==0){
      message("No SNPs shared accross all populations, we recommend using \"include_called_allpops=F\"")
    }
    message(m); write.log(m, paste(working_dir, project_name, "logs/analysis.log", sep="/"))
    
    if(include_called_allpops){

    #export sites called in all populations
    #popl.allpops.exports <- data.frame(lapply(called_allpops, as.character), stringsAsFactors=FALSE)
    write.csv(called_allpops, paste0(working_dir,"/", project_name,"/output/", fst.calc, "/", project_name, "_calledAllPools_post",fst.calc,".csv"), row.names=FALSE)
    m <- paste0("Exported variable sites called in all pools to ","output/", fst.calc, "/", project_name,"_calledAllPools_post",fst.calc,".csv")
    message(m); write.log(m, paste(working_dir, project_name, "logs/analysis.log", sep="/"))
    
    }
    if (".fst" %in% filetype){
      
      #extract sites that show strong differentiation
      popl.appfx <-  postpop.master.wide %>% filter_at(vars(contains(".fst")), any_vars(. > strong_diff))
      
      #popl.appfx <- na.omit(postpop.master.wide[apply(postpop.master.wide[,grep(".fst", names(postpop.master.wide))], 1, function(x) any(x >= strong_diff)),])
      m <- paste0(nrow(popl.appfx), " SNPs are strongly differentiated (FST>=", strong_diff, ") in at least one comparison.")
      message(m); write.log(m, paste(working_dir, project_name, "logs/analysis.log", sep="/"))
      
      write.csv(popl.appfx, paste0(working_dir,"/", project_name,"/output/", fst.calc, "/", project_name, "_strongly_differentiated_sites.csv"), row.names=FALSE)
      m <- paste0("Exported strongly differentiated sites to ","output/", fst.calc, "/", project_name,"_strongly_differentiated_sites.csv")
      message(m); write.log(m, paste(working_dir, project_name, "logs/analysis.log", sep="/"))
      
      #extract alternatively fixed sites 
      popl.fixed <- postpop.master.wide %>% filter_at(vars(contains(".fst")), any_vars(. >= 0.99))
      #popl.fixed <- na.omit(postpop.master.wide[apply(postpop.master.wide[,grep(".fst", names(postpop.master.wide))], 1, function(x) any(x == 1.0)),])
      m <- paste(nrow(popl.fixed), "sites are alternatively fixed (FST=1) in an least one comparison.")
      message(m); write.log(m, paste(working_dir, project_name, "logs/analysis.log", sep="/"))
      
      write.csv(popl.fixed, paste0(working_dir,"/", project_name,"/output/", fst.calc, "/", project_name, "_fixed_snps.csv"), row.names=FALSE)
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
    
    #add additional rows for each MinCov passed by a SNP
    cov_tmp_df <- list()
    
    for(i in 1:length(covs)){
      cov_tmp_df[[i]] <- postpop.master.long %>% filter(CovStepCutoff > (covs[i]-1)) %>% mutate(MinCoverage=covs[i])
    }
    
    postpop.master.long <- bind_rows(cov_tmp_df)
    
    #get summaries for all SNPs
    total.summary.out <- summarizePopoolation(postpop.master.long)
    
    cov.allpairs.table.total <- total.summary.out$cov.allpairs.table
    cov.perpair.table.total <- total.summary.out$cov.perpair.table
    #cov.perpair.table.total$pair <- gsub(paste(project_name,"_",sep=""), "", cov.perpair.table.total$pair)
    
      if(include_called_allpops){
      #get summaries for SNPs called in all pools 
      postpop.master.long.allpools <- postpop.master.long %>% filter(snpid %in% called_allpops$snpid)
    
      cov.allpairs.table.allpools <- total.summary.out$cov.allpairs.table.allpools
      cov.perpair.table.allpools <- total.summary.out$cov.perpair.table.allpools
      #cov.perpair.table.allpools$pair <- gsub(paste(project_name,"_",sep=""), "", cov.perpair.table.allpools$pair)
      
      return(list("cov.allpairs.table.total"=cov.allpairs.table.total, 
                  "cov.perpair.table.total"=cov.perpair.table.total, 
                  "cov.allpairs.table.allpools"=cov.allpairs.table.allpools, 
                  "cov.perpair.table.allpools"=cov.perpair.table.allpools,
                  "postpop.master.long"=postpop.master.long,
                  "postpop.master.long.allpools"=postpop.master.long.allpools))
      } else {
      
      return(list("cov.allpairs.table.total"=cov.allpairs.table.total, 
                  "cov.perpair.table.total"=cov.perpair.table.total, 
                  "postpop.master.long"=postpop.master.long))
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
                              lowP_to_zeros=lowP_to_zeros,
                              lowDP_to_zeros=lowDP_to_zeros,
                              include_called_allpops=include_called_allpops,
                              min.pairwise.prop=min.pairwise.prop,
                              assessPool_thinning=assessPool_thinning)
    
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
                              lowP_to_zeros=lowP_to_zeros,
                              lowDP_to_zeros=lowDP_to_zeros,
                              include_called_allpops=include_called_allpops,
                              min.pairwise.prop=min.pairwise.prop,
                              assessPool_thinning=assessPool_thinning)
    
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
    
    #postpop.master.long_popoolation$pair <<- gsub(paste0(project_name,   "_"), "", postpop.master.wide_pp$pair) 
    postpop.master.long_popoolation <<- separate(postpop.master.long_popoolation, pair, c("popA","popB"), sep="_", remove=F)
    
    cov.perpair.table.total_popoolation <<- cov.perpair.table.total_popoolation[order(cov.perpair.table.total_popoolation$MinCoverage),]
    cov.perpair.table.total_popoolation <<- separate(cov.perpair.table.total_popoolation, pair, c("popA","popB"), sep="_", remove=F)
    cov.perpair.table.total_popoolation$NumSNPs <<- as.numeric(cov.perpair.table.total_popoolation$NumSNPs)
    cov.perpair.table.total_popoolation$NumContigs <<- as.numeric(cov.perpair.table.total_popoolation$NumContigs)
    cov.perpair.table.total_popoolation$pair <<- as.character(cov.perpair.table.total_popoolation$pair)
    
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
    
    
    #postpop.master.long_poolfstat$pair <<- gsub(paste0(project_name,   "_"), "", postpop.master.wide_pp$pair) 
    postpop.master.long_poolfstat <<- separate(postpop.master.long_poolfstat, pair, c("popA","popB"), sep="_", remove=F)
    
    cov.perpair.table.total_poolfstat <<- cov.perpair.table.total_poolfstat[order(cov.perpair.table.total_poolfstat$MinCoverage),]
    cov.perpair.table.total_poolfstat <<- separate(cov.perpair.table.total_poolfstat, pair, c("popA","popB"), sep="_", remove=F)
    cov.perpair.table.total_poolfstat$NumSNPs <<- as.numeric(cov.perpair.table.total_poolfstat$NumSNPs)
    cov.perpair.table.total_poolfstat$NumContigs <<- as.numeric(cov.perpair.table.total_poolfstat$NumContigs)
    cov.perpair.table.total_poolfstat$pair <<- as.character(cov.perpair.table.total_poolfstat$pair)

        }
}







