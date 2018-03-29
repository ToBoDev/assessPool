#cessPOOL_syncSetup.R
#Split by site types (biallelic, multiallelic, intertions, deletions)

syncSetup <-function(project_name, as.st, POPS, popcomb, include.multiallelic, include.indels, perform_snpfreq, perform_fst, perform_fet, min_count, min_cov, max_cov, min_covered_fract, window_size, step_size, pool_size) { 
    
    #clear out old sync files if they exist
    file.remove(dir(  
      paste(working_dir, project_name, "popoolation", sep="/"), 
      pattern = "*.sync", 
      full.names = TRUE
    ))
  
    #######################Biallelic SNPS###################################
  
    #create SYNC columns by allele (biallelic)
    extract.allele.ba.columns <- function(dataframe, allele_list){
      for (allele in allele_list){
        dataframe[allele] <- ifelse(dataframe$REF==allele, dataframe$RO, ifelse(dataframe$ALT==allele, dataframe$AO, 0))
      }
      dataframe #returns dataframe
    }
    
    #Biallelic SNPS
    as.snp <- as.st[which(as.st$TYPE=="bi-snp" & as.st$AN<=1 & as.st$Rlen<=1),]
    
    if (nrow(as.snp) != 0){
      as.snp[is.na(as.snp)] <- 0 #AO must be integer for this to work
      
      #Biallelic SNPS: create SYNC columns A:T:C:G:N:DEL
      as.snp <- extract.allele.ba.columns(as.snp, c("A","T","C","G","N","DEL"))
      as.snp <- within(as.snp,  sync <- paste(A,T,C,G,N,DEL, sep=":")) #Concatenate into single SYNC column
      as.sync.st <- as.snp
    }
    
    m <- paste("Found ", length(unique(as.snp$snpid)), " biallelic SNPs.", sep="")
    message(m); write.log(paste(m,"\n"), paste(working_dir, project_name, "logs/setup.log", sep="/"))
    
    #######################Multiallelic SNPS####################################
    extract.allele.ma.columns <- function(dataframe, allele_list){
      for (allele in allele_list){
        #dataframe[allele] <- ifelse(dataframe$REF==allele, dataframe$RO, ifelse(dataframe$ALT==allele, dataframe$AO, 0))
        dataframe[allele] <- ifelse(dataframe$REF==allele, dataframe$RO, 
                                    ifelse(dataframe$alt1==allele, dataframe$AO.1, 
                                           ifelse(dataframe$alt2==allele, dataframe$AO.2, 
                                                  ifelse(dataframe$alt3==allele, dataframe$AO.3, 0))))
      }
      dataframe #returns dataframe
    }
    
    #Subset MultAllele SNPs (>1 ALT allele)
    as.mst <- as.st[which(as.st$Rlen==1 & as.st$AN>1),] 
    
    if (nrow(as.mst) != 0){
      #split up alternate alleles
      as.mst$ALT <- strsplit(as.mst$ALT, split="\\,")
      as.mst$AO <- strsplit(as.mst$AO, split="\\,")
      as.mst$alt1 <- NA; as.mst$alt2 <- NA; as.mst$alt3 <- NA
      as.mst$AO.1 <- NA; as.mst$AO.2 <- NA; as.mst$AO.3 <- NA
      
      for (x in 1:nrow(as.mst)){
        as.mst$alt1[x] <- as.mst$ALT[x][[1]][1]
        as.mst$alt2[x] <- as.mst$ALT[x][[1]][2]
        as.mst$AO.1[x] <- as.mst$AO[x][[1]][1]
        as.mst$AO.2[x] <- as.mst$AO[x][[1]][2]
        
        if (as.mst$AN[x]==3){
          as.mst$alt3[x] <- as.mst$ALT[x][[1]][3]
          as.mst$AO.3[x] <- as.mst$AO[x][[1]][3]
        }
      }
      
      as.mst[is.na(as.mst)] <- 0   #AO must be integer for this to work
      
      #extract multiallele columns
      as.mst <- extract.allele.ma.columns(as.mst, c("A","T","C","G","N","DEL"))
      
      #Concatenate into single SYNC column
      as.mst <- within(as.mst,  sync <- paste(A,T,C,G,N,DEL, sep=":"))
      as.mst.merge <- as.mst[,-which(names(as.mst) %in% c("alt1","alt2","alt3","AO.1","AO.2","AO.3"))]
      as.sync.st <- rbind(as.sync.st, as.mst.merge)
    }
  
    m <- paste("Found ", length(unique(as.mst$snpid)), " multiallelic SNPs.", sep="")
    message(m); write.log(paste(m,"\n"), paste(working_dir, project_name, "logs/setup.log", sep="/"))
    
    
    #Popoolation's SYNC format does not have a call for Insertions - so we need to code them as Deletions.
    #So the ALT becomes the REF and then its a DELETION
    #######################Deletions####################################
    extract.allele.del.columns <- function(dataframe, ref_col, allele_list){
      for (allele in allele_list){
        dataframe[allele] <- ifelse(dataframe[ref_col]==allele, dataframe$RO,0)
      }
      dataframe #returns dataframe
    }
    
    as.del <- as.st[which(as.st$TYPE=="del" & as.st$DEL.len<=1),] #restrict to single base deletions
    
    if (nrow(as.del)!=0){
      as.del$Ref.del <- stri_sub(as.del$REF, from=2, to=2)  #reference base that is deleted (2nd base of REF) 
      as.del$AO <- as.integer(as.del$AO); as.del$AO <- as.integer(as.del$AO)
      as.del[is.na(as.del)] <- 0   #AO must be integer for this to work
      #DEL: create SYNC columns A:T:C:G:N:DEL
      as.del <- extract.allele.del.columns(as.del, "Ref.del", c("A","T","C","G","N","DEL"))
      as.del$DEL <- as.del$AO
      #Concatenate into single SYNC column
      as.del <- within(as.del,  sync <- paste(A,T,C,G,N,DEL, sep=":"))
      as.del.merge <- as.del[,-which(names(as.del) %in% c("Ref.del"))]
      as.sync.st <- rbind(as.sync.st, as.del.merge)
      
    }
    
    m <- paste("Found ", length(unique(as.del$snpid)), " single-base deletions.", sep="")
    message(m); write.log(paste(m,"\n"), paste(working_dir, project_name, "logs/setup.log", sep="/"))
    
    #######################Insertions####################################
    as.ins <- as.st[which(as.st$TYPE=="ins" & as.st$INS.len<=1),] #restrict to single base insertions
    
    if (nrow(as.ins) != 0){
      as.ins$Alt.ins <- stri_sub(as.ins$ALT, from=2, to=2)  #ALT allele that is inserted (2nd base of ALT)
      as.ins[is.na(as.ins)] <- 0   #AO must be integer for this to work
      
      #because this is an insertion, need to switch alt and ref
      tmp1 <- as.ins$REF; tmp2 <- as.ins$RO
      as.ins$REF <- as.ins$ALT; as.ins$RO <- as.ins$AO
      as.ins$ALT <- tmp1; as.ins$AO <- tmp2
      rm(tmp1, tmp2)
      
      # INSERTION: create SYNC columns A:T:C:G:N:DEL
      as.ins <- extract.allele.del.columns(as.ins, "Alt.ins", c("A","T","C","G","N","DEL"))
      as.ins$DEL <- as.ins$AO
      #Concatenate into single SYNC column
      as.ins <- within(as.ins,  sync <- paste(A,T,C,G,N,DEL, sep=":"))
      as.ins.merge <- as.ins[,-which(names(as.ins) %in% c("Alt.ins"))]
      as.sync.st <- rbind(as.sync.st, as.ins.merge)
    }
    
    m <- paste("Found ", length(unique(as.ins$snpid)), " single-base insertions.", sep="")
    message(m); write.log(paste(m,"\n"), paste(working_dir, project_name, "logs/setup.log", sep="/"))
    
    ###########################################################
    
    #JOIN rows
    #as.sync.st <- rbind(as.snp, as.mst.merge, as.del.merge, as.ins.merge)
    
    #Now, unstack merged file by population (long -> wide)
    as.sync.st <- as.sync.st[,which(names(as.sync.st) %in% c("snpid","CHROM","POS","REF","sync","Pop","Rlen","TYPE"))]
    colnames(as.sync.st)[which(colnames(as.sync.st)=="Pop")] <- "sync."
    as.sync <- spread(as.sync.st, sync., sync, sep="")
    
    #Export All Variable Sites
    as.sync.export <- as.sync[,-which(names(as.sync) %in% c("snpid","Rlen","TYPE"))]
    write.table(as.sync.export,paste(working_dir, "/", project_name, "/output/", project_name, "_", length(POPS),"pools_allvar.txt",sep=""),sep="\t",row.names=F, col.names=T, quote=F)
    m <- paste("Exported all variable sites to output/", project_name, "_", length(POPS),"pools_allvar.txt", sep="")
    message(m); write.log(paste(m,"\n"), paste(working_dir, project_name, "logs/setup.log", sep="/"))
    
    #Export All SNPs
    as.sync.snps <- as.sync[which(as.sync$TYPE=="bi-snp" | as.sync$TYPE=="tri-snp" | as.sync$TYPE=="quad-snp"),-which(names(as.sync) %in% c("snpid","Rlen","TYPE"))]
    write.table(as.sync.export,paste(working_dir, "/", project_name, "/output/", project_name, "_", length(POPS),"pools_allsnps.txt", sep=""),sep="\t",row.names=F, col.names=T, quote=F)
    m <- paste("Exported all SNPs to output/", project_name, "_", length(POPS),"pools_allsnps.txt", sep="")
    message(m); write.log(paste(m,"\n"), paste(working_dir, project_name, "logs/setup.log", sep="/"))
    
    #Exports only biallelic SNPs
    as.sync.bsnps <- as.sync[which(as.sync$TYPE=="bi-snp"),-which(names(as.sync) %in% c("snpid","Rlen","TYPE"))]
    write.table(as.sync.export,paste(working_dir, "/", project_name, "/output/", project_name, "_", length(POPS),"pools_bisnps.txt", sep=""),sep="\t",row.names=F, col.names=T, quote=F)
    m <- paste("Exported all bi-SNPs to output/", project_name, "_", length(POPS),"pools_bisnps.txt", sep="")
    message(m); write.log(paste(m,"\n"), paste(working_dir, project_name, "logs/setup.log", sep="/"))
    
    #Export INDELS and add TDP to as.indels
    #as.sync.indelsnps <- as.sync[which(as.sync$TYPE=="del" | as.sync$TYPE=="ins" | as.sync$TYPE=="bi-snp"),-which(names(as.sync) %in% c("snpid","Rlen","TYPE"))]
    as.indels <- as.sync[which(as.sync$TYPE=="del" | as.sync$TYPE=="ins"),]
    as.indels <- merge(as.indels, as[,which(names(as) %in% c("TDP","snpid"))], by="snpid", all.x=TRUE, all.y=FALSE )
    write.table(as.indels, paste(working_dir, "/", project_name, "/output/", project_name, "_", length(POPS),"pools_allindels.txt", sep=""), sep="\t",row.names=F, col.names=T, quote=F)
    m <- paste("Exported all INDELs to output/", project_name, "_", length(POPS),"pools_allindels.txt", sep="")
    message(m); write.log(paste(m,"\n"), paste(working_dir, project_name, "logs/setup.log", sep="/"))
    
    if (include.multiallelic & !include.indels) { 
      vcf2popool_sync_out <- as.sync.snps #include multiallelic but not indels
      m <- paste("Writing ",ncol(popcomb)," pairwise .sync files containing all SNPs.\n", sep="")
    } else if (include.indels & !include.multiallelic) {
      vcf2popool_sync_out <- as.sync.indelsnps #include indels but not multiallelic
      m <- paste("Writing ",ncol(popcomb)," pairwise .sync files containing biallelic SNPs and indels.\n", sep="")
    } else if (!include.indels & !include.multiallelic) {
      vcf2popool_sync_out <- as.sync.bsnps #only include biallelic snps
      m <- paste("Writing ",ncol(popcomb)," pairwise .sync files containing biallelic SNPs.\n", sep="")
    } else {
      vcf2popool_sync_out <- as.sync.export #include all variable sites
      m <- paste("Writing ",ncol(popcomb)," pairwise .sync files containing all variable sites.\n", sep="")
    }
    
    message(m); write.log(m, paste(working_dir, project_name, "logs/popoolation2.log", sep="/"))
    
    rm(list=ls(pattern="*as.")[-which(ls(pattern="*as.")=="as.st")])
    colnames(vcf2popool_sync_out) <- c("CHROM","POS", "REF", POPS)
    vcf2popool_sync_out <- arrange(vcf2popool_sync_out, CHROM,POS)
    
    #create pairwise files for all populations
    for(i in 1:ncol(popcomb))  write.table(bind_cols(vcf2popool_sync_out[,c(1:3)], select(vcf2popool_sync_out, popcomb[1,i]), select(vcf2popool_sync_out, popcomb[2,i])), paste(working_dir, "/", project_name, "/popoolation/", project_name, '_', popcomb[1,i], '_', popcomb[2,i], '.sync', sep=""), quote = F, row.names = F, col.names = F)
    
    #define desired tests and parameters
    setwd(paste(working_dir, "/", project_name, "/popoolation", sep=""))
    pop_path <- pop_path <- paste(working_dir,"/scripts/p2/",sep="")
    test <- c("snp-frequency-diff.pl", "fst-sliding.pl", "fisher-test.pl")
    
    #configure parameters to be used by Popoolation2
    min_count <- paste("--min-count", as.character(min_count))
    min_cov <- paste("--min-coverage", as.character(min_cov))
    max_cov <- paste("--max-coverage", as.character(max_cov))
    min_covered_fract <- paste("--min-covered-fraction", as.character(min_covered_fract))
    window_size <- paste("--window-size", as.character(window_size))
    step_size <- paste("--step-size", as.character(step_size))
    pool_size <- paste("--pool-size", as.character(pool_size))
    
    #write parameters to logfile
    write.log("PoPoolation2 Parameters:", paste(working_dir, project_name, "logs/popoolation2.log", sep="/"))
    write.log(min_count, paste(working_dir, project_name, "logs/popoolation2.log", sep="/"))
    write.log(min_cov, paste(working_dir, project_name, "logs/popoolation2.log", sep="/"))
    write.log(max_cov, paste(working_dir, project_name, "logs/popoolation2.log", sep="/"))
    write.log(min_covered_fract, paste(working_dir, project_name, "logs/popoolation2.log", sep="/"))
    write.log(window_size, paste(working_dir, project_name, "logs/popoolation2.log", sep="/"))
    write.log(step_size, paste(working_dir, project_name, "logs/popoolation2.log", sep="/"))
    write.log(pool_size, paste(working_dir, project_name, "logs/popoolation2.log", sep="/"))
    
    write("\n", file="popool2_run.sh", append= F)
    
    if (perform_snpfreq) {
      for(i in 1:ncol(popcomb)){  
        write(paste('perl', paste(pop_path, test[1],sep=''), "--input", paste(project_name, '_', popcomb[1,i], '_', popcomb[2,i], '.sync', sep=""), "--output-prefix", paste(project_name, '_', popcomb[1,i], '_', popcomb[2,i], sep=""), min_count, min_cov, max_cov,  sep=" "), file="popool2_run.sh", append = T)
      }
      write(paste("\n", paste(replicate(100, "#"), collapse=''),"\n"), file="popool2_run.sh", append = T)
      write.log("Performing SNP frequency runs on SYNC files", paste(working_dir, project_name, "logs/popoolation2.log", sep="/"))
    }
    
    #FST scripts
    if (perform_fst){
      for(i in 1:ncol(popcomb)){            
        write(paste('perl', paste(pop_path, test[2],sep=''), "--input", paste(project_name, '_', popcomb[1,i], '_', popcomb[2,i], '.sync', sep=""), "--output", paste(project_name, '_', popcomb[1,i], '_', popcomb[2,i], '.fst', sep=""), "--suppress-noninformative", min_count, min_cov, max_cov, min_covered_fract, window_size, step_size, pool_size,  sep=" "), file="popool2_run.sh", append = T)
      }
      write(paste("\n",paste(replicate(100, "#"), collapse=''),"\n"), file="popool2_run.sh", append = T)
      write.log("Performing FST runs on SYNC files", paste(working_dir, project_name, "logs/popoolation2.log", sep="/"))
    }
    
    #FET scripts 
    if (perform_fet){
      for(i in 1:ncol(popcomb)){            
        write(paste('perl', paste(pop_path, test[3],sep=''), "--input", paste(project_name, '_', popcomb[1,i], '_', popcomb[2,i], '.sync', sep=""), "--output", paste(project_name, '_', popcomb[1,i], '_', popcomb[2,i], '.fet', sep=""), "--suppress-noninformative", min_count, min_cov, max_cov,  sep=" "),  file="popool2_run.sh", append = T)
      }
      write(paste("\n",paste(replicate(100, "#"), collapse=''),"\n"), file="popool2_run.sh", append = T)
      write.log("Performing FET runs on SYNC files", paste(working_dir, project_name, "logs/popoolation2.log", sep="/"))
    }
    
}
