#define new sumlog function to deal with NAs and 0s
sumlog2 <- function (p) {
  keep <- (p > 0) & (p <= 1)
  lnp <- log(p[keep])
  chisq <- (-2) * sum(lnp)
  df <- 2 * length(lnp)
  if (sum(1L * keep) < 2) 
    k <- NA
  else( k <- 1)
  if (length(lnp) != length(p)) {
    warning("Some studies omitted")
  }
  res <- list(chisq = chisq*k, df = df*k, p = pchisq(chisq, df, lower.tail = FALSE)*k, validp = p[keep])
  class(res) <- c("sumlog2", "metap")
  res
}

poolHeatmaps <- function(heatmap_cov, postpop){
  
  #extracts necessary coverage level from parameters given
  #if invalid coverage, throw error
  pos_cov <- unique(postpop$MinCoverage)
  postpop <- postpop[which(postpop$MinCoverage==heatmap_cov),]
  if (nrow(postpop)==0) {
    stop(paste("Invalid coverage level. Please pick from following list:\n", paste(pos_cov, collapse=", "), sep=""))
  }
  
  #summarize FST/FET
  sum.postpop <- spread(postpop, analysis, value) %>% group_by(pair) %>% summarise(Fst=mean(.fst), chisq=unlist(sumlog2(.fet)[1]), df=unlist(sumlog(.fet)[2]), Pvals=unlist(sumlog(.fet)[3]))
  tmp.cols <- colsplit(sum.postpop$pair, project_name, c("pn", "pop"))[,-1]
  tmp.cols <- colsplit(tmp.cols, "_", c("pop1", "pop2"))[,-1]; colnames(tmp.cols) <- c("pop1", "pop2")
  sum.postpop$pop1 <- tmp.cols$pop1; sum.postpop$pop2 <- tmp.cols$pop2; rm(tmp.cols)
  
  #set up FST/Chi-Square matrices
  an <- with(sum.postpop, sort(unique(c(as.character(pop1),as.character(pop2)))))
  fst.matrix <- array(NA, c(length(an), length(an)), list(an, an))
  chsq.matrix <- array(NA, c(length(an), length(an)), list(an, an))
  f <- match(sum.postpop$pop1, an)
  j <- match(sum.postpop$pop2, an)
  fst.matrix[cbind(f,j)] <- fst.matrix[cbind(j,f)] <- sum.postpop$Fst
  chsq.matrix[cbind(f,j)] <- chsq.matrix[cbind(j,f)] <- sum.postpop$chisq
  
  #export FST/FET matrices
  maxpools <- rownames(which(fst.matrix == max(fst.matrix, na.rm=TRUE), arr.ind = TRUE))
  minpools <- rownames(which(fst.matrix == min(fst.matrix, na.rm=TRUE), arr.ind = TRUE))
  
  message("Largest mean FST at this coverage is between ", maxpools[1], " and ", maxpools[2], " at ", round(max(fst.matrix, na.rm=TRUE), 4))
  message("Smallest mean FST at this coverage is between ", minpools[1], " and ", minpools[2], " at ", round(min(fst.matrix, na.rm=TRUE), 4))
  
  message("\nExporting FST matrix to ", "output/",project_name,"_Fst.matrix.", heatmap_cov,"x.csv")
  message("\nExporting Chi-Square matrix to ", "output/",project_name,"_Fet.matrix.", heatmap_cov,"x.csv")
  write.csv(fst.matrix, paste(working_dir, "/", project_name, "/output/",project_name,"_Fst.matrix.", heatmap_cov,"x.csv", sep=""))
  write.csv(chsq.matrix, paste(working_dir, "/", project_name, "/output/",project_name,"_Fet.matrix.", heatmap_cov,"x.csv", sep=""))
  
  return(list("fst"=fst.matrix, "chisq"=chsq.matrix))
}