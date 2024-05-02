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

poolHeatmaps <- function(data_in, heatmap_cov, pool_order, fst.c, pools_to_remove){

    #extracts necessary coverage level from parameters given
    #if invalid coverage, throw error
    if (nrow(data_in)==0) {
      stop(cat(paste("Invalid coverage level. Please pick from following list:", paste(covs, collapse=", "), sep="\n")))
    }
    
    #summarize FST/FET

  sum.postpop <- data_in %>% filter(fst.dp > heatmap_cov-1) %>% group_by(pair) %>% select(-MinCoverage) %>% distinct_all() %>% 
    dplyr::summarise(fst=mean(.fst), 
                     ttest=t.test(.fst, mu = 0 , alternative="greater")[3]$p.value,
                     SNPs=n_distinct(snpid),
                     mean.coverage=mean(fst.dp))

    # sum.postpop <- postpop %>% pivot_wider(values_from = value, names_from = analysis) %>% group_by(pair) %>%
    #   dplyr::summarise(Fst=mean(.fst), 
    #             chisq=unlist(sumlog2(.fet))[1], 
    #             df=unlist(sumlog2(.fet)[2]), 
    #             Pvals=unlist(sumlog2(.fet)[3]))
    # 
    # tmp.cols <- colsplit(sum.postpop$pair, paste0(project_name,"_"), c("pn", "pop"))
    tmp.cols <- colsplit(sum.postpop$pair, "_", c("pop1", "pop2")); colnames(tmp.cols) <- c("pop1", "pop2")
    sum.postpop$pop1 <- tmp.cols$pop1; sum.postpop$pop2 <- tmp.cols$pop2; rm(tmp.cols)
    
    #set up FST/Chi-Square matrices
    an <- with(sum.postpop, sort(unique(c(as.character(pop1),as.character(pop2)))))
    fst.matrix <- array(NA, c(length(an), length(an)), list(an, an))
    ttest.matrix <- array(NA, c(length(an), length(an)), list(an, an))
    SNPs.matrix <- array(NA, c(length(an), length(an)), list(an, an))
    cov.matrix <- array(NA, c(length(an), length(an)), list(an, an))
   # pval.matrix <- array(NA, c(length(an), length(an)), list(an, an))
    f <- match(sum.postpop$pop1, an)
    j <- match(sum.postpop$pop2, an)
    fst.matrix[cbind(f,j)] <- fst.matrix[cbind(j,f)] <- sum.postpop$fst
    ttest.matrix[cbind(f,j)] <- ttest.matrix[cbind(j,f)] <- sum.postpop$ttest
    SNPs.matrix[cbind(f,j)] <- SNPs.matrix[cbind(j,f)] <- sum.postpop$SNPs
    cov.matrix[cbind(f,j)] <- cov.matrix[cbind(j,f)] <- sum.postpop$mean.coverage
    
    #export FST/FET matrices
    maxpools <- rownames(which(fst.matrix == max(fst.matrix, na.rm=TRUE), arr.ind = TRUE))
    minpools <- rownames(which(fst.matrix == min(fst.matrix, na.rm=TRUE), arr.ind = TRUE))
    
    fst <- as.matrix(fst.matrix)[pool_order,rev(pool_order)]
    ttest <- as.matrix(ttest.matrix)[pool_order,rev(pool_order)]
    snps <- as.matrix(SNPs.matrix)[pool_order,rev(pool_order)]
    coverage <- as.matrix(cov.matrix)[pool_order,rev(pool_order)]
    
    
    if(!any(pools_to_remove %in% "none")){
      pool.index <- which(colnames(fst) %in% pools_to_remove)
      
      fst <- fst[-pool.index,-pool.index]
      ttest <- ttest[-pool.index,-pool.index]
      snps <- snps[-pool.index,-pool.index]
      coverage <- coverage[-pool.index,-pool.index]
    }
    
    message("Largest mean FST at this coverage is between ", maxpools[1], " and ", maxpools[2], " at ", round(max(fst.matrix, na.rm=TRUE), 4))
    message("Smallest mean FST at this coverage is between ", minpools[1], " and ", minpools[2], " at ", round(min(fst.matrix, na.rm=TRUE), 4))
    
    message("\nExporting FST matrix to ", "output/",project_name,"_Fst.matrix.", heatmap_cov,"x_",fst.c,".csv")
    message("\nExporting t-test matrix to ", "output/",project_name,"_ttest.matrix.", heatmap_cov,"x_",fst.c,".csv")
    message("\nExporting numSNPs matrix to ", "output/",project_name,"_numSNPs.matrix.", heatmap_cov,"x_",fst.c,".csv")
    message("\nExporting coverage matrix to ", "output/",project_name,"_coverage.matrix.", heatmap_cov,"x_",fst.c,".csv")
    #message("\nExporting Chi-Square matrix to ", "output/",project_name,"_pval_sumlog.matrix.", heatmap_cov,"x.csv")
    write.csv(fst, paste0(working_dir, "/", project_name, "/output/",project_name,"_Fst.matrix.", heatmap_cov,"x_",fst.c,".csv")[1])
    write.csv(ttest, paste0(working_dir, "/", project_name, "/output/",project_name,"_ttest.matrix.", heatmap_cov,"x_",fst.c,".csv")[1])
    write.csv(snps, paste0(working_dir, "/", project_name, "/output/",project_name,"_SNPs.matrix.", heatmap_cov,"x_",fst.c,".csv")[1])
    write.csv(coverage, paste0(working_dir, "/", project_name, "/output/",project_name,"_covs.matrix.", heatmap_cov,"x_",fst.c,".csv")[1])
   # write.csv(pval.matrix, paste(working_dir, "/", project_name, "/output/",project_name,"_pval_sumlog.matrix.", heatmap_cov,"x.csv", sep=""))
    
    return(list("fst"=fst, "ttest"=ttest, "snps"=snps, "coverage"=coverage))
  }


generate_heatmap <- function(data_in, metric, fst.c, palette, main_title, theme){

if(fst.c=="popoolation") {sub_title="PoPoolation2"}
if(fst.c=="poolfstat") {sub_title="{poolfstat}"}
  
  
if(theme=="light"){
  main_cols <- colorRamp(natparks.pals(palette))
  text_cols <- "darkslategrey4"
  bg_cols <- c('black', "#f9f5fc")
}
if(theme=="white"){
  main_cols <- colorRamp(natparks.pals(palette))
  text_cols <- "darkslategrey4"
  bg_cols <- c('white', "#white")
}
if(theme=="dark"){
  main_cols <- colorRamp(natparks.pals(palette))
  text_cols <- "#abc7c4"
  bg_cols <- c('#0b1211', "#0b1211")
}

    
#individual interactive plots

if(metric == "fst"){
  #Fst plot
  fst_plot <- plot_ly(z = data_in$fst, colors = main_cols, 
          x = colnames(data_in$fst), y=rownames(data_in$fst), type = "heatmap",
          colorbar = list(title = list(text="<i>F<sub>ST<sub><i>",
                                       font=list(color=text_cols, family = "Times New Roman")),
                          tickfont = list(size = 12, color=text_cols, family = "Times New Roman"),
                          len = 0.75),
          hovertemplate = paste('%{y}','<br>%{x}','<br>Fst:%{z:.3f}')) %>% 
    layout(margin=list(l = 100, r = 20, b = 80, t = 45),
           yaxis = list(tickangle=0, showgrid = F,
                        tickfont=list(size = 16, color = text_cols, family = "Times New Roman")),
           xaxis = list(tickangle=45, showgrid = F,
                        tickfont=list(size = 16, color = text_cols, family = "Times New Roman")),
           title = list(text=paste0("<i>",main_title,"</i>"," (", sub_title, ")"),
                        y=.97,
                        font = list(size = 20, color = text_cols,  family = "Times New Roman")),
           plot_bgcolor = bg_cols[1], 
           paper_bgcolor = bg_cols[2])
  
  return(fst_plot)
}
 
  if(metric == "ttest"){ 
  #ttest
  ttest_plot <- plot_ly(z = data_in$ttest, colors = main_cols, 
          x = colnames(data_in$fst), y=rownames(data_in$fst), type = "heatmap",
          colorbar = list(title = list(text="<i>pvalue<i>",
                                       font=list(color=text_cols, family = "Times New Roman")),
                          tickfont = list(size = 12, color=text_cols, family = "Times New Roman"),
                          len = 0.75)) %>% 
    layout(margin=list(l = 100, r = 20, b = 80, t = 45),
           yaxis = list(tickangle=0, showgrid = F,
                        tickfont=list(size = 16, color = text_cols, family = "Times New Roman")),
           xaxis = list(tickangle=45, showgrid = F,
                        tickfont=list(size = 16, color = text_cols, family = "Times New Roman")),
           title = list(text=paste0("<i>",main_title,"</i>"," (", sub_title, ")"),
                        y=.97,
                        font = list(size = 20, color = text_cols,  family = "Times New Roman")),
           plot_bgcolor = bg_cols[1], 
           paper_bgcolor = bg_cols[2])
  
  return(ttest_plot)
  }
  
  if(metric == "snps"){
  #snps
  snps_plot <- plot_ly(z = data_in$snps, colors = main_cols, 
          x = colnames(data_in$fst), y=rownames(data_in$fst), type = "heatmap",
          colorbar = list(title = list(text="SNPs",
                                       font=list(color=text_cols, family = "Times New Roman")),
                          tickfont = list(size = 12, color=text_cols, family = "Times New Roman"),
                          len = 0.75)) %>% 
    layout(margin=list(l = 100, r = 20, b = 80, t = 45),
           yaxis = list(tickangle=0, showgrid = F,
                        tickfont=list(size = 16, color = text_cols, family = "Times New Roman")),
           xaxis = list(tickangle=45, showgrid = F,
                        tickfont=list(size = 16, color = text_cols, family = "Times New Roman")),
           title = list(text=paste0("<i>",main_title,"</i>"," (", sub_title, ")"),
                        y=.97,
                        font = list(size = 20, color = text_cols,  family = "Times New Roman")),
           plot_bgcolor = bg_cols[1], 
           paper_bgcolor = bg_cols[2])
    
  return(snps_plot)
    
  }
    
    
  if(metric == "coverage"){
  #coverage
  coverage_plot <- plot_ly(z = data_in$coverage, colors = main_cols, 
          x = colnames(data_in$fst), y=rownames(data_in$fst), type = "heatmap",
          colorbar = list(title = list(text="Coverage",
                                       font=list(color=text_cols, family = "Times New Roman")),
                          tickfont = list(size = 12, color=text_cols, family = "Times New Roman"),
                          len = 0.75)) %>% 
    layout(margin=list(l = 100, r = 20, b = 80, t = 45),
           yaxis = list(tickangle=0, showgrid = F,
                        tickfont=list(size = 16, color = text_cols, family = "Times New Roman")),
           xaxis = list(tickangle=45, showgrid = F,
                        tickfont=list(size = 16, color = text_cols, family = "Times New Roman")),
           title = list(text=paste0("<i>",main_title,"</i>"," (", sub_title, ")"),
                        y=.97,
                        font = list(size = 20, color = text_cols,  family = "Times New Roman")),
           plot_bgcolor = bg_cols[1], 
           paper_bgcolor = bg_cols[2])
    
    return(coverage_plot)
  }
}


#static ggplot outputs
 
  # data_in$fst %>% as.data.frame() %>% rownames_to_column("f_id") %>% pivot_longer(-c(f_id), names_to = "samples", values_to = "Fst") %>% 
  #   ggplot(aes(x=samples, y=f_id, fill=Fst)) + 
  #   geom_raster() +
  #   scale_fill_viridis_c()

  


#tiled summary plot
    


#wrapper for plotting
plot_heatmaps <- function(fst.c,include_called_allpops, pools_to_remove){
  
  if(include_called_allpops){
    
    data_in <- get(paste0("postpop.master.long.allpools_",fst.c))
}
  if(include_called_allpops==F){  
      
    data_in <- get(paste0("postpop.master.long_",fst.c))
}

hm_list_out <- poolHeatmaps(data_in=data_in,
                            heatmap_cov=heatmap_cov,
                            pool_order=pool_order,
                            fst.c=fst.c,
                            pools_to_remove=pools_to_remove)

hm.fst <- generate_heatmap(data_in = hm_list_out, metric = "fst", fst.c=fst.c,
                 main_title = paste(spp_name),
                 palette = palette, theme = theme)

hm.ttest <- generate_heatmap(data_in = hm_list_out, metric = "ttest", fst.c=fst.c,
                 main_title = paste(spp_name),
                 palette = palette, theme = theme)

hm.snps <- generate_heatmap(data_in = hm_list_out, metric = "snps", fst.c=fst.c,
                 main_title = paste(spp_name),
                 palette = palette, theme = theme)

hm.cov <- generate_heatmap(data_in = hm_list_out, metric = "coverage", fst.c=fst.c,
                 main_title = paste(spp_name),
                 palette = palette, theme = theme)

return(list( hm.fst,  hm.ttest,  hm.snps, hm.cov))
}
