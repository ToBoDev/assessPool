allpairs_summary_plot <-
  function(cov.allpairs.table.total,
           cov.allpairs.table.allpools,
           title.text) {
    
    if(exists("cov.allpairs.table.allpools")){
      
      allpairs_summary_plot <-
        
        plot_ly(cov.allpairs.table.total,
          x = ~ MinCoverage,
          mode = "lines+markers",
          alpha = 0.75) %>%
        
        add_trace(y = ~ NumSNPs,
          name = "NumSNPs",
          mode = "lines+markers") %>%
        
        add_trace(y = ~ NumContigs,
          name = "NumContigs",
          visible = F,
          mode = "lines+markers") %>%
        
        add_trace(y = ~ MeanSNPsPerContig,
          name = "MeanSNPsPerContig",
          visible = F,
          mode = "lines+markers") %>%
        
        add_trace(y = ~ MeanFST,
          name = "MeanFST",
          visible = F,
          mode = "lines+markers") %>%
        
        add_trace(y = ~ SdFST,
          name = "SdFST",
          visible = F,
          mode = "lines+markers")
      
      all.pairs.total  <-   list(
        type = "buttons",
        direction = "right",
        xanchor = 'center',
        yanchor = "top",
        pad = list('r' = 0, 't' = 10, 'b' = 10),
        x = 0.5,
        y = 1.27,
        buttons = list(
          list(method = "update",
            args = list(list(visible = list(
              T, F, F, F, F, F, F, F, F, F
            )),
            list(yaxis = list(title = "Number of SNPs"))),
            label = "NumSNPs"),
          
          list(method = "update",
            args = list(list(visible = list(
              F, T, F, F, F, F, F, F, F, F
            )),
            list(yaxis = list(title = "Number of Contigs"))),
            label = "NumContigs"),
          
          list(method = "update",
            args = list(list(visible = list(
              F, F, T, F, F, F, F, F, F
            )),
            list(yaxis = list(title = "Mean SNPs Per Contig"))),
            label = "MeanSNPsPerContig"),
          
          list(method = "update",
            args = list(list(visible = list(
              F, F, F, T, F, F, F, F, F
            )),
            list(yaxis = list(title = "MeanFST"))),
            label = "MeanFST"),
          
          list(method = "update",
            args = list(list(visible = list(
              F, F, F, F, T, F, F, F, F
            )),
            list(yaxis = list(title = "Standard Deviation FST"))),
            label = "SdFST")
        )
      )
      
      annot <-
        list(
          list(
            text = "Total",
            x = .1,
            y = 1.22,
            xref = 'paper',
            yref = 'paper',
            showarrow = FALSE
          )
        )
      
      allpairs_summary_plot <- allpairs_summary_plot %>% layout(
        yaxis = list(title = "Num SNPs"),
        xaxis = list(title = "Minimun Coverage"),
        title = list(text=title.text, yref= 'paper', y=1),
        showlegend = F,
        updatemenus = list(all.pairs.total),
        annotations = annot
      )
      allpairs_summary_plot
      
      
    } else {
      
      allpairs_summary_plot <-
        plot_ly(
          cov.allpairs.table.total,
          x = ~ MinCoverage,
          mode = "lines+markers",
          alpha = 0.75
        ) %>%
        add_trace(y = ~ NumSNPs,
                  name = "NumSNPs",
                  mode = "lines+markers") %>%
        add_trace(
          y = ~ NumContigs,
          name = "NumContigs",
          visible = F,
          mode = "lines+markers"
        ) %>%
        add_trace(
          y = ~ MeanSNPsPerContig,
          name = "MeanSNPsPerContig",
          visible = F,
          mode = "lines+markers"
        ) %>%
        add_trace(
          y = ~ MeanFST,
          name = "MeanFST",
          visible = F,
          mode = "lines+markers"
        ) %>%
        add_trace(
          y = ~ SdFST,
          name = "SdFST",
          visible = F,
          mode = "lines+markers"
        ) %>%
        add_trace(
          data = cov.allpairs.table.allpools,
          x = ~ MinCoverage,
          y = ~ NumSNPs,
          name = "NumSNPs",
          visible = F,
          mode = "lines+markers"
        ) %>%
        add_trace(
          data = cov.allpairs.table.allpools,
          x = ~ MinCoverage,
          y = ~ NumContigs,
          name = "NumContigs",
          visible = F,
          mode = "lines+markers"
        ) %>%
        add_trace(
          data = cov.allpairs.table.allpools,
          x = ~ MinCoverage,
          y = ~ MeanSNPsPerContig,
          name = "MeanSNPsPerContig",
          visible = F,
          mode = "lines+markers"
        ) %>%
        add_trace(
          data = cov.allpairs.table.allpools,
          x = ~ MinCoverage,
          y = ~ MeanFST,
          name = "MeanFST",
          visible = F,
          mode = "lines+markers"
        ) %>%
        add_trace(
          data = cov.allpairs.table.allpools,
          x = ~ MinCoverage,
          y = ~ SdFST,
          name = "SdFST",
          visible = F,
          mode = "lines+markers"
        )
      
      all.pairs.total  <-   list(
        type = "buttons",
        direction = "right",
        xanchor = 'center',
        yanchor = "top",
        pad = list('r' = 0, 't' = 10, 'b' = 10),
        x = 0.5,
        y = 1.27,
        buttons = list(
          list(
            method = "update",
            args = list(list(visible = list(
              T, F, F, F, F, F, F, F, F, F
            )),
            list(yaxis = list(title = "Number of SNPs"))),
            label = "NumSNPs"
          ),
          list(
            method = "update",
            args = list(list(visible = list(
              F, T, F, F, F, F, F, F, F, F
            )),
            list(yaxis = list(title = "Number of Contigs"))),
            label = "NumContigs"
          ),
          list(
            method = "update",
            args = list(list(visible = list(
              F, F, T, F, F, F, F, F, F
            )),
            list(yaxis = list(title = "Mean SNPs Per Contig"))),
            label = "MeanSNPsPerContig"
          ),
          list(
            method = "update",
            args = list(list(visible = list(
              F, F, F, T, F, F, F, F, F
            )),
            list(yaxis = list(title = "MeanFST"))),
            label = "MeanFST"
          ),
          list(
            method = "update",
            args = list(list(visible = list(
              F, F, F, F, T, F, F, F, F
            )),
            list(yaxis = list(title = "Standard Deviation FST"))),
            label = "SdFST"
          )
        )
      )
      all.pairs.all.pools  <-   list(
        type = "buttons",
        direction = "right",
        xanchor = 'center',
        yanchor = "top",
        pad = list('r' = 0, 't' = 10, 'b' = 10),
        x = 0.5,
        y = 1.17,
        buttons = list(
          list(
            method = "update",
            args = list(list(visible = list(
              F, F, F, F, F, T, F, F, F, F
            )),
            list(yaxis = list(title = "Number of SNPs"))),
            label = "NumSNPS"
          ),
          list(
            method = "update",
            args = list(list(visible = list(
              F, F, F, F, F, F, T, F, F, F
            )),
            list(yaxis = list(title = "Number of Contigs"))),
            label = "NumContigs"
          ),
          list(
            method = "update",
            args = list(list(visible = list(
              F, F, F, F, F, F, F, T, F, F
            )),
            list(yaxis = list(title = "Mean SNPs Per Contig"))),
            label = "MeanSNPsPerContig"
          ),
          list(
            method = "update",
            args = list(list(visible = list(
              F, F, F, F, F, F, F, F, T, F
            )),
            list(yaxis = list(title = "Mean FST"))),
            label = "MeanFST"
          ),
          list(
            method = "update",
            args = list(list(visible = list(
              F, F, F, F, F, F, F, F, F, T
            )),
            list(yaxis = list(title = "Standard Deviation FST"))),
            label = "SdFST"
          )
        )
      )
      
      annot <-
        list(
          list(
            text = "Total",
            x = .1,
            y = 1.22,
            xref = 'paper',
            yref = 'paper',
            showarrow = FALSE
          ),
          list(
            text = "All Pools",
            x = 0.09,
            y = 1.12,
            xref = 'paper',
            yref = 'paper',
            showarrow = FALSE
          )
        )
      
      allpairs_summary_plot <- allpairs_summary_plot %>% layout(
        yaxis = list(title = "Num SNPs"),
        xaxis = list(title = "Minimun Coverage"),
        title = list(text=title.text, yref= 'paper', y=1),
        showlegend = F,
        updatemenus = list(all.pairs.total, all.pairs.all.pools),
        annotations = annot
      )
      allpairs_summary_plot
    }

  }

all_loci_distribution_plot <- function(data_in, title.text){
  data_in %>%
    plot_ly(
      x = ~.fst, 
      y = ~pair, 
      color = ~popA,
      frame = ~MinCoverage, 
      #text = ~pair, 
      hoverinfo = ~.fst,
      type = 'box'
      #mode = 'markers'
    ) %>%
    layout(margin=list(l = 200, r = 20, b = 80, t = 40),
           yaxis = list(tickangle=-35,tickfont=list(size=18)),
           xaxis = list(tickfont=list(size=18)),
           title = list(text=title.text, y=0.95))
}


pairwise_fst_by_cov <- function(data_in, title.text){
  data_in %>%
    plot_ly(
      x = ~MeanFST, 
      y = ~pair, 
      color = ~NumSNPs,
      frame = ~MinCoverage, 
      text = ~NumContigs, 
      hoverinfo = ~pair,
      type = 'scatter',
      mode = 'markers',
      marker = list(size = 12)
    ) %>%
    layout(margin=list(l = 100, r = 20, b = 10, t = 30),
           yaxis = list(tickangle=-35,tickfont=list(size=12)),
           xaxis = list(tickfont=list(size=12)),
           title = list(text=title.text, y=0.95))
}
