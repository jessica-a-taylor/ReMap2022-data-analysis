library(ggplot2)
library(data.table)
library(grid)
library(readr)

# For each chromatin modification, plot a bar graph of the enrichment in each R-gene.
dataToUse <- resultsProportions[[analysis]][grepl(tissue, resultsProportions[[analysis]]$dataToAnalyse) & 
                                              grepl("NLR", resultsProportions[[analysis]]$dataToAnalyse) & 
                                              !grepl("luster", resultsProportions[[analysis]]$dataToAnalyse),]


for (mod in c("H3K9me2","H3K27me3","H2A-Z","H2AK121ub","H3K4me3","H3K36me3","H3K27ac","H3K9ac")) {
  df <- dataToUse[dataToUse$Modification==mod,]
  
  for (r in df$Region) {
    df1 <- df[df$Region == r,]
    
    levelsIncluded <- unique(df1$Expression)
    
    if (length(levelsIncluded) == 2) {
      
    plot <- ggplot(df1, aes(x = Gene, y = Proportion, fill = Expression)) + 
      scale_y_continuous(limits = c(0,1), expand = c(0,0)) + 
      scale_fill_manual(values = c("darkorange3","darkslategray4")) +
      geom_bar(stat = "identity", position = "dodge") + theme_minimal() +
      labs(x = "", y = "Average proportion of gene region", title = paste(mod, "-", r, sep = " ")) +
      geom_vline(xintercept=0, color="grey", linewidth=1) +
      coord_cartesian(ylim= c(0,1), clip = "off") + theme(plot.margin = unit(c(1,1,.3,1), "lines"), plot.title = element_text(hjust = 0.5, size = 16),
                                                          axis.text.x = element_text(angle = 90, size = 10, vjust = 0.5), axis.text.y = element_text(size = 14),
                                                          axis.title.y = element_text(size = 16, vjust = 2)) 
      
    ggsave(paste("Graphs\\Enrichment\\", analysis, "\\", tissue, "\\", mod, " Enrichment per gene", ".pdf", sep = ""), plot = plot, width = 24, height = 6)  
    
    } else if (length(levelsIncluded) == 3) {
      
      plot <- ggplot(df1, aes(x = Gene, y = Proportion, fill = Expression)) + 
        scale_y_continuous(limits = c(0,1), expand = c(0,0)) + 
        scale_fill_manual(values = c("darkolivegreen","darkorange3","darkslategray4")) +
        geom_bar(stat = "identity", position = "dodge") + theme_minimal() +
        labs(x = "", y = "Average proportion of gene region", title = paste(mod, "-", r, sep = " ")) +
        geom_vline(xintercept=0, color="grey", linewidth=1) +
        coord_cartesian(ylim= c(0,1), clip = "off") + theme(plot.margin = unit(c(1,1,.3,1), "lines"), plot.title = element_text(hjust = 0.5, size = 16),
                                                            axis.text.x = element_text(angle = 90, size = 10, vjust = 0.5), axis.text.y = element_text(size = 14),
                                                            axis.title.y = element_text(size = 16, vjust = 2)) 
      
      ggsave(paste("Graphs\\Enrichment\\", analysis, "\\", tissue, "\\", mod, " Enrichment per gene", ".pdf", sep = ""), plot = plot, width = 24, height = 6)  
    
    } else if (length(levelsIncluded) == 4) {
      
      plot <- ggplot(df1, aes(x = Gene, y = Proportion, fill = Expression)) + 
        scale_y_continuous(limits = c(0,1), expand = c(0,0)) + 
        scale_fill_manual(values = c("deeppink4","darkolivegreen","darkorange3","darkslategray4")) +
        geom_bar(stat = "identity", position = "dodge") + theme_minimal() +
        labs(x = "", y = "Average proportion of gene region", title = paste(mod, "-", r, sep = " ")) +
        geom_vline(xintercept=0, color="grey", linewidth=1) +
        coord_cartesian(ylim= c(0,1), clip = "off") + theme(plot.margin = unit(c(1,1,.3,1), "lines"), plot.title = element_text(hjust = 0.5, size = 16),
                                                            axis.text.x = element_text(angle = 90, size = 10, vjust = 0.5), axis.text.y = element_text(size = 14),
                                                            axis.title.y = element_text(size = 16, vjust = 2)) 
      
      ggsave(paste("Graphs\\Enrichment\\", analysis, "\\", tissue, "\\", mod, " Enrichment per gene", ".pdf", sep = ""), plot = plot, width = 24, height = 6)  
      
    }
  }
}
