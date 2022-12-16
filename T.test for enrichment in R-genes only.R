library(ggplot2)
library(data.table)
library(grid)
library(readr)

# T-Test - is there a significant difference in the average proportion of coverage of  
# each gene region by a particular modification between R-genes at each expression level?
# Plots comparing the average proportion of coverage of each gene region by a particular 
# modification in R-genes at each expression level.

levelsIncluded <- c()

geneCount <- as.data.frame(read_csv(paste("Data\\", analysis, "\\", tissue, "\\Gene count.txt", sep = "")))
geneCount <- geneCount[,-1]

NoExSampleSize <- geneCount[grepl("NLR", geneCount$GeneSet) & grepl("No Expression", geneCount$GeneSet) &
                              !grepl("luster", geneCount$GeneSet), "GeneCount"]
if (NoExSampleSize > 10) {levelsIncluded <- append(levelsIncluded, "No Expression")}

LowExSampleSize <- geneCount[grepl("NLR", geneCount$GeneSet) & grepl("Low Expression", geneCount$GeneSet) &
                               !grepl("luster", geneCount$GeneSet), "GeneCount"]
if (LowExSampleSize > 10) {levelsIncluded <- append(levelsIncluded, "Low Expression")}

InterExSampleSize <- geneCount[grepl("NLR", geneCount$GeneSet) & grepl("Intermediate Expression", geneCount$GeneSet) &
                                 !grepl("luster", geneCount$GeneSet), "GeneCount"]
if (InterExSampleSize > 10) {levelsIncluded <- append(levelsIncluded, "Intermediate Expression")}

HighExSampleSize <- geneCount[grepl("NLR", geneCount$GeneSet) & grepl("High Expression", geneCount$GeneSet) &
                                !grepl("luster", geneCount$GeneSet), "GeneCount"]
if (HighExSampleSize > 10) {levelsIncluded <- append(levelsIncluded, "High Expression")}

df <- df[c(which(df$Expression %in% levelsIncluded)),]

if (length(levelsIncluded) == 2) {
  for (mod in unique(df$Modification)) {
    df1 <- df[df$Modification==mod,]
    
    plot <- ggplot(df1, aes(x = axisGroup, y = Proportion, color = Expression)) +
      geom_line(aes(group = Expression),linewidth = 1) +
      scale_x_continuous(limits = c(-60, 140), breaks = seq(-60, 140, 20), labels = axisText) +
      scale_y_continuous(limits = c(0,1), expand = c(0,0)) + 
      labs(x = "", y = "Average proportion of gene region", color = "Expression level", title = mod) +
      geom_vline(xintercept=0, color="grey", linewidth=1) + theme_minimal() +
      scale_colour_manual(values = c("darkorange3", "darkslategray4"), labels = c(paste("Low Expression (n = ", LowExSampleSize, ")", sep = ""),
                                                                                  paste("No Expression (n = ", NoExSampleSize, ")", sep = ""))) +
      coord_cartesian(ylim= c(0,1), clip = "off") + theme(plot.margin = unit(c(1,1,.3,1), "lines")) +
      annotation_custom(textGrob("% of gene length from TSS", gp=gpar(fontsize=16, col = "grey33")),xmin=0,xmax=100,ymin=-.17,ymax=-.17) + 
      annotation_custom(textGrob("Gene region", gp=gpar(fontsize=16)),xmin=0,xmax=100,ymin=-.25,ymax=-.25) +
      theme(axis.text.x = element_text(size = 14, colour = "black", angle = 45, vjust = 1, hjust = 1), axis.text.y = element_text(size = 14,colour = "black"), 
            axis.title.y = element_text(size = 16, vjust = 2), plot.title = element_text(hjust = .5, size = 16),
            legend.text = element_text(size = 12), legend.title = element_text(size = 14), axis.line = element_line(linewidth = .6))
    
    ggsave(paste("Graphs\\Enrichment\\", analysis, "\\", tissue, "\\R-gene enrichment for ", mod, ".pdf", sep = ""), plot = plot, width = 12, height = 6)
    
  }
} else if (length(levelsIncluded) == 3) {
  for (mod in unique(df$Modification)) {
    df1 <- df[df$Modification==mod,]
    
    geneCountThreshold <- c()
    
    for (level in exLevel) {
      if (nrow(resultsAverageProportions[[analysis]][grepl("NLRs", resultsAverageProportions[[analysis]]$Tissue) &
                                                     grepl(tissue, resultsAverageProportions[[analysis]]$Tissue) &
                                                     !grepl("luster", resultsAverageProportions[[analysis]]$Tissue) &
                                                     grepl(level, resultsAverageProportions[[analysis]]$Expression),]) <= 1) {
        geneCountThreshold <- append(geneCountThreshold, TRUE)
      } else geneCountThreshold <- append(geneCountThreshold, FALSE)
    }
    
    levelsIncluded <- exLevel[c(which(geneCountThreshold==FALSE))]
    
    if (length(levelsIncluded) == 0) {
      df <- resultsAverageProportions[[analysis]][grepl("NLRs", resultsAverageProportions[[analysis]]$Tissue) &
                                                    grepl(tissue, resultsAverageProportions[[analysis]]$Tissue) &
                                                    !grepl("luster", resultsAverageProportions[[analysis]]$Tissue),]
    } else if (length(levelsIncluded) != 0) {
      df <- resultsAverageProportions[[analysis]][grepl("NLRs", resultsAverageProportions[[analysis]]$Tissue) &
                                                    grepl(tissue, resultsAverageProportions[[analysis]]$Tissue) &
                                                    !grepl("luster", resultsAverageProportions[[analysis]]$Tissue),]
      df <- df[c(which(df$Expression %in% levelsIncluded)),]
    }
    
    plot <- ggplot(df1, aes(x = axisGroup, y = Proportion, color = Expression)) +
      geom_line(aes(group = Expression),linewidth = 1) +
      scale_x_continuous(limits = c(-60, 140), breaks = seq(-60, 140, 20), labels = axisText) +
      scale_y_continuous(limits = c(0,1), expand = c(0,0)) + 
      labs(x = "", y = "Average proportion of gene region", color = "Expression level", title = mod) +
      geom_vline(xintercept=0, color="grey", linewidth=1) + theme_minimal() +
      scale_colour_manual(values = c("darkolivegreen", "darkorange3", "darkslategray4"), labels = c(paste("Intermediate Expression (n = ", InterExSampleSize, ")", sep = ""),
                                                                                                    paste("Low Expression (n = ", LowExSampleSize, ")", sep = ""),
                                                                                                    paste("No Expression (n = ", NoExSampleSize, ")", sep = ""))) +
      coord_cartesian(ylim= c(0,1), clip = "off") + theme(plot.margin = unit(c(1,1,.3,1), "lines")) +
      annotation_custom(textGrob("% of gene length from TSS", gp=gpar(fontsize=16, col = "grey33")),xmin=0,xmax=100,ymin=-.17,ymax=-.17) + 
      annotation_custom(textGrob("Gene region", gp=gpar(fontsize=16)),xmin=0,xmax=100,ymin=-.25,ymax=-.25) +
      theme(axis.text.x = element_text(size = 14, colour = "black", angle = 45, vjust = 1, hjust = 1), axis.text.y = element_text(size = 14,colour = "black"), 
            axis.title.y = element_text(size = 16, vjust = 2), plot.title = element_text(hjust = .5, size = 16),
            legend.text = element_text(size = 12), legend.title = element_text(size = 14), axis.line = element_line(linewidth = .6))
    
    ggsave(paste("Graphs\\Enrichment\\", analysis, "\\", tissue, "\\R-gene enrichment for ", mod, ".pdf", sep = ""), plot = plot, width = 12, height = 6)
    
  }
} else if (length(levelsIncluded) == 4) {
  for (mod in unique(df$Modification)) {
    df1 <- df[df$Modification==mod,]
    
    plot <- ggplot(df1, aes(x = axisGroup, y = Proportion, color = Expression)) +
      geom_line(aes(group = Expression),linewidth = 1) +
      scale_x_continuous(limits = c(-60, 140), breaks = seq(-60, 140, 20), labels = axisText) +
      scale_y_continuous(limits = c(0,1), expand = c(0,0)) + 
      labs(x = "", y = "Average proportion of gene region", color = "Expression level", title = mod) +
      geom_vline(xintercept=0, color="grey", linewidth=1) + theme_minimal() +
      scale_colour_manual(values = c("deeppink4","darkolivegreen", "darkorange3", "darkslategray4"), labels = c(paste("High Expression (n = ", HighExSampleSize, ")", sep = ""),
                                                                                                                paste("Intermediate Expression (n = ", InterExSampleSize, ")", sep = ""),
                                                                                                                paste("Low Expression (n = ", LowExSampleSize, ")", sep = ""),
                                                                                                                paste("No Expression (n = ", NoExSampleSize, ")", sep = ""))) +
      coord_cartesian(ylim= c(0,1), clip = "off") + theme(plot.margin = unit(c(1,1,.3,1), "lines")) +
      annotation_custom(textGrob("% of gene length from TSS", gp=gpar(fontsize=16, col = "grey33")),xmin=0,xmax=100,ymin=-.17,ymax=-.17) + 
      annotation_custom(textGrob("Gene region", gp=gpar(fontsize=16)),xmin=0,xmax=100,ymin=-.25,ymax=-.25) +
      theme(axis.text.x = element_text(size = 14, colour = "black", angle = 45, vjust = 1, hjust = 1), axis.text.y = element_text(size = 14,colour = "black"), 
            axis.title.y = element_text(size = 16, vjust = 2), plot.title = element_text(hjust = .5, size = 16),
            legend.text = element_text(size = 12), legend.title = element_text(size = 14), axis.line = element_line(linewidth = .6))
    
    ggsave(paste("Graphs\\Enrichment\\", analysis, "\\", tissue, "\\R-gene enrichment for ", mod, ".pdf", sep = ""), plot = plot, width = 12, height = 6)
    
  }
}