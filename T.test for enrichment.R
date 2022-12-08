library(ggplot2)
library(data.table)
library(grid)

# Wilcox & Kolmogorov-Smirnov Tests - is there a significant difference in the average proportion of coverage of  
# each gene region by a particular modification between R-genes and controls?

# Compare R-genes to each control gene set.
dataToUse <- allResultsProportions[grepl(tissue, allResultsProportions$SampleGenes) & 
                              !grepl("luster", allResultsProportions$SampleGenes),]

t.testDF <- data.frame()


for (mod in unique(allResultsProportions$Modification)) {
  df <- dataToUse[dataToUse$Modification==mod,]
  
  for (r in unique(allResultsProportions$Region)) {
    df1 <- df[df$Region==r,]
    
    for (level in c("No Expression", "Low Expression")) {
      df2 <- df1[df1$Expression == level,]
      
      statTest <- t.test(df2[c(which(df2$SampleGenes == paste("NLRs_", tissue, sep = ""))), "Proportion"], mu = mean(df2[grepl("control", df2$SampleGenes) & 
                                                                                                                           grepl(tissue, df2$SampleGenes), "Proportion"]))
      
      if (is.na(statTest$p.value)==TRUE) {
        t.testDF <- t.testDF
      } else if (statTest$p.value <= 0.05) {
        t.testDF <- rbind(t.testDF, data.frame(Modification = mod,
                                               Region = r,
                                               Expression = level,
                                               W.statistic = statTest$statistic,
                                               p.value = statTest$p.value,
                                               Significance = "*"))
      } else  t.testDF <- rbind(t.testDF, data.frame(Modification = mod,
                                                     Region = r,
                                                     Expression = level,
                                                     W.statistic = statTest$statistic,
                                                     p.value = statTest$p.value,
                                                     Significance = " "))
    }
  }
}

print("Test done!")
write.csv(t.testDF, file = paste("Tests\\", tissue,"\\T.test_proportions.csv", sep = ""))


# Plots comparing the average proportion of coverage of each gene region by a particular modification for R-genes and controls.

dataToUse <- allResultsAverageProportions[grepl(tissue, allResultsAverageProportions$Tissue) & 
                                            !grepl("luster", allResultsAverageProportions$Tissue),]

for (mod in epiMods) {
  df <- dataToUse[dataToUse$Modification==mod,]
  
  for (level in df$Expression) {
    
    df1 <- df[df$Expression==level,]
    
    controlSampleSize <- sum(unique(df1[grepl("control", df1$Tissue), "SampleSize"]))
    RgeneSampleSize <- sum(unique(df1[grepl("NLRs", df1$Tissue), "SampleSize"]))
    
    if (RgeneSampleSize > 10) {
      plot <- ggplot(df1, aes(x = axisGroup, y = Proportion, color = Tissue)) + 
        scale_x_continuous(limits = c(-60, 140), breaks = seq(-60, 140, 20), labels = axisText) +
        scale_y_continuous(limits = c(0,1), expand = c(0,0)) + 
        scale_colour_manual("Gene set", limits = c(paste("control1_", tissue, sep = ""), paste("NLRs_", tissue, sep = "")), 
                            values=c("grey43", "black"), labels = c(paste("Controls (n = ", controlSampleSize, ")", sep = ""), 
                                                                    paste("R-genes (n = ", RgeneSampleSize, ")", sep = ""))) +
        geom_line(aes(group = Tissue),linewidth = 1) +
        geom_point(aes(group = Tissue), size = 1.5) + theme_minimal() +
        labs(x = "", y = "Average proportion of gene region", title = paste(mod, "-", level, sep = " ")) +
        geom_vline(xintercept=0, color="grey", linewidth=1) +
        coord_cartesian(ylim= c(0,1), clip = "off") + theme(plot.margin = unit(c(1,1,.3,1), "lines")) +
        annotation_custom(textGrob("% of gene length from TSS", gp=gpar(fontsize=16, col = "grey33")),xmin=0,xmax=100,ymin=-.17,ymax=-.17) + 
        annotation_custom(textGrob("Gene region", gp=gpar(fontsize=16)),xmin=0,xmax=100,ymin=-.25,ymax=-.25) +
        theme(axis.text.x = element_text(size = 14, colour = "black", angle = 45, vjust = 1, hjust = 1), axis.text.y = element_text(size = 14,colour = "black"), 
              axis.title.y = element_text(size = 16, vjust = 2), plot.title = element_text(hjust = .5, size = 16),
              legend.text = element_text(size = 12), legend.title = element_text(size = 14), axis.line = element_line(linewidth = .6))
      
      ggsave(paste("Graphs\\Enrichment\\", tissue, "\\", mod, " ", level, " Enrichment.pdf", sep = ""), plot = plot, width = 12, height = 6)  
    }
  }
}