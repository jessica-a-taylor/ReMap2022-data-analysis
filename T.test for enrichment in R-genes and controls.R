library(ggplot2)
library(data.table)
library(grid)
library(readr)

# T-Test - is there a significant difference in the average proportion of coverage of  
# each gene region by a particular modification between R-genes and controls?
# Plots comparing the average proportion of coverage of each gene region by a particular modification in R-genes and controls.

# Compare R-genes to each control gene set.
dataToUse <- resultsProportions[[analysis]][grepl(tissue, resultsProportions[[analysis]]$dataToAnalyse) & 
                              !grepl("luster", resultsProportions[[analysis]]$dataToAnalyse),]

t.testDF <- data.frame()


for (mod in unique(resultsProportions[[analysis]]$Modification)) {
  df <- dataToUse[dataToUse$Modification==mod,]
  
  for (r in unique(resultsProportions[[analysis]]$Region)) {
    df1 <- df[df$Region==r,]
    
    for (level in c("No Expression", "Low Expression")) {
      df2 <- df1[df1$Expression == level,]
      
        statTest <- t.test(df2[c(which(df2$dataToAnalyse == paste("NLRs_", tissue, sep = ""))), "Proportion"], mu = mean(df2[grepl("control", df2$dataToAnalyse) & 
                                                                                                                               grepl(tissue, df2$dataToAnalyse), "Proportion"]))
      
      if (is.na(statTest$p.value)==TRUE) {
        t.testDF <- t.testDF
      } else if (0.05 >= statTest$p.value & statTest$p.value > 0.01) {
        t.testDF <- rbind(t.testDF, data.frame(Modification = mod,
                                               Region = r,
                                               Expression = level,
                                               W.statistic = statTest$statistic,
                                               p.value = statTest$p.value,
                                               Significance = "*"))
        
      } else if (0.01 >= statTest$p.value & statTest$p.value > 0.001) {
        t.testDF <- rbind(t.testDF, data.frame(Modification = mod,
                                               Region = r,
                                               Expression = level,
                                               W.statistic = statTest$statistic,
                                               p.value = statTest$p.value,
                                               Significance = "**"))
        
      } else if (statTest$p.value <= 0.001) {
        t.testDF <- rbind(t.testDF, data.frame(Modification = mod,
                                               Region = r,
                                               Expression = level,
                                               W.statistic = statTest$statistic,
                                               p.value = statTest$p.value,
                                               Significance = "***"))
        
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
write.csv(t.testDF, file = paste("Tests\\", analysis, "\\", tissue,"\\T.test_proportions.csv", sep = ""))


# Plots comparing the average proportion of coverage of each gene region by a particular modification for R-genes and controls.
geneCount <- as.data.frame(read_csv(paste("Data\\", analysis, "\\", tissue, "\\Gene count.txt", sep = "")))
geneCount <- geneCount[,-1]

dataToUse <- resultsAverageProportions[[analysis]][grepl(tissue, resultsAverageProportions[[analysis]]$Tissue) & 
                                            !grepl("luster", resultsAverageProportions[[analysis]]$Tissue),]

for (mod in epiMods) {
  df <- dataToUse[dataToUse$Modification==mod,]
  
  for (level in unique(df$Expression)) {
    
    df1 <- df[df$Expression==level,]
    
    RgeneSampleSize <- geneCount[grepl("NLR", geneCount$GeneSet) & grepl(level, geneCount$GeneSet) &
                                   !grepl("luster", geneCount$GeneSet), "GeneCount"]
    controlSampleSize <- sum(geneCount[grepl("control", geneCount$GeneSet) & grepl(level, geneCount$GeneSet), "GeneCount"])
    
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
      
      ggsave(paste("Graphs\\Enrichment\\", analysis, "\\", tissue, "\\", mod, " ", level, " Enrichment.pdf", sep = ""), plot = plot, width = 12, height = 6)  
    }
  }
}