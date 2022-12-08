library(ggplot2)
library(data.table)
library(grid)

# Fisher's Exact Test - are R-genes enriched amongst those that possess a particular chromatin modification?

df <- allResultsFrequencies[grepl(tissue, allResultsFrequencies$SampleGenes),]

hypergeometricTest <- data.frame()

for (mod in unique(allResultsFrequencies$Modification)) {
  df1 <- df[df$Modification==mod,]
  
  for (r in unique(allResultsFrequencies$Region)) {
    df2 <- df1[df1$Region==r,]
    
    for (level in unique(allResultsFrequencies$Expression)) {
      df3 <- df2[df2$Expression==level,]
      df3 <- df3[!grepl("luster", df3$SampleGenes),]
      
      # Create a list of genes from each each sample set.
      geneList <- c()
      
      for (row in 1:nrow(df3)) {
        geneList <- append(geneList, rep(df3[row,"SampleGenes"], times = df3[row, "n"]))
      }
      
      meanNLRs <- c()
      for (n in 1:10) {
        geneSample <- sample(geneList, length(geneList)*0.1)
        meanNLRs <- append(meanNLRs, length(geneSample[grepl("NLRs", geneSample)]))
      }
      
      meanNLRs <- mean(meanNLRs)
      
      if (length(geneSample) > 1) {
        
        statTest <- fisher.test(matrix(c(meanNLRs, length(geneList[grepl("NLRs", geneList)])-meanNLRs,
                                         length(geneSample) - meanNLRs, length(geneList[grepl("control", geneList)])-length(geneSample) - meanNLRs), 2,2), 
                                alternative = "less")
        
        hypergeometricTest <- rbind(hypergeometricTest, data.frame(Expression = level,
                                                                   Modification = mod,
                                                                   Region = r,
                                                                   Sample.R.genes = meanNLRs,
                                                                   Sample.all.genes = length(geneSample),
                                                                   All.R.genes = length(geneList[grepl("NLRs", geneList)]),
                                                                   All.genes = length(geneList),
                                                                   p.value = statTest$p.value))
      }
    }
  }
  print("Test done")
}

write.csv(hypergeometricTest, file = paste("Tests\\", tissue, "_Fisher.Test_frequencies.csv", sep = ""))

# Plots comparing the occurrence of chromatin modifications in the seedlings of R-genes and controls.
axisText <- c("Intergenic", "Promotor \n(1kb)", "Promotor \n(500bp)", "TSS",
              "20%", "40%", "60%", "80%", "100%", 
              "Downstream \n(200bp)", "Intergenic")

dataToUse <- allResultsFrequencies[grepl(tissue, allResultsFrequencies$SampleGenes) & 
                                     !grepl("luster", allResultsFrequencies$SampleGenes),]
for (mod in epiMods) {
  df <- dataToUse[dataToUse$Modification==mod,]
  
  for (level in df$Expression) {
    
    df1 <- df[df$Expression==level,]
    
    controlSampleSize <- sum(unique(df1[grepl("control", df1$SampleGenes), "n"]))
    RgeneSampleSize <- sum(unique(df1[grepl("NLRs", df1$SampleGenes), "n"]))
    
    if (RgeneSampleSize > 10) {
      plot <- ggplot(df1, aes(x = axisGroup, y = Frequency, color = SampleGenes)) + 
        scale_x_continuous(limits = c(-60, 140), breaks = seq(-60, 140, 20), labels = axisText) +
        scale_y_continuous(limits = c(0,100), expand = c(0,0)) + 
        geom_line(aes(group = SampleGenes),linewidth = 1) +
        geom_point(aes(group = SampleGenes), size = 1.5) + theme_minimal() + 
        scale_colour_manual("Gene set", limits = c(paste("control1_", tissue, sep = ""), paste("NLRs_", tissue, sep = "")), 
                            values=c("grey43", "black"), labels = c(paste("Controls (n = ", controlSampleSize, ")", sep = ""), 
                                                                    paste("R-genes (n = ", RgeneSampleSize, ")", sep = ""))) +
        labs(x = "", y = "% Genes modified", title = paste(mod, "-", level, sep = " ")) +
        geom_vline(xintercept=0, color="grey", size=1) +
        coord_cartesian(ylim= c(0,100), clip = "off") + theme(plot.margin = unit(c(1,1,.5,1), "lines")) +
        annotation_custom(textGrob("% of gene length from TSS", gp=gpar(fontsize=16, col = "grey33")),xmin=0,xmax=100,ymin=-22,ymax=-22) + 
        annotation_custom(textGrob("Gene region", gp=gpar(fontsize=16)),xmin=0,xmax=100,ymin=-30,ymax=-30) +
        theme(axis.text.x = element_text(size = 14, colour = "black", angle = 45, vjust = 1, hjust = 1), axis.text.y = element_text(size = 14,colour = "black"), 
              axis.title.y = element_text(size = 16, vjust = 2), plot.title = element_text(hjust = .5, size = 16),
              legend.text = element_text(size = 12), legend.title = element_text(size = 14), axis.line = element_line(linewidth = .6)) 
      
      ggsave(paste("Graphs\\Enrichment\\", tissue, "\\Percentage genes associated with ", mod, "_", level, ".pdf", sep = ""), plot = plot, width = 12, height = 6)
    }
  }
}
  