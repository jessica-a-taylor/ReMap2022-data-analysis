library(ggplot2)
library(data.table)
library(grid)
library(readr)


# Fisher's Exact Test - are R-genes enriched amongst those that possess a particular chromatin modification?
geneCount <- as.data.frame(read_csv("Data\\Seedling\\Gene count.txt"))
geneCount <- geneCount[,-1]

df <- allResultsFrequencies[grepl(tissue, allResultsFrequencies$SampleGenes),]

hypergeometricTest <- data.frame()

for (mod in unique(allResultsFrequencies$Modification)) {
  df1 <- df[df$Modification==mod,]
  
  for (r in unique(allResultsFrequencies$Region)) {
    df2 <- df1[df1$Region==r,]
    
    for (level in c("No Expression", "Low Expression")) {
      df3 <- df2[df2$Expression==level,]
      df3 <- df3[!grepl("luster", df3$SampleGenes),]
      
      if (nrow(df3) > 1) {
        
        # Calculate the number of modified and unmodified R-genes and control genes.
        
        RgenesModified <- geneCount[grepl("NLRs", geneCount$GeneSet) & !grepl("luster", geneCount$GeneSet) &
                                      grepl(level, geneCount$GeneSet), "GeneCount"] * df3[grepl("NLRs", df3$SampleGenes), "Frequency"]/100
        
        RgenesUnmodified <- geneCount[grepl("NLRs", geneCount$GeneSet) & !grepl("luster", geneCount$GeneSet) &
                                        grepl(level, geneCount$GeneSet), "GeneCount"] - RgenesModified
        
        allControls <- c()
        controlGenesModified <- c()
        
        for (row in 1:nrow(df3[grepl("control", df3$SampleGenes),])) {
          allControls <- append(allControls, geneCount[grepl(df3[row, "SampleGenes"], geneCount$GeneSet) &
                                                         grepl(level, geneCount$GeneSet), "GeneCount"])
          
          controlGenesModified <- append(controlGenesModified,
                                         geneCount[grepl(df3[row,"SampleGenes"], geneCount$GeneSet) & 
                                                     grepl(level, geneCount$GeneSet), "GeneCount"] * df3[row, "Frequency"]/100)
        }
        allControls <- sum(allControls)
        controlGenesModified <- sum(controlGenesModified)
        
        statTest <- fisher.test(data.frame(NLR = c(RgenesModified, RgenesUnmodified),
                                           Control = c(controlGenesModified, allControls-controlGenesModified),
                                           row.names = c("Modified", "Unmodified")), alternative = "less")
        
        if (statTest$p.value <= 0.05 & statTest$p.value > 0.01) {
          hypergeometricTest <- rbind(hypergeometricTest, data.frame(Expression = level,
                                                                     Modification = mod,
                                                                     Region = r,
                                                                     Modified.Rgenes = RgenesModified,
                                                                     Unmodified.Rgenes = RgenesUnmodified,
                                                                     Modified.controls = controlGenesModified,
                                                                     Unmodified.controls = allControls-controlGenesModified,
                                                                     p.value = statTest$p.value,
                                                                     Significance = "*"))
        } else if (statTest$p.value <= 0.01 & statTest$p.value > 0.001) {
          hypergeometricTest <- rbind(hypergeometricTest, data.frame(Expression = level,
                                                                     Modification = mod,
                                                                     Region = r,
                                                                     Modified.Rgenes = RgenesModified,
                                                                     Unmodified.Rgenes = RgenesUnmodified,
                                                                     Modified.controls = controlGenesModified,
                                                                     Unmodified.controls = allControls-controlGenesModified,
                                                                     p.value = statTest$p.value,
                                                                     Significance = "**"))
        } else if (statTest$p.value <= 0.001){
          hypergeometricTest <- rbind(hypergeometricTest, data.frame(Expression = level,
                                                                     Modification = mod,
                                                                     Region = r,
                                                                     Modified.Rgenes = RgenesModified,
                                                                     Unmodified.Rgenes = RgenesUnmodified,
                                                                     Modified.controls = controlGenesModified,
                                                                     Unmodified.controls = allControls-controlGenesModified,
                                                                     p.value = statTest$p.value,
                                                                     Significance = "***"))
        }
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
  
  for (level in c("No Expression", "Low Expression")) {
    
    df1 <- df[df$Expression==level,]
    
    RgeneSampleSize <- geneCount[grepl("NLR", geneCount$GeneSet) & grepl(level, geneCount$GeneSet) &
                                     !grepl("luster", geneCount$GeneSet), "GeneCount"]
    controlSampleSize <- sum(geneCount[grepl("control", geneCount$GeneSet) & grepl(level, geneCount$GeneSet), "GeneCount"])
    
  
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
  
