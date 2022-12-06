library(ggplot2)
library(data.table)
library(grid)
library(hash)

# Wilcox & Kolmogorov-Smirnov Tests - is there a significant difference in the average proportion of coverage of  
# each gene region by a particular modification between R-genes and controls?

# Compare R-genes to each control gene set individually.
controlSets <- c("control1","control2","control3","control4","control5",
                 "control6","control7","control8","control9","control10")

wilcoxHash <- hash()
ksHash <- hash()

dataToUse <- allResultsProportions[grepl(tissue, allResultsProportions$SampleGenes) & 
                              !grepl("luster", allResultsProportions$SampleGenes),]

for (set in controlSets) {
  wilcoxDF <- data.frame()
  ksDF <- data.frame()
  
  
  for (mod in unique(allResultsProportions$Modification)) {
    df <- dataToUse[dataToUse$Modification==mod,]
    
    for (r in unique(allResultsProportions$Region)) {
      df1 <- df[df$Region==r,]
      
      for (level in df$Expression) {
        df2 <- df1[df1$Expression == level,]
        
        wilcoxTest <- wilcox.test(Proportion~SampleGenes, df2[c(which(df2$SampleGenes == paste(set, "_", tissue, sep = "")), which(df2$SampleGenes == paste("NLRs_", tissue))),])
        
        if (is.na(wilcoxTest$p.value)==TRUE) {
          wilcoxDF <- wilcoxDF
        } else if (wilcoxTest$p.value <= 0.05) {
          wilcoxDF <- rbind(wilcoxDF, data.frame(Modification = mod,
                                                 Region = r,
                                                 W.statistic = wilcoxTest$statistic,
                                                 p.value = wilcoxTest$p.value,
                                                 Significance = "*"))
        } else  wilcoxDF <- rbind(wilcoxDF, data.frame(Modification = mod,
                                                       Region = r,
                                                       W.statistic = wilcoxTest$statistic,
                                                       p.value = wilcoxTest$p.value,
                                                       Significance = " "))
        
        ksTest <- ks.test(df2[c(which(df2$SampleGenes == paste(set, "_", tissue, sep = ""))),"Proportion"], df2[c(which(df2$SampleGenes == paste("NLRs_", tissue))),"Proportion"])
        
        if (is.na(ksTest$p.value)==TRUE) {
          ksDF <- ksDF
        } else if (ksTest$p.value <= 0.05) {
          ksDF <- rbind(ksDF, data.frame(Modification = mod,
                                         Region = r,
                                         W.statistic = ksTest$statistic,
                                         p.value = ksTest$p.value,
                                         Significance = "*"))
        } else ksDF <- rbind(ksDF, data.frame(Modification = mod,
                                              Region = r,
                                              W.statistic = ksTest$statistic,
                                              p.value = ksTest$p.value,
                                              Significance = " "))
      }
    }
  }
  wilcoxHash[[set]] <- wilcoxDF
  ksHash[[set]] <- ksDF
}

wilcoxDF <- wilcoxHash[[set]][,c(1:2)]
ksDF <- ksHash[[set]][,c(1:2)]


for (set in names(wilcoxHash)) {
  wilcoxDF <- cbind(wilcoxDF, wilcoxHash[[set]][,c(3:5)])
  ksDF <- cbind(ksDF, ksHash[[set]][,c(3:5)])
}

write.csv(wilcoxDF, file = paste("Tests\\", tissue,"\\Wilcox.test_proportions_", level, ".csv", sep = ""))
write.csv(ksDF, file = paste("Tests\\", tissue,"\\Kolmogorov-Smirnov.test_proportions_", level, ".csv", sep = ""))



# Compare R-genes to all control genes.

dataToUse <- allResultsProportions[grepl(tissue, allResultsProportions$SampleGenes) & 
                                     !grepl("luster", allResultsProportions$SampleGenes),]

for (level in df$Expression) {
  df <- dataToUse[dataToUse$Expression==level,]
  
  wilcoxDF <- data.frame()
  ksDF <- data.frame()
  
  for (mod in unique(allResultsProportions$Modification)) {
    df1 <- df[df$Modification==mod,]
    
    for (r in unique(allResultsProportions$Region)) {
      df2 <- df1[df1$Region == r,]
      
      sampleGenesList <- c()
      for (row in 1:nrow(df2)) {
        if (grepl("control", df2[row, "SampleGenes"])) {
          sampleGenesList <- append(sampleGenesList, "control")
        } else sampleGenesList <- append(sampleGenesList, "NLR")
      }
      
      df2 <- cbind(df2, sampleGenesList)
      
      wilcoxTest <- wilcox.test(Proportion~sampleGenesList, df2)
      
      if (is.na(wilcoxTest$p.value)==TRUE) {
        wilcoxDF <- wilcoxDF
      } else if (wilcoxTest$p.value <= 0.05) {
        wilcoxDF <- rbind(wilcoxDF, data.frame(Modification = mod,
                                               Region = r,
                                               W.statistic = wilcoxTest$statistic,
                                               p.value = wilcoxTest$p.value,
                                               Significance = "*"))
      } else  wilcoxDF <- rbind(wilcoxDF, data.frame(Modification = mod,
                                                     Region = r,
                                                     W.statistic = wilcoxTest$statistic,
                                                     p.value = wilcoxTest$p.value,
                                                     Significance = " "))
      
      ksTest <- ks.test(df2[grepl("control", df2$sampleGenesList), "Proportion"], df2[grepl("NLR", df2$sampleGenesList), "Proportion"])
      
      if (is.na(ksTest$p.value)==TRUE) {
        ksDF <- ksDF
      } else if (ksTest$p.value <= 0.05) {
        ksDF <- rbind(ksDF, data.frame(Modification = mod,
                                       Region = r,
                                       W.statistic = ksTest$statistic,
                                       p.value = ksTest$p.value,
                                       Significance = "*"))
      } else ksDF <- rbind(ksDF, data.frame(Modification = mod,
                                            Region = r,
                                            W.statistic = ksTest$statistic,
                                            p.value = ksTest$p.value,
                                            Significance = " "))
    }
  }
  write.csv(wilcoxDF, file = paste("Tests\\", tissue,"\\Wilcox.test_proportions_", level, ".csv", sep = ""))
  write.csv(ksDF, file = paste("Tests\\", tissue,"\\Kolmogorov-Smirnov.test_proportions_", level, ".csv", sep = ""))
}


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