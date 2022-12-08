# Calculate the difference in the proportion of each R-gene region covered by each modification between R-genes of the same cluster.

modificationDifference_Clusters <- data.frame()

for (mod in modList[[modMiniList]]) {
  df <- allResultsProportions_ClusterAnalysis[allResultsProportions_ClusterAnalysis$Modification == mod,]
  
  for (r in unique(df$Region)) {
    df1 <- df[df$Region == r,]
    
    for (cluster in unique(df1[grepl("cluster", df1$Expression),"Expression"])) {
      df2 <- df1[df1$Expression == cluster,]
      
      # For all pairwise comparisons,
      for (i in 1:nrow(df2)) {
        for (j in 1:nrow(df2)) {
          if (i != j & j > i) {
            # calculate the % difference in the proportion of each gene region covered by each modification,
            # using the abs() function to make all values positive.
           modificationDifference_Clusters <- rbind(modificationDifference_Clusters, 
                                                            data.frame(Region = r,
                                                                       Modification = mod,
                                                                       axisGroup = df2$axisGroup[1],
                                                                       ModificationDifference = abs(df2[i, "Proportion"]-df2[j, "Proportion"]),
                                                                       Cluster = cluster))
          }
          else next
        }
      }
    }
  }
  print(mod)
}

write.csv(modificationDifference_Clusters, file = paste(modMiniList, " - Cluster modification comparisons.csv", sep = ""))

# Calculate the difference in the proportion of each gene region covered by each modification between all R-genes.
modificationDifference_all <- data.frame()


for (mod in modList[[modMiniList]]) {
  df <- allResultsProportions_ClusterAnalysis[allResultsProportions_ClusterAnalysis$Modification == mod,]
  
  for (r in unique(df$Region)) {
    df1 <- df[df$Region == r,]
    
    for (i in 1:nrow(df1)) {
      for (j in 1:nrow(df1)) {
        if (i != j & j > i) {
          modificationDifference_all <- rbind(modificationDifference_all, 
                                              data.frame(Region = r,
                                                         Modification = mod,
                                                         axisGroup = df1$axisGroup[1],
                                                         ModificationDifference = abs(df1[i, "Proportion"]-df1[j, "Proportion"])))
        }
      }
    }
    print(r)
  }
  print(mod)
}

write.csv(modificationDifference_all, file = paste(modMiniList, " - All R-gene modification comparisons.csv", sep = ""))

# Create a new new dataframe merging `modificationDifference_Clusters` and `modificationDifference_all`.
modificationDifference <- data.frame(Comparison = rep("Between clustered R-genes", times = length(modificationDifference_Clusters)),
                                     ModificationDifference = modificationDifference_Clusters[,c(1:4)])

modificationDifference <- rbind(modificationDifference, data.frame(Comparison = rep("Between all R-genes", times = length(modificationDifference_all)),
                                                                   ModificationDifference = modificationDifference_all))

colnames(modificationDifference) <- c("Comparison", "Region", "Modification", "axisGroup", "Difference")


# Wilcox Test - is there a significant difference in the degree of variation among R-genes within a cluster
# compared to all R-genes?

statTestDF <- data.frame(Modification = character(),
                         Region = character(),
                         W.statistic = numeric(),
                         p.value = numeric())

for (mod in modList[[modMiniList]]) {
  df <- modificationDifference[modificationDifference$Modification==mod,]
  
  for (r in unique(df$Region)) {
    df1 <- df[df$Region==r,]
    statTest <- wilcox.test(Difference~Comparison, df1)
    
    statTestDF <- rbind(statTestDF, data.frame(Modification= mod,
                                               Region= r,
                                               W.statistic = statTest$statistic,
                                               p.value = statTest$p.value))
  }
  
  # Plots comparing the degree of variation among R-genes within a cluster compared to all R-genes.
  
  plot <- ggplot(df, aes(x = axisGroup, y = Difference)) + 
    scale_x_continuous(limits = c(-70, 150), breaks = seq(-60, 140, 20), labels = axisText) +
    geom_boxplot(aes(group = Region)) + theme_minimal() + coord_cartesian(ylim= c(0,1), clip = "off") +
    labs(x = "", y = "Difference in proportion of gene region modified") +
    geom_vline(xintercept=0, color="grey", linewidth=1) + theme(plot.margin = unit(c(1,1,1,1), "lines")) +
    annotation_custom(textGrob("% of gene length from TSS", gp=gpar(fontsize=12, col = "grey33")),xmin=-10,xmax=100,ymin=-.23,ymax=-.23) + 
    annotation_custom(textGrob("Gene region", gp=gpar(fontsize=14)),xmin=-10,xmax=100,ymin=-.3,ymax=-.3) +
    theme(axis.text.x = element_text(size = 11, colour = "black", angle = 45, vjust = 1, hjust = 1), axis.text.y = element_text(size = 12,colour = "black"), 
          axis.title.y = element_text(size = 14, vjust = 2), strip.text = element_text(size = 16)) + 
    facet_wrap(~Comparison) + 
    stat_summary(fun="mean", geom="point", color="black", size=2) +
    stat_summary(fun="mean", geom="line", color="black", size=1)
  
  
  ggsave(paste(mod, "_", "Clustered vs Unclustered.pdf", sep = ""), plot = plot, width = 16, height = 6)
}


write.csv(statTestDF, file = paste(modMiniList, " - Difference in proportion of gene region modified.csv", sep = ""))