# Determine whether the enrichment of each chromatin mark is significantly more similar between genes within a cluster 
# than between all genes.

if (analysis == "PlantExp data") {
  dataToAnalyse <- sampleGenesPlantExp
} else if (analysis == "RNA-seq data") {
  dataToAnalyse <- sampleGenesRNAseq
}

# Choose ecotype and tissue for dataToAnalyse.
# Options: leafGenes, rootGenes, seedlingGenes
if (tissue == "leaves") {
  tissueForAnalysis <- "leafGenes"
} else if (tissue == "root") {
  tissueForAnalysis <- "rootGenes"
} else if (tissue == "seedlings") {
  tissueForAnalysis <- "seedlingGenes"
}

# Create a hash for storing the proportion of each R-gene region modified with each mark.
dataToAnalyseProportions_ClusterAnalysis <- hash()

for (test in c("clusteredNLRs", "notClusteredNLRs")) {
  
  for (cluster in unique(dataToAnalyse[[test]]$Clustering)) {
    dataToUse <- dataToAnalyse[[test]][dataToAnalyse[[test]]$Clustering == cluster,]
    
    # Create a hash with the ReMap data in a particular tissue for the current set of genes. 
    allModifications <- ReMapPerGene(dataToUse, tissueForAnalysis)
    
    # For each gene in the current set of genes, create a new hash with the occurrences of each chromatin modification.
    geneModifications <- modificationOccurrences(allModifications)
    
    rm(allModifications)
    
    # For each gene in the current set of genes, merge the overlapping occurrences of each modification.
    allOverlaps <- mergeOverlappingModifications(geneModifications)
    
    
    # Determine the % R-genes with a chromatin mark in each gene region (frequency)
    # and the proportion of each gene region with that mark.
    geneRegions <- getGeneCoordinates(dataToUse)
    
    modProportionPerRegion <- modProportionsFunction(geneRegions, allOverlaps, epiMods)
    
    # Add a column to modProportionPerRegion with the numbers for 
    # each gene region that will correspond with their position on the x axis.
    modProportionPerRegion <- geneRegionAxisLocations(modProportionPerRegion, geneRegions)
    
    # Add a column to modProportionPerRegion with the current expression level.
    modProportionPerRegion <- expressionColumn(modProportionPerRegion, cluster)
    
    # Store final results on the appropriate hash.
    dataToAnalyseProportions_ClusterAnalysis[[test]][[cluster]] <- modProportionPerRegion
  }
}


# Merge all data from all sample gene sets into one big dataframe.
allResultsProportions_ClusterAnalysis <- data.frame()

for (test in names(dataToAnalyseProportions_ClusterAnalysis)) {
  for (cluster in names(dataToAnalyseProportions_ClusterAnalysis[[test]])) {
    
    df2 <- dataToAnalyseProportions_ClusterAnalysis[[test]][[cluster]]
    df2 <- cbind(df2, data.frame(dataToAnalyse = rep(test, times = nrow(df2))))
    
    allResultsProportions_ClusterAnalysis <- rbind(allResultsProportions_ClusterAnalysis, df2)
  }
}


# Divide the list of chromatin modifications into 3 lists.
modList <- list(miniList1 = unique(allResultsProportions$Modification)[1:7],
                miniList2 = unique(allResultsProportions$Modification)[8:14],
                miniList3 = unique(allResultsProportions$Modification)[15:21])

# Calculate the difference in the proportion of each R-gene region covered by each modification between R-genes of the same cluster.
for (modMiniList in names(modList)) {
  jobRunScript("Pairwise comparisons of enrichment.R", name = modMiniList, importEnv = TRUE)
}

modificationDifference <- data.frame()

for (comparison in c("Cluster modification comparisons", "All R-gene modification comparisons")) {
  for (file in list.files(path = "Data\\",pattern = comparison)) {
    df <- as.data.frame(read_xlsx(paste("Data\\", file, sep = "")))
    df <- df[,c(1:4)]
    df <- cbind(df, data.frame(Comparison = comparison))
    
    modificationDifference <- rbind(modificationDifference, df)
  }
}

# T-Test - is there a significant difference in the degree of variation among R-genes within a cluster
# compared to all R-genes?

statTestDF <- data.frame()

for (mod in unique(allResultsProportions$Modification)) {
  df <- modificationDifference[modificationDifference$Modification==mod,]
  
  compareAverages <- data.frame()
  
  for (r in unique(df$Region)) {
    df1 <- df[df$Region==r,]
    statTest <- t.test(df1[df1$Comparison=="Cluster modification comparisons","ModificationDifference"],
                       mu = mean(df1[df1$Comparison=="All R-gene modification comparisons","ModificationDifference"]))
    
    compareAverages <- rbind(compareAverages, data.frame(Region = r,
                                                         axisGroup = df1$axisGroup[1],
                                                         ModificationDifference = mean(df1[df1$Comparison=="Cluster modification comparisons","ModificationDifference"]),
                                                         Comparison = "Cluster modification comparisons"))
    
    compareAverages <- rbind(compareAverages, data.frame(Region = r,
                                                         axisGroup = df1$axisGroup[1],
                                                         ModificationDifference = mean(df1[df1$Comparison=="All R-gene modification comparisons","ModificationDifference"]),
                                                         Comparison = "All R-gene modification comparisons"))
    
    if (is.na(statTest$p.value)==TRUE) {
      statTestDF <- statTestDF
    } else if (0.05 >= statTest$p.value & statTest$p.value > 0.01) {
      statTestDF <- rbind(statTestDF, data.frame(Modification = mod,
                                                 Region = r,
                                                 Statistic = statTest$statistic,
                                                 p.value = statTest$p.value,
                                                 Significance = "*"))
      
    } else if (0.01 >= statTest$p.value & statTest$p.value > 0.001) {
      statTestDF <- rbind(statTestDF, data.frame(Modification = mod,
                                                 Region = r,
                                                 Statistic = statTest$statistic,
                                                 p.value = statTest$p.value,
                                                 Significance = "**"))
      
    } else if (statTest$p.value <= 0.001) {
      statTestDF <- rbind(statTestDF, data.frame(Modification = mod,
                                                 Region = r,
                                                 Statistic = statTest$statistic,
                                                 p.value = statTest$p.value,
                                                 Significance = "***"))
      
    } else  statTestDF <- rbind(statTestDF, data.frame(Modification = mod,
                                                       Region = r,
                                                       Statistic = statTest$statistic,
                                                       p.value = statTest$p.value,
                                                       Significance = " "))
  }
  
  # Plots comparing the degree of variation among R-genes within a cluster compared to all R-genes.
  facetLabels <- c(paste(mod, " - Between R-genes in the same cluster", sep = ""), paste(mod, " - Between all R-genes", sep = ""))
  names(facetLabels) <- unique(df$Comparison)
  
  boxplot <- ggplot(df, aes(x = axisGroup, y = ModificationDifference)) + 
    scale_x_continuous(limits = c(-70, 150), breaks = seq(-60, 140, 20), labels = axisText) +
    scale_y_continuous(limits = c(0,1), expand = c(0,0)) + 
    geom_boxplot(aes(group = Region)) + theme_minimal() + coord_cartesian(ylim= c(0,1), clip = "off") +
    labs(x = "", y = "Difference in proportion of gene region modified") +
    geom_vline(xintercept=0, color="grey", linewidth=1) + theme(plot.margin = unit(c(1,1,.3,1), "lines")) +
    annotation_custom(textGrob("% of gene length from TSS", gp=gpar(fontsize=14, col = "grey33")),xmin=-10,xmax=100,ymin=-.15,ymax=-.15) + 
    annotation_custom(textGrob("Gene region", gp=gpar(fontsize=14)),xmin=-10,xmax=100,ymin=-.22,ymax=-.22) +
    theme(axis.text.x = element_text(size = 12, colour = "black", angle = 45, vjust = 1, hjust = 1), axis.text.y = element_text(size = 12,colour = "black"), 
          axis.title.y = element_text(size = 14, vjust = 2), strip.text = element_text(size = 16, vjust = 4),
          legend.text = element_text(size = 12), legend.title = element_text(size = 14), axis.line = element_line(linewidth = .6)) + 
    facet_wrap(~Comparison, labeller = labeller(Comparison = facetLabels)) + 
    stat_summary(fun="mean", geom="point", color="black", size=2) +
    stat_summary(fun="mean", geom="line", color="black", linewidth=1)
  
  ggsave(paste("Graphs\\Clusters vs. Background\\", analysis, "\\", tissue, "\\", mod, "_boxplot.pdf", sep = ""), 
         plot = boxplot, width = 16, height = 6)
  
  
  plot <- ggplot(compareAverages, aes(axisGroup, y = ModificationDifference, color = Comparison)) +
    scale_x_continuous(limits = c(-70, 150), breaks = seq(-60, 140, 20), labels = axisText) +
    scale_y_continuous(limits = c(0,1), expand = c(0,0)) + 
    geom_point(size=2) + geom_line(linewidth = 1) + theme_minimal() + coord_cartesian(ylim= c(0,1), clip = "off") +
    geom_vline(xintercept=0, color="grey", linewidth=1) + theme(plot.margin = unit(c(1,1,.3,1), "lines")) +
    annotation_custom(textGrob("% of gene length from TSS", gp=gpar(fontsize=14, col = "grey33")),xmin=-10,xmax=100,ymin=-.15,ymax=-.15) + 
    annotation_custom(textGrob("Gene region", gp=gpar(fontsize=14)),xmin=-10,xmax=100,ymin=-.22,ymax=-.22) +
    theme(axis.text.x = element_text(size = 12, colour = "black", angle = 45, vjust = 1, hjust = 1), axis.text.y = element_text(size = 12,colour = "black"), 
          axis.title.y = element_text(size = 14, vjust = 2), strip.text = element_text(size = 16, vjust = 4),
          legend.text = element_text(size = 12), legend.title = element_text(color = "white"), axis.line = element_line(linewidth = .6)) +
    scale_color_manual(values = c("coral2", "darkslategrey"), labels = c("Between all R-genes", "Between R-genes in \nthe same cluster")) +
    labs(x = "", y = "Average difference in enrichment")
  
  ggsave(paste("Tests\\Clusters vs. Background\\", analysis, "\\", tissue, "\\", mod, "_lineplot.pdf", sep = ""), 
         plot = plot, width = 10, height = 6)
}

write.csv(statTestDF, file = paste("Tests\\Clusters vs. Background\\", analysis, "\\", tissue, "\\", "T_test.csv", sep = ""))
