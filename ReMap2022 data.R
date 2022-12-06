if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggplot2")

library(readxl)
library(karyoploteR)
library(rtracklayer)
library(dplyr)
library(stringr)
library(hash)
library(sets)
library(TxDb.Athaliana.BioMart.plantsmart28)
library(ggplot2)
library(data.table)
library(grid)
library(readr)
library(rstudioapi)

source("Functions\\Get range - merge gene coordinates.R")


# Import all Arabidopsis genes.
Atgenes <- as.data.frame(transcriptsBy(TxDb.Athaliana.BioMart.plantsmart28, by="gene"))
colnames(Atgenes)[2] <- "Gene"

# Remove duplicate genes (different versions).
Atgenes <- Atgenes[-c(which(Atgenes$tx_name == str_match(Atgenes$tx_name, "^([0-9a-zA-Z]+)([.])([2-9])$")[,1])),]


# Remove genes within the centromeric and pericentromeric geneRegions.
pericentromericgeneRegions <- data.frame(Chromosome = c(1:5),
                                         Start = c("11500000", "1100000", "10300000", "1500000", "9000000"),
                                         End = c("17700000", "7200000", "17300000", "6300000", "16000000"))

euchromaticgeneRegions <- data.frame()

for (row in 1:nrow(pericentromericgeneRegions)) {
  df <- Atgenes[c(which(Atgenes$seqnames==row & Atgenes$start < as.numeric(pericentromericgeneRegions[row, "Start"]))),]
  df <- rbind(df, Atgenes[c(which(Atgenes$seqnames==row & Atgenes$end > as.numeric(pericentromericgeneRegions[row, "End"]))),])
  
  euchromaticgeneRegions <- rbind(euchromaticgeneRegions, df)
}

rm(pericentromericgeneRegions, df)

euchromaticgeneRegions <- euchromaticgeneRegions[,-c(1,8,9)]
euchromaticgeneRegions$ranges <- paste(euchromaticgeneRegions$start,"-",euchromaticgeneRegions$end, sep = "")

# Remove duplicate genes.
newEuchromaticgeneRegions <- euchromaticgeneRegions
euchromaticgeneRegions <- data.frame()

for (gene in unique(newEuchromaticgeneRegions$Gene)) {
  euchromaticgeneRegions <- rbind(euchromaticgeneRegions, newEuchromaticgeneRegions[newEuchromaticgeneRegions$Gene==gene,][1,])
}

rm(newEuchromaticgeneRegions)

# Remove TEs from the euchromaticgeneRegions dataframe.
transposableElements <- as.data.frame(read_xlsx("Data\\Arabidopsis TE genes.xlsx"))

withoutTEs <- euchromaticgeneRegions[-c(which(euchromaticgeneRegions$Gene %in% transposableElements$Locus)),]

rm(transposableElements)


# Get 10 sets of random genes and store in a hash from gene dataset of interest.
source("Functions\\Sample random genes.R")

dataToUse <- withoutTEs

sampleGenes <- geneSets(dataToUse)


# Import list of R-genes.
ArabidopsisNLRs <- as.data.frame(read_xlsx("Data\\Arabidopsis NLRs.xlsx", sheet = 1))
clusteredNLRs <- ArabidopsisNLRs[grepl("cluster", ArabidopsisNLRs$Clustering),]
notClusteredNLRs <- ArabidopsisNLRs[c(which(ArabidopsisNLRs$Clustering =="single")),]

NLRgenes <- dataToUse[which(dataToUse$Gene %in% ArabidopsisNLRs$Gene),]
NLRgenes <- cbind(NLRgenes, 
                  data.frame(Clustering = ArabidopsisNLRs[which(ArabidopsisNLRs$Gene %in% NLRgenes$Gene),"Clustering"]))


clusteredNLRgenes <- dataToUse[which(dataToUse$Gene %in% clusteredNLRs$Gene),]
clusteredNLRgenes <- cbind(clusteredNLRgenes, 
                           data.frame(Clustering = clusteredNLRs[which(clusteredNLRs$Gene %in% clusteredNLRgenes$Gene),"Clustering"]))

notClusteredNLRgenes <- dataToUse[which(dataToUse$Gene %in% notClusteredNLRs$Gene),]
notClusteredNLRgenes <- cbind(notClusteredNLRgenes, 
                              data.frame(Clustering = notClusteredNLRs[which(notClusteredNLRs$Gene %in% notClusteredNLRgenes$Gene),"Clustering"]))


# Add R-genes to sampleGenes.
sampleGenes[["NLRs"]] <- NLRgenes
sampleGenes[["clusteredNLRs"]] <- clusteredNLRgenes
sampleGenes[["notClusteredNLRs"]] <- notClusteredNLRgenes

rm(ArabidopsisNLRs, NLRgenes, Atgenes)


# Get filtered expression data for each set of sample genes in each tissue. 
# Add dataframes to sampleGenes for gene sets with particular expression levels.
source("Functions\\PlantExp.R")
exLevel <- c("No Expression", "Low Expression", "Intermediate Expression",
             "High Expression", "V.High Expression")

sampleGenes <- PlantExp(sampleGenes, exLevel)


# Use ReMap2022 data to analyse the enrichment of chromatin marks on the R-genes and controls.
source("Functions\\Modifications per gene.R")
source("Functions\\Coordinates per gene region.R")
source("Functions\\Modification frequencies & proportions.R")

# Import filtered ReMap2022 data.
ReMap <- as.data.frame(read_xlsx("Data\\Filtered ReMap data.xlsx"))

# Create list of chromatin modifications.
epiMods <- unique(ReMap$epiMod)

# For each gene set (R-genes and controls), tissue, and expression level, 
# find the percentage of genes with each chromatin modification (frequency)
# and the enrichment of the modification (proportion) in each gene region.

# Run parallel analyses as background jobs.
for (tissue in c("Leaf", "Root", "Seedling")) {
  jobRunScript("ReMap analysis.R", name = tissue, importEnv = TRUE)
}

# Import the results.
allResultsFrequencies <- data.frame()
allResultsProportions <- data.frame()

for (tissue in c("Leaf", "Root", "Seedling")) {
  allResultsFrequencies <- rbind(allResultsFrequencies, as.data.frame(read_xlsx(paste("Data\\", tissue,"_allResultsFrequencies.xlsx", sep = ""))))
  allResultsProportions <- rbind(allResultsProportions, as.data.frame(read_xlsx(paste("Data\\", tissue,"_allResultsProportions.xlsx", sep = ""))))
}


# Fisher's Exact Test - are R-genes enriched amongst those that possess a particular chromatin modification?
# Plots comparing the occurrence of chromatin modifications in the seedlings of R-genes and controls.
for (tissue in c("Leaf", "Root", "Seedling")) {
  jobRunScript("Fisher's test for enrichment.R", name = paste("Enrichment_", tissue, sep = ""), importEnv = TRUE)
}


# Calculate the mean proportion of overlap and add as a new column to the dataframe.
allResultsAverageProportions <- data.frame()

for (test in unique(allResultsProportions$SampleGenes)) {
  df <- allResultsProportions[allResultsProportions$SampleGenes==test,]
  
    for (level in unique(df$Expression)) { 
      df1 <- df[df$Expression==level,]
      
      if (nrow(df1) >= 1) {
        for (mod in unique(df1$Modification)) {
          df2 <- df1[df1$Modification==mod,]
          
          for (r in unique(df2$Region)) {
            df3 <- df2[df2$Region==r,]
            
            if (nrow(df3) >= 10) {
              allResultsAverageProportions <- rbind(allResultsAverageProportions, data.frame(Region = r,
                                                                                             Modification = mod,
                                                                                             Proportion = mean(df3$Proportion),
                                                                                             Tissue = df3$SampleGenes[1],
                                                                                             axisGroup = df3$axisGroup[1],
                                                                                             Expression = level,
                                                                                             SampleSize = nrow(df3)))
            }
            else allResultsAverageProportions <- allResultsAverageProportions
          }
        }
      } else allResultsAverageProportions <- allResultsAverageProportions
    }
}


# Wilcox & Kolmogorov-Smirnov Tests - is there a significant difference in the average proportion of coverage of  
# each gene region by a particular modification between R-genes and controls?

# Plots comparing the average proportion of coverage of each gene region by a particular modification in the 
# seedlings of R-genes and controls.

for (tissue in c("Leaf", "Root", "Seedling")) {
  jobRunScript("Wilcox test for enrichment.R", name = paste("Enrichment_", tissue, sep = ""), importEnv = TRUE)
}



# Plot the expression of all R-genes and controls.
geneExpression <- data.frame()

for (test in names(sampleGenes[["NLRs_Seedling"]])) {
  if (nrow(sampleGenes[["NLRs_Seedling"]][[test]]) >= 1) {
    
  geneExpression <- rbind(geneExpression, data.frame(Gene = sampleGenes[["NLRs_Seedling"]][[test]]$Gene,
                                                   Expression = sampleGenes[["NLRs_Seedling"]][[test]]$FPKM,
                                                   Level = sampleGenes[["NLRs_Seedling"]][[test]]$ExpressionLevel,
                                                   SampleSize = paste(sampleGenes[["NLRs_Seedling"]][[test]]$ExpressionLevel, 
                                                                      paste("(n = ", length(sampleGenes[["NLRs_Seedling"]][[test]]$Gene), ")", sep = ""), sep = " ")))
  }
  else geneExpression <- geneExpression
}

geneExpression <- geneExpression[order(geneExpression$Gene),] 

plot <- ggplot(geneExpression[c(125:155),], aes(x = Gene, y = Expression, fill = SampleSize)) +
  geom_bar(stat = "identity") + theme_minimal() + labs(x = "R-gene", y = "Expression (FPKM)") +
  theme(axis.text.x = element_text(angle = 90, size = 10), axis.text.y = element_text(size = 10), 
        axis.title = element_text(size = 12), panel.background = element_rect(fill = "white", color = "white")) +
  scale_fill_brewer(palette = "Reds", direction=-1, name = "Expression Level", labels = c(unique(geneExpression$SampleSize))) + 
  coord_cartesian(ylim = c(0,65))

ggsave("R-gene Seedling Expression 5.png", plot = plot, width = 15, height = 6)


for (test in names(sampleGenes)[grepl("Seedling", names(sampleGenes)) & grepl("control", names(sampleGenes))]) {
  geneExpression <- data.frame()
  
  for (n in exLevel) {
    if (nrow(sampleGenes[[test]][[n]]) >= 1) {
      geneExpression <- rbind(geneExpression, data.frame(Gene = sampleGenes[[test]][[n]]$Gene,
                                                         Expression = sampleGenes[[test]][[n]]$FPKM,
                                                         Level = sampleGenes[[test]][[n]]$ExpressionLevel,
                                                         SampleSize = paste(sampleGenes[[test]][[n]]$ExpressionLevel, 
                                                                            paste("(n = ", length(sampleGenes[[test]][[n]]$Gene), ")", sep = ""), sep = " ")))
    }
    else geneExpression <- geneExpression
  }
  
  plot <- ggplot(geneExpression, aes(x = Gene, y = Expression, fill = factor(Level, levels = exLevel))) +
    geom_bar(stat = "identity") + theme_minimal() + labs(x = "R-gene", y = "Expression (FPKM)") +
    theme(axis.text.x = element_text(angle = 45, size = 8, hjust = 1)) +
    scale_fill_brewer(palette = "Reds", direction=1, name = "Expression Level", labels = c(unique(geneExpression$SampleSize)))
  
  ggsave(paste(test,"Expression.pdf", sep = " "), plot = plot, width = 36, height = 6)
}
  
# Expression boxplot.
geneExpression <- data.frame()

for (test in names(sampleGenes)[c(4,52)]) {
  
  for (n in names(sampleGenes[[test]])) {
    geneExpression <- rbind(geneExpression, data.frame(GeneSet = rep(test, times = nrow(sampleGenes[[test]][[n]])),
                                                       Expression = sampleGenes[[test]][[n]]$FPKM))
  }
}

geneExpressionMean <- data.frame(GeneSet = "Controls",
                                 Expression = geneExpression[c(grep("control", geneExpression$GeneSet)), "Expression"])

geneExpressionMean <- rbind(geneExpressionMean, data.frame(GeneSet = "R-genes",
                                                            Expression = geneExpression[c(grep("NLR", geneExpression$GeneSet)), "Expression"]))

geneExpressionCut <- geneExpression[c(which(geneExpression$Expression <= 50)),]

plot <- ggplot(geneExpression, aes(x = GeneSet, y = Expression)) +
                 geom_boxplot(aes(group = GeneSet)) + theme_minimal() + labs(x = "Gene set", y = "Expression (FPKM)") +
                 theme(axis.text.x = element_text(angle = 45, size = 8, hjust = 1))


ggsave("Expression boxplot.pdf", plot = plot, width = 12, height = 6)

statTest <- wilcox.test(Expression~GeneSet, geneExpressionMean)

# Determine whether the expression is more similar between R-genes within a cluster than between all R-genes in seedlings.
clusteredExpression <- data.frame()

for (test in names(sampleGenes[["clusteredNLRs_Seedling"]])) {
  clusteredExpression <- rbind(clusteredExpression, data.frame(Gene = sampleGenes[["clusteredNLRs_Seedling"]][[test]]$Gene,
                                                               Cluster = sampleGenes[["clusteredNLRs_Seedling"]][[test]]$Clustering,
                                                               Expression = sampleGenes[["clusteredNLRs_Seedling"]][[test]]$FPKM))
}

expressionDifference_Clusters <- c()

for (cluster in unique(clusteredExpression$Cluster)) {
  df <- clusteredExpression[clusteredExpression$Cluster == cluster,]
  
  for (i in 1:nrow(df)) {
    for (j in 1:nrow(df)) {
      if (i != j & j > i) {
        expressionDifference_Clusters <- append(expressionDifference_Clusters, abs(df[i, "Expression"]-df[j, "Expression"]))
      }
    }
  }
}

allExpression <- data.frame()

for (test in names(sampleGenes[["NLRs_Seedling"]])) {
  allExpression <- rbind(allExpression, data.frame(Gene = sampleGenes[["NLRs_Seedling"]][[test]]$Gene,
                                                   Expression = sampleGenes[["NLRs_Seedling"]][[test]]$FPKM))
}

expressionDifference_all <- c()

for (i in 1:nrow(allExpression)) {
  for (j in 1:nrow(allExpression)) {
    if (i != j & j > i) {
      expressionDifference_all <- append(expressionDifference_all, abs(allExpression[i, "Expression"]-allExpression[j, "Expression"]))
    }
  }
}

expressionDifference <- data.frame(Comparison = rep("Between clustered genes", times = length(expressionDifference_Clusters)),
                                   ExpressionDifference = expressionDifference_Clusters)

expressionDifference <- rbind(expressionDifference, data.frame(Comparison = rep("Between all genes", times = length(expressionDifference_all)),
                                                               ExpressionDifference = expressionDifference_all))

plot <- ggplot(expressionDifference, aes(x = Comparison, y = ExpressionDifference, fill = Comparison)) +
  geom_boxplot(aes(group = Comparison)) + theme_minimal() + labs(x = "Comparison", y = "Expression Difference (FPKM)") +
  theme(axis.text.x = element_text(size = 10), axis.title.x = element_blank(), legend.position = "none")



# Determine whether the enrichment of each chromatin mark is significantly more similar between genes within a cluster 
# than between all genes.

# Choose ecotype and tissue for analysis.
# Options: ColLeaf, ColRoot
tissueForAnalysis <- "ColLeaf"

# Create a hash for storing the proportion of each R-gene region modified with each mark.
sampleGenesProportions_ClusterAnalysis <- hash()

for (test in names(sampleGenes)[c(1,49)]) {
  
  for (cluster in unique(sampleGenes[[test]]$Clustering)) {
    dataToUse <- sampleGenes[[test]][sampleGenes[[test]]$Clustering == cluster,]

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
    sampleGenesProportions_ClusterAnalysis[[test]][[cluster]] <- modProportionPerRegion
  }
}


# Merge all data from all sample gene sets into one big dataframe.
allResultsProportions_ClusterAnalysis <- data.frame()

for (test in names(sampleGenesProportions_ClusterAnalysis)) {
  for (cluster in names(sampleGenesProportions_ClusterAnalysis[[test]])) {
    
    df2 <- sampleGenesProportions_ClusterAnalysis[[test]][[cluster]]
    df2 <- cbind(df2, data.frame(SampleGenes = rep(test, times = nrow(df2))))
    
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


