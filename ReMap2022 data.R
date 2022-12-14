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


# Import coordinates of the genomic regions of interest.
# Options: euchromaticRegions, heterochromaticRegions, euchromaticWithoutTEs

dataToUse <- as.data.frame(read_csv("Data\\euchromaticWithoutTEs.csv"))

# Get 10 sets of random genes and store in a hash from gene dataset of interest.
source("Functions\\Sample random genes.R")

if (TRUE %in% grepl("control", list.files(path = paste("Data\\RNA-seq data\\"),pattern="*.txt"))) {
  sampleGenes <- existingSets(dataToUse)
} else sampleGenes <- geneSets(dataToUse)


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
# Add dataframes to new sampleGenes hashes for gene sets with particular expression levels.
source("Functions\\PlantExp.R")

exLevel <- c("No Expression", "Low Expression", "Intermediate Expression",
             "High Expression")

sampleGenesPlantExp <- PlantExp(sampleGenes[c("control1","control10","control2","control3","control4",
                                              "control5","control6","control7","control8","control9","NLRs")], exLevel)

source("Functions\\RNA-seq data.R")
sampleGenesRNAseq <- RNA_seqAnalysis(sampleGenes[c("control1","control10","control2","control3","control4",
                                                   "control5","control6","control7","control8","control9","NLRs")], exLevel)

rm(gene, row, test)

# Use ReMap2022 data to analyse the enrichment of chromatin marks on the R-genes and controls.
source("Functions\\Modifications per gene.R")
source("Functions\\Coordinates per gene region.R")
source("Functions\\Modification frequencies & proportions.R")
source("Functions\\Get range - merge gene coordinates.R")

# Import filtered ReMap2022 data.
ReMap <- as.data.frame(read_xlsx("Data\\Filtered ReMap data.xlsx"))

# Create list of chromatin modifications.
epiMods <- unique(ReMap$epiMod)

# For each gene set (R-genes and controls), tissue, and expression level, 
# find the percentage of genes with each chromatin modification (frequency)
# and the enrichment of the modification (proportion) in each gene region.

# Run parallel analyses as background jobs.
for (analysis in c("PlantExp data", "RNA-seq data")) {
  for (tissue in c("leaves", "root", "seedlings")) {
  jobRunScript("ReMap analysis.R", name = paste(analysis, tissue, sep = "_"), importEnv = TRUE)
  }
}

rm(betweenFunction, expressionColumn, findItem, geneRegionAxisLocations, geneSets, getGeneCoordinates,
   mergeCoordinates, mergeOverlappingModifications, modFrequenciesFunction, modificationOccurrences,
   modProportionsFunction, overlapsFunction, ReMapPerGene, PlantExp, RNA_seqAnalysis, newOverlapsFunction)

# Import the results.
resultsFrequencies <- hash()
resultsProportions <- hash()
resultsAverageProportions <- hash()

for (analysis in c("PlantExp data", "RNA-seq data")) {
  allResultsFrequencies <- data.frame()
  allResultsProportions <- data.frame()
  allResultsAverageProportions <- data.frame()
  
  for (tissue in c("leaves", "root", "seedlings")) {
    allResultsFrequencies <- rbind(allResultsFrequencies, as.data.frame(read_csv(paste("Data\\", analysis, "\\", tissue,"\\allResultsFrequencies.csv", sep = ""))))
    allResultsProportions <- rbind(allResultsProportions, as.data.frame(read_csv(paste("Data\\", analysis, "\\", tissue,"\\allResultsProportions.csv", sep = ""))))
    allResultsAverageProportions <- rbind(allResultsAverageProportions, as.data.frame(read_csv(paste("Data\\", analysis, "\\", tissue,"\\allResultsAverageProportions.csv", sep = ""))))
  }
  resultsFrequencies[[analysis]] <- allResultsFrequencies
  resultsProportions[[analysis]] <- allResultsProportions
  resultsAverageProportions[[analysis]] <- allResultsAverageProportions
}

axisText <- c("intergenic", "Promotor \n(1kb)", "Promotor \n(500bp)", "TSS", "20%",
              "40%", "60%", "80%", "100%", "Downstream \n(200bp)", "Intergenic")

# Fisher's Exact Test - are R-genes enriched amongst those that possess a particular chromatin modification?
# Plots comparing the occurrence of chromatin modifications in the seedlings of R-genes and controls.
for (tissue in c("leaves", "root", "seedlings")) {
  jobRunScript("Fisher's test for enrichment.R", name = paste("Enrichment_", tissue, sep = ""), importEnv = TRUE)
}


# T-Test - is there a significant difference in the average proportion of coverage of  
# each gene region by a particular modification between R-genes and controls?
# Plots comparing the average proportion of coverage of each gene region by a particular modification in the 
# seedlings of R-genes and controls.
for (tissue in c("leaves", "root", "seedlings")) {
    jobRunScript("T.test for enrichment.R", name = paste("Enrichment_", tissue, sep = ""), importEnv = TRUE)
}


# For each chromatin modification, plot the enrichment in each R-gene for each expression level.
for (analysis in c("PlantExp data", "RNA-seq data")) {
  for (tissue in c("leaves", "root", "seedlings")) {
    df <- resultsAverageProportions[[analysis]][grepl("NLRs", resultsAverageProportions[[analysis]]$Tissue) &
                                                  grepl(tissue, resultsAverageProportions[[analysis]]$Tissue) &
                                                  !grepl("luster", resultsAverageProportions[[analysis]]$Tissue),]
    
    for (mod in unique(df$Modification)) {
      df1 <- df[df$Modification==mod,]
      
      plot <- ggplot(df1, aes(x = axisGroup, y = Proportion, color = Expression)) +
        geom_line(aes(group = Expression),linewidth = 1) +
        scale_x_continuous(limits = c(-60, 140), breaks = seq(-60, 140, 20), labels = axisText) +
        scale_y_continuous(limits = c(0,1), expand = c(0,0)) + 
        labs(x = "", y = "Average proportion of gene region", color = "Expression level", title = mod) +
        geom_vline(xintercept=0, color="grey", linewidth=1) +
        scale_color_manual(values = c("coral2", "darkslategrey")) +        
        coord_cartesian(ylim= c(0,1), clip = "off") + theme(plot.margin = unit(c(1,1,.3,1), "lines")) +
        annotation_custom(textGrob("% of gene length from TSS", gp=gpar(fontsize=16, col = "grey33")),xmin=0,xmax=100,ymin=-.17,ymax=-.17) + 
        annotation_custom(textGrob("Gene region", gp=gpar(fontsize=16)),xmin=0,xmax=100,ymin=-.25,ymax=-.25) +
        theme(axis.text.x = element_text(size = 14, colour = "black", angle = 45, vjust = 1, hjust = 1), axis.text.y = element_text(size = 14,colour = "black"), 
              axis.title.y = element_text(size = 16, vjust = 2), plot.title = element_text(hjust = .5, size = 16),
              legend.text = element_text(size = 12), legend.title = element_text(size = 14), axis.line = element_line(linewidth = .6))
      
      ggsave(paste("Graphs\\Enrichment\\", analysis, "\\", tissue, "\\R-gene enrichment for ", mod, ".pdf", sep = ""), plot = plot, width = 12, height = 6)
      
    }
  }
}

# For each chromatin modification, plot a bar graph of the enrichment in each R-gene.




# Determine whether the enrichment of each chromatin mark is significantly more similar between genes within a cluster 
# than between all genes.
for (analysis in c("PlantExp data", "RNA-seq data")) {
  for (tissue in c("leaves", "root", "seedlings")) {
    jobRunScript("Enrichment variability within clusters.R", name = tissue, importEnv = TRUE)
  }
}







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




