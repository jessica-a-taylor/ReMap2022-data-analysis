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
library(openxlsx)


# Import coordinates of the genomic regions of interest.
# Options: euchromaticRegions, heterochromaticRegions, euchromaticWithoutTEs

genomicData <- as.data.frame(read_csv("Data\\Protein coding genes.csv"))
genomicData <- genomicData[,-1]

# Get coordinates for 10 sets of random genes and store in a hash.
source("Functions\\Sample random genes.R")

if (TRUE %in% grepl("control", list.files(path = paste("Data\\RNA-seq data\\"),pattern="*.txt"))) {
  sampleGenes <- existingSets(genomicData)
} else sampleGenes <- geneSets(genomicData)

rm(geneSets, existingSets)


# Import list of R-genes.
ArabidopsisNLRs <- as.data.frame(read_xlsx("Data\\Arabidopsis NLRs.xlsx", sheet = 1))
clusteredNLRs <- ArabidopsisNLRs[grepl("cluster", ArabidopsisNLRs$Clustering),]
notClusteredNLRs <- ArabidopsisNLRs[c(which(ArabidopsisNLRs$Clustering =="single")),]

NLRgenes <- genomicData[which(genomicData$Gene %in% ArabidopsisNLRs$Gene),]
NLRgenes <- cbind(NLRgenes, 
                  data.frame(Clustering = ArabidopsisNLRs[which(ArabidopsisNLRs$Gene %in% NLRgenes$Gene),"Clustering"]))


clusteredNLRgenes <- genomicData[which(genomicData$Gene %in% clusteredNLRs$Gene),]
clusteredNLRgenes <- cbind(clusteredNLRgenes, 
                           data.frame(Clustering = clusteredNLRs[which(clusteredNLRs$Gene %in% clusteredNLRgenes$Gene),"Clustering"]))

notClusteredNLRgenes <- genomicData[which(genomicData$Gene %in% notClusteredNLRs$Gene),]
notClusteredNLRgenes <- cbind(notClusteredNLRgenes, 
                              data.frame(Clustering = notClusteredNLRs[which(notClusteredNLRs$Gene %in% notClusteredNLRgenes$Gene),"Clustering"]))

# Add R-genes to sampleGenes.
sampleGenes[["NLRs"]] <- NLRgenes
sampleGenes[["clusteredNLRs"]] <- clusteredNLRgenes
sampleGenes[["notClusteredNLRs"]] <- notClusteredNLRgenes

rm(ArabidopsisNLRs, NLRgenes)

# Get filtered expression data for each set of sample genes in each tissue. 
# Add dataframes to new sampleGenes hashes for gene sets with particular expression levels.
exLevel <- c("No Expression", "Low Expression", "Intermediate Expression",
             "High Expression")

source("Functions\\PlantExp.R")

sampleGenesPlantExp <- PlantExp(sampleGenes[c("control1","control10","control2","control3","control4",
                                              "control5","control6","control7","control8","control9","NLRs")], exLevel)

source("Functions\\RNA-seq data.R")
sampleGenesRNAseq <- RNA_seqAnalysis(sampleGenes[c("control1","control10","control2","control3","control4",
                                                   "control5","control6","control7","control8","control9","NLRs")], exLevel)


rm(PlantExp, RNA_seqAnalysis)

# Use ReMap2022 data to analyse the enrichment of chromatin marks on the R-genes and controls.
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

rm(sampleGenesPlantExp, sampleGenesRNAseq)


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

rm(allResultsFrequencies, allResultsProportions, allResultsAverageProportions)

axisText <- c("intergenic", "Promotor \n(1kb)", "Promotor \n(500bp)", "TSS", "20%",
              "40%", "60%", "80%", "100%", "Downstream \n(200bp)", "Intergenic")

# Fisher's Exact Test - are R-genes enriched amongst those that possess a particular chromatin modification?
# Plots comparing the occurrence of chromatin modifications in the seedlings of R-genes and controls.
for (tissue in c("leaves", "root", "seedlings")) {
  jobRunScript("Fisher's test for enrichment.R", name = paste("Enrichment_", tissue, sep = ""), importEnv = TRUE)
}


# T-Test - is there a significant difference in the average proportion of coverage of  
# each gene region by a particular modification between R-genes and controls?
# Plots comparing the average proportion of coverage of each gene region by a particular modification in R-genes and controls.
for (tissue in c("leaves", "root", "seedlings")) {
    jobRunScript("T.test for enrichment in R-genes and controls.R", name = paste("Enrichment_", tissue, sep = ""), importEnv = TRUE)
}



# T-Test - is there a significant difference in the average proportion of coverage of  
# each gene region by a particular modification between R-genes at each expression level?
# Plots comparing the average proportion of coverage of each gene region by a particular 
# modification in R-genes at each expression level.
for (analysis in c("PlantExp data", "RNA-seq data")) {
  for (tissue in c("leaves", "root", "seedlings")) {
    jobRunScript("T.test for enrichment in R-genes only.R", name = paste("Enrichment_", tissue, sep = ""), importEnv = TRUE)
  }
}

# For each chromatin modification, plot a bar graph of the enrichment in each R-gene.
for (analysis in c("PlantExp data", "RNA-seq data")) {
  for (tissue in c("leaves", "root", "seedlings")) {
  jobRunScript("Enrichment per gene.R", name = paste("Enrichment_", analysis, "_", tissue, sep = ""), importEnv = TRUE)
  }
}


# Add sheets to 'Arabidopsis NLRs' spreadsheet giving a summary of the enrichment of chromatin modifications in 
# each gene region for each tissue.
wb <- loadWorkbook("Data\\Arabidopsis NLRs.xlsx")

for (mod in c("H3K9me2","H3K27me3","H2A-Z","H2AK121ub","H3K4me3","H3K36me3","H3K27ac","H3K9ac")) {
  template <- as.data.frame(read_xlsx("Data\\Arabidopsis NLRs.xlsx", sheet = 4))
  
  for (tissue in c("leaves", "root", "seedlings")) {
    df <- resultsProportions[["PlantExp data"]][grepl("NLRs", resultsProportions[["PlantExp data"]]$dataToAnalyse) &
                                           grepl(tissue, resultsProportions[["PlantExp data"]]$dataToAnalyse) &
                                           grepl(mod, resultsProportions[["PlantExp data"]]$Modification),]
    
    df <- df[order(df$Gene),]
    
    Leaf_NEAs <- c()
    Root_NEAs <- c()
    
    for (gene in unique(df$Gene)) {
      overlappingLeafNEAs <- Leaf_NE_data[Leaf_NE_data$Gene==gene,]
      
      if (nrow(overlappingLeafNEAs) != 0) {
        Leaf_NEAs <- append(Leaf_NEAs, "Yes")
      } else Leaf_NEAs <- append(Leaf_NEAs, "No")
      
      overlappingRootNEAs <- Root_NE_data[Root_NE_data$Gene==gene,]
      
      if (nrow(overlappingRootNEAs) != 0) {
        Root_NEAs <- append(Root_NEAs, "Yes")
      } else Root_NEAs <- append(Root_NEAs, "No")
    }
    
    for (r in unique(df$Region)) {
      
      df1 <- df[df$Region == r,]
      
      control_ACR <- c()
      ETI_ACR <- c()
      
      for (row in 1:nrow(df1)) {
        
        overlappingACRs <- Ding_Control_ACR[which(Ding_Control_ACR$Gene==df1[row,"Gene"] & Ding_Control_ACR$Region==df1[row,"Region"]),]
        
        if (nrow(overlappingACRs) != 0) {
          control_ACR <- append(control_ACR, "Yes")
        } else control_ACR <- append(control_ACR, "No")
        
        overlappingACRs <- Ding_ETI_ACR[which(Ding_ETI_ACR$Gene==df1[row,"Gene"] & Ding_ETI_ACR$Region==df1[row,"Region"]),]
        
        if (nrow(overlappingACRs) != 0) {
          ETI_ACR <- append(ETI_ACR, "Yes")
        } else ETI_ACR <- append(ETI_ACR, "No")
      }
      
      template[,which(grepl(r, names(template)) & grepl(tissue, names(template)))] <- df1$Proportion
      template[,which(grepl("Control_ACR", names(template)) & grepl(r, names(template)) & grepl(tissue, names(template)))] <- control_ACR
      template[,which(grepl("ETI_ACR", names(template)) & grepl(r, names(template)) & grepl(tissue, names(template)))] <- ETI_ACR
    }
    template[,which(grepl("Leaf_NE_association", names(template)) & grepl(tissue, names(template)))] <- Leaf_NEAs
    template[,which(grepl("Root_NE_association", names(template)) & grepl(tissue, names(template)))] <- Root_NEAs
  }
  addWorksheet(wb,mod)
  writeData(wb,mod,template)
  saveWorkbook(wb,"Data\\Arabidopsis NLRs.xlsx", overwrite = TRUE)
}

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




