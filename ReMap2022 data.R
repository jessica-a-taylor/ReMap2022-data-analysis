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
for (analysis in c("PlantExp data", "RNA-seq data")) {
  for (tissue in c("leaves", "root", "seedlings")) {
    jobRunScript("T.test for enrichment in R-genes and controls.R", name = paste("Enrichment_", tissue, sep = ""), importEnv = TRUE)
  }
}



# T-Test - is there a significant difference in the average proportion of coverage of  
# each gene region by a particular modification between R-genes at each expression level?
# Plots comparing the average proportion of coverage of each gene region by a particular 
# modification in R-genes at each expression level.
for (analysis in c("PlantExp data", "RNA-seq data")) {
  for (tissue in c("leaves", "root", "seedlings")) {
    jobRunScript("T.test for enrichment in R-genes only.R", name = paste("Enrichment_", analysis, "_", tissue, sep = ""), importEnv = TRUE)
  }
}

# For each chromatin modification, plot a bar graph of the enrichment in each R-gene.
for (analysis in c("PlantExp data", "RNA-seq data")) {
  for (tissue in c("leaves", "root", "seedlings")) {
  jobRunScript("Enrichment per gene.R", name = paste("Enrichment_", analysis, "_", tissue, sep = ""), importEnv = TRUE)
  }
}


source("Functions\\Overlaps functions.R")
source("Functions\\Modifications per gene.R")
source("Functions\\Coordinates per gene region.R")
source("Functions\\Modification frequencies & proportions.R")
source("Functions\\Get range - merge gene coordinates.R")

# Which R-genes overlap with ACRs?

# Import the expression data from the ACRs paper (Ding et al. 2021).
Ding_Control_ACR <- as.data.frame(read_xlsx("Data\\ACRs data Ding et al., 2021.xlsx", sheet = 2))
Ding_ETI_ACR <- as.data.frame(read_xlsx("Data\\ACRs data Ding et al., 2021.xlsx", sheet = 3))

# Filter for R-genes.
Ding_Control_ACR <- Ding_Control_ACR[which(Ding_Control_ACR$geneId %in% sampleGenes[["NLRs"]]$Gene),-c(1,4,5,8,9,10,13,14)]
Ding_ETI_ACR <- Ding_ETI_ACR[which(Ding_ETI_ACR$geneId %in% sampleGenes[["NLRs"]]$Gene),-c(1,4,5,8,9,10,13,14)]

# Merge start and end coordinates columns to create a ranges column.
source("Functions\\Get range - merge gene coordinates.R")

Ding_Control_ACR$ranges <- mergeCoordinates(Ding_Control_ACR)
Ding_Control_ACR$width <- Ding_Control_ACR$end - Ding_Control_ACR$start
colnames(Ding_Control_ACR) <- c("start", "end", "annotation", "seqnames", "strand", "Gene", "ranges", "width")

for (row in 1:nrow(Ding_Control_ACR)) {
  if (Ding_Control_ACR[row, "strand"] == 1) {
    Ding_Control_ACR[row, "strand"] <- "+"
  } else if (Ding_Control_ACR[row, "strand"] == 2) {
    Ding_Control_ACR[row, "strand"] <- "-"
  }
}

Ding_ETI_ACR$ranges <- mergeCoordinates(Ding_ETI_ACR)
Ding_ETI_ACR$width <- Ding_ETI_ACR$end - Ding_ETI_ACR$start
colnames(Ding_ETI_ACR) <- c("start", "end", "annotation", "seqnames", "strand", "Gene", "ranges", "width")

for (row in 1:nrow(Ding_ETI_ACR)) {
  if (Ding_ETI_ACR[row, "strand"] == 1) {
    Ding_ETI_ACR[row, "strand"] <- "+"
  } else if (Ding_ETI_ACR[row, "strand"] == 2) {
    Ding_ETI_ACR[row, "strand"] <- "-"
  }
}


# Determine which ACRs overlap with the R-genes.
Ding_Control_ACR <- cbind(Ding_Control_ACR, data.frame(Condition = rep("Control", times = nrow(Ding_Control_ACR))))
Ding_ETI_ACR <- cbind(Ding_ETI_ACR, data.frame(Condition = rep("ETI", times = nrow(Ding_ETI_ACR))))

ACR_data <- Ding_Control_ACR
ACR_data <- rbind(ACR_data, Ding_ETI_ACR)

source("Functions\\Coordinates per gene region.R")
ACR_GeneRegions <- getGeneCoordinates(sampleGenes[["NLRs"]])
`%!in%` <- Negate(`%in%`)

ACR_regions <- c() 
for (row in 1:nrow(ACR_data)) {
  overlappingACRs <- data.frame()
  
  for (n in names(ACR_GeneRegions)) {
    df <- ACR_GeneRegions[[n]][ACR_GeneRegions[[n]]$Gene == ACR_data[row, "Gene"],]
    
    if (nrow(df) != 0) {
      if (overlapsFunction(ACR_data[row, "start"], ACR_data[row, "end"], 
                           df[, "start"], df[, "end"])==TRUE) {
        
        overlappingACRs <- rbind(overlappingACRs, data.frame(Region = n,
                                                             Overlap = newOverlapsFunction(as.numeric(ACR_data[row, "start"]), as.numeric(ACR_data[row, "end"]),
                                                                                           as.numeric(df[, "start"]), as.numeric(df[, "end"]))))
        
      } else overlappingACRs <- overlappingACRs
    } else overlappingACRs <- overlappingACRs
  }
  if ((nrow(overlappingACRs)==0 & row == 84)==TRUE) {
    ACR_regions <- append(ACR_regions, "Downstream (within neighbouring gene)")
  }
  if ((nrow(overlappingACRs)==0 & row != 84 & row %in% c(17,149,203))==TRUE) {
    ACR_regions <- append(ACR_regions, "UpstreamIntergenic")
  }
  if ((nrow(overlappingACRs)==0 & row %!in% c(17,149,203, 84) & row %in% c(122,262,376,377,410,433:435))==TRUE) {
    ACR_regions <- append(ACR_regions, "UpstreamIntergenic (within neighbouring gene)")
  }
  ACR_regions <- append(ACR_regions, overlappingACRs[c(which(overlappingACRs$Overlap == max(overlappingACRs$Overlap)))[1], "Region"])
}

ACR_data <- cbind(ACR_data, data.frame(Region = ACR_regions))
Ding_Control_ACR <- ACR_data[which(ACR_data$Condition=="Control"),]
Ding_ETI_ACR <- ACR_data[which(ACR_data$Condition=="ETI"),]

# Which R-genes are associated with the nuclear envelope?
Leaf_NE_data <- as.data.frame(read_xlsx("Data\\Bi et al. (2017) NE associations.xlsx", sheet = 2))
Root_NE_data <- as.data.frame(read_xlsx("Data\\Bi et al. (2017) NE associations.xlsx", sheet = 1))

Leaf_NE_data <- cbind(Leaf_NE_data, data.frame(Tissue = rep("Leaf", times = nrow(Leaf_NE_data))))
Root_NE_data <- cbind(Root_NE_data, data.frame(Tissue = rep("Root", times = nrow(Root_NE_data))))

NE_data <- Leaf_NE_data
NE_data <- rbind(NE_data, Root_NE_data)

geneRegions <- getGeneCoordinates(sampleGenes[["NLRs"]])

NE_associations <- data.frame()

for (r in names(geneRegions)) {
  for (i in 1:nrow(NE_data)) {
    for (j in 1:nrow(geneRegions[[r]]))
      
      if (NE_data[i, "seqnames"]==geneRegions[[r]][j, "seqnames"] & 
          overlapsFunction(NE_data[i, "start"], NE_data[i, "end"],
                           geneRegions[[r]][j, "start"], geneRegions[[r]][j, "end"])==TRUE){
        
        NE_associations <- rbind(NE_associations, data.frame(seqnames = geneRegions[[r]][j, "seqnames"],
                                                             Gene = geneRegions[[r]][j, "Gene"],
                                                             start = NE_data[i, "start"],
                                                             end = NE_data[i, "end"],
                                                             Region = r,
                                                             Tissue = NE_data[i, "Tissue"]))
      }
    else NE_associations <- NE_associations
  }
}

NE_associations$ranges <- mergeCoordinates(NE_associations)

Leaf_NE_data <- NE_associations[which(NE_associations$Tissue == "Leaf"),]
Root_NE_data <- NE_associations[which(NE_associations$Tissue == "Root"),]

# Add sheets to 'Arabidopsis NLRs' spreadsheet giving a summary of the enrichment of chromatin modifications in 
# each gene region for each tissue.
wb <- loadWorkbook("Data\\Arabidopsis NLRs.xlsx")

for (mod in c("H3K9me2","H3K27me3","H2A-Z","H2AK121ub","H3K4me3","H3K36me3","H3K27ac","H3K9ac")) {
  template <- as.data.frame(read_xlsx("Data\\Arabidopsis NLRs.xlsx", sheet = 2))
  
  for (tissue in c("leaves", "root", "seedlings")) {
    df <- resultsProportions[["RNA-seq data"]][grepl("NLRs", resultsProportions[["RNA-seq data"]]$dataToAnalyse) &
                                           grepl(tissue, resultsProportions[["RNA-seq data"]]$dataToAnalyse) &
                                           grepl(mod, resultsProportions[["RNA-seq data"]]$Modification),]
    
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
      
      template[,which(grepl(r, names(template)) & grepl(tissue, names(template)) & grepl("Enrichment", names(template)))] <- df1$Proportion
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




