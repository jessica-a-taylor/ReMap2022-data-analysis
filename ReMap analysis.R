if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("karyoploteR","rtracklayer", "TxDb.Athaliana.BioMart.plantsmart28", "ggplot2"))

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

# Create hashes for storing the % R-genes with a chromatin mark in each gene region (frequency)
# and the proportion of each gene region with that mark.
sampleGenesFrequencies <- hash()
sampleGenesProportions <- hash()

# Choose ecotype and tissue for analsis.
# Options: ColLeaf, ColRoot
tissueForAnalysis <- "ColLeaf"

genesForAnalysis <- c("AT1G72840","AT1G72850","AT1G72852","AT1G72860","AT1G72870","AT1G72890",
                      "AT1G72900","AT1G72910","AT1G72920","AT1G72930", "AT1G72940","AT1G72950")


for (test in names(sampleGenes)[grepl("_",names(sampleGenes))]) {
  
  for (level in exLevel) {
    dataToUse <- sampleGenes[[test]][[level]]
    # dataToUse <- sampleGenes[[test]][[level]][c(which(sampleGenes[[test]][[level]]$Gene %in% genesForAnalysis)),]
    
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
    
    modFrequencyPerRegion <- modFrequenciesFunction(geneRegions, allOverlaps, epiMods)
    modProportionPerRegion <- modProportionsFunction(geneRegions, allOverlaps, epiMods)
    
    # Collect all hashes for modFrequencyPerRegion and modProportionPerRegion into single dataframes.
    modFrequencyPerRegion <- mergeResults(modFrequencyPerRegion)
    modProportionPerRegion <- mergeResults(modProportionPerRegion)
    
    # Add a column to modFrequencyPerRegion and modProportionPerRegion with the numbers for 
    # each gene region that will correspond with their position on the x axis.
    modFrequencyPerRegion <- geneRegionAxisLocations(modFrequencyPerRegion, geneRegions)
    modProportionPerRegion <- geneRegionAxisLocations(modProportionPerRegion, geneRegions)
    
    # Add a column to modFrequencyPerRegion and modProportionPerRegion with the current expression level.
    modFrequencyPerRegion <- expressionColumn(modFrequencyPerRegion, level)
    modProportionPerRegion <- expressionColumn(modProportionPerRegion, level)
    
    # Store final results on the appropriate hash.
    sampleGenesFrequencies[[test]][[level]] <- modFrequencyPerRegion
    sampleGenesProportions[[test]][[level]] <- modProportionPerRegion
  }
}

# Merge all data from all sample gene sets into one big dataframe.
allResultsFrequencies <- data.frame()
allResultsProportions <- data.frame()

for (test in names(sampleGenesFrequencies)) {
  for (level in exLevel) {
    df1 <- sampleGenesFrequencies[[test]][[level]]
    df1 <- cbind(df1, data.frame(SampleGenes = rep(test, times = nrow(df1))))
    
    allResultsFrequencies <- rbind(allResultsFrequencies, df1)
    
    df2 <- sampleGenesProportions[[test]][[level]]
    df2 <- cbind(df2, data.frame(SampleGenes = rep(test, times = nrow(df2))))
    
    allResultsProportions <- rbind(allResultsProportions, df2)
  }
}

rm(test, df1, df2, tissueForAnalysis, allOverlaps, modFrequencyPerRegion, modProportionPerRegion, dataToUse)

write.csv(allResultsFrequencies, file = "allResultsFrequencies.csv")
write.csv(allResultsProportions, file = "allResultsProportions.csv")