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

# Create hashes for storing the % R-genes with a chromatin mark in each gene region (frequency)
# and the proportion of each gene region with that mark.
LeafFrequencies <- hash()
LeafProportions <- hash()

# Choose ecotype and tissue for analsis.
# Options: ColLeaf, ColRoot
tissueForAnalysis <- "ColLeaf"


for (test in names(sampleGenes)[grepl("_Leaf",names(sampleGenes))]) {
  
  for (level in exLevel) {
    dataToUse <- sampleGenes[[test]][[level]]

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
    
    # Add a column to modFrequencyPerRegion and modProportionPerRegion with the numbers for 
    # each gene region that will correspond with their position on the x axis.
    modFrequencyPerRegion <- geneRegionAxisLocations(modFrequencyPerRegion, geneRegions)
    modProportionPerRegion <- geneRegionAxisLocations(modProportionPerRegion, geneRegions)
    
    # Add a column to modFrequencyPerRegion and modProportionPerRegion with the current expression level.
    modFrequencyPerRegion <- expressionColumn(modFrequencyPerRegion, level)
    modProportionPerRegion <- expressionColumn(modProportionPerRegion, level)
    
    # Store final results on the appropriate hash.
    LeafFrequencies[[test]][[level]] <- modFrequencyPerRegion
    LeafProportions[[test]][[level]] <- modProportionPerRegion
  }
}

# Merge all data from all sample gene sets into one big dataframe.
allLeafFrequencies <- data.frame()
allLeafProportions <- data.frame()

for (test in names(LeafFrequencies)) {
  for (level in exLevel) {
    df1 <- LeafFrequencies[[test]][[level]]
    df1 <- cbind(df1, data.frame(SampleGenes = rep(test, times = nrow(df1))))
    
    allLeafFrequencies <- rbind(allLeafFrequencies, df1)
    
    df2 <- LeafProportions[[test]][[level]]
    df2 <- cbind(df2, data.frame(SampleGenes = rep(test, times = nrow(df2))))
    
    allLeafProportions <- rbind(allLeafProportions, df2)
  }
}

rm(test, df1, df2, tissueForAnalysis, allOverlaps, modFrequencyPerRegion, modProportionPerRegion, dataToUse)
