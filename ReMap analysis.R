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

source("Functions\\Overlaps functions.R")
source("Functions\\Modifications per gene.R")
source("Functions\\Coordinates per gene region.R")
source("Functions\\Modification frequencies & proportions.R")

if (analysis == "PlantExp data") {
  dataToAnalyse <- sampleGenesPlantExp
} else if (analysis == "RNA-seq data") {
  dataToAnalyse <- sampleGenesRNAseq
}

# Create hashes for storing the % R-genes with a chromatin mark in each gene region (frequency)
# and the proportion of each gene region with that mark.
dataToAnalyseFrequencies <- hash()
dataToAnalyseProportions <- hash()

# Generate a list of the number of genes in each set.
geneCount <- data.frame()

# Choose ecotype and tissue for dataToAnalyse.
# Options: leafGenes, rootGenes, seedlingGenes
if (tissue == "leaves") {
  tissueForAnalysis <- "leafGenes"
} else if (tissue == "root") {
  tissueForAnalysis <- "rootGenes"
} else if (tissue == "seedlings") {
  tissueForAnalysis <- "seedlingGenes"
}


for (test in names(dataToAnalyse)[grepl(paste("_", tissue, sep = ""),names(dataToAnalyse))]) {
  
  for (level in exLevel) {
    dataToUse <- dataToAnalyse[[test]][[level]]

    # Create a hash with the ReMap data in a particular tissue for the current set of genes. 
    allModifications <- ReMapPerGene(dataToUse, tissueForAnalysis)
    
    # For each gene in the current set of genes, create a new hash with the occurrences of each chromatin modification.
    geneModifications <- modificationOccurrences(allModifications)
    
    rm(allModifications)
    
    # For each gene in the current set of genes, merge the overlapping occurrences of each modification.
    allOverlaps <- mergeOverlappingModifications(geneModifications)
    
    geneCount <- rbind(geneCount, data.frame(GeneSet = paste(test, "_", level, sep = ""),
                                             GeneCount = length(names(allOverlaps))))
    print(length(names(allOverlaps)))
    
    rm(geneModifications)
    
    # Determine the % R-genes with a chromatin mark in each gene region (frequency)
    # and the proportion of each gene region with that mark.
    geneRegions <- getGeneCoordinates(dataToUse)
    
    modFrequencyPerRegion <- modFrequenciesFunction(geneRegions, allOverlaps, epiMods)
    modProportionPerRegion <- modProportionsFunction(geneRegions, allOverlaps, epiMods)
    
    # Add a column to modFrequencyPerRegion and modProportionPerRegion with the numbers for 
    # each gene region that will correspond with their position on the x axis.
    modFrequencyPerRegion <- geneRegionAxisLocations(modFrequencyPerRegion, geneRegions)
    modProportionPerRegion <- geneRegionAxisLocations(modProportionPerRegion, geneRegions)
    
    rm(geneRegions)
    
    # Add a column to modFrequencyPerRegion and modProportionPerRegion with the current expression level.
    modFrequencyPerRegion <- expressionColumn(modFrequencyPerRegion, level)
    modProportionPerRegion <- expressionColumn(modProportionPerRegion, level)
    
    # Store final results on the appropriate hash.
    dataToAnalyseFrequencies[[test]][[level]] <- modFrequencyPerRegion
    dataToAnalyseProportions[[test]][[level]] <- modProportionPerRegion
    
    print(level)
  }
  print(test)
}

write.csv(geneCount, file = paste("Data\\", analysis, "\\", tissue, "\\Gene count.txt", sep=""))

# Merge all data from all sample gene sets into one big dataframe.
allResultsFrequencies <- data.frame()
allResultsProportions <- data.frame()

for (test in names(dataToAnalyseFrequencies)) {
  for (level in exLevel) {
    df1 <- dataToAnalyseFrequencies[[test]][[level]]
    df1 <- cbind(df1, data.frame(dataToAnalyse = rep(test, times = nrow(df1))))
    
    allResultsFrequencies <- rbind(allResultsFrequencies, df1)
    
    df2 <- dataToAnalyseProportions[[test]][[level]]
    df2 <- cbind(df2, data.frame(dataToAnalyse = rep(test, times = nrow(df2))))
    
    allResultsProportions <- rbind(allResultsProportions, df2)
  }
}

write.csv(allResultsFrequencies, file = paste("Data\\", analysis, "\\", tissue,"\\allResultsFrequencies.csv", sep = ""))
write.csv(allResultsProportions, file = paste("Data\\", analysis, "\\", tissue, "\\allResultsProportions.csv", sep = ""))

# Calculate the mean proportion of overlap and add as a new column to the dataframe.
allResultsAverageProportions <- data.frame()

for (test in unique(allResultsProportions$dataToAnalyse)) {
  df <- allResultsProportions[allResultsProportions$dataToAnalyse==test,]
  
  for (level in unique(df$Expression)) { 
    df1 <- df[df$Expression==level,]
    
    if (nrow(df1) >= 1) {
      for (mod in unique(df1$Modification)) {
        df2 <- df1[df1$Modification==mod,]
        
        for (r in unique(df2$Region)) {
          df3 <- df2[df2$Region==r,]
          
            allResultsAverageProportions <- rbind(allResultsAverageProportions, data.frame(Region = r,
                                                                                           Modification = mod,
                                                                                           Proportion = mean(df3$Proportion),
                                                                                           Tissue = df3$dataToAnalyse[1],
                                                                                           axisGroup = df3$axisGroup[1],
                                                                                           Expression = level,
                                                                                           SampleSize = nrow(df3)))

        }
      }
    } else allResultsAverageProportions <- allResultsAverageProportions
  }
}

write.csv(allResultsAverageProportions, file = paste("Data\\", analysis, "\\", tissue, "\\allResultsAverageProportions.csv", sep = ""))


rm(test, df1, df2, tissueForAnalysis, allOverlaps, modFrequencyPerRegion, modProportionPerRegion, dataToUse)