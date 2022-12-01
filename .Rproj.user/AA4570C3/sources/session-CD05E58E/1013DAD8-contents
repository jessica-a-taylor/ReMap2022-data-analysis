# Create hashes for storing the % R-genes with a chromatin mark in each gene region (frequency)
# and the proportion of each gene region with that mark.
sampleGenesFrequencies <- hash()
sampleGenesProportions <- hash()

# Choose ecotype and tissue for analsis.
# Options: ColLeaf, ColRoot
tissueForAnalysis <- "ColLeaf"

genesForAnalysis <- c("AT1G72840","AT1G72850","AT1G72852","AT1G72860","AT1G72870","AT1G72890",
                      "AT1G72900","AT1G72910","AT1G72920","AT1G72930", "AT1G72940","AT1G72950")


for (test in names(sampleGenes)[grepl(paste("_", tissue, sep = ""),names(sampleGenes))]) {
  
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


write.csv(allResultsFrequencies, file = paste(tissue,"_allResultsFrequencies.csv", sep = ""))
write.csv(allResultsProportions, file = paste(tissue,"_allResultsProportions.csv", sep = ""))