source("Functions\\Overlaps functions.R")
source("Functions\\Modifications per gene.R")
source("Functions\\Coordinates per gene region.R")
source("Functions\\Modification frequencies & proportions.R")
source("Functions\\Get range - merge gene coordinates.R")


# Import the expression data from the ACRs paper (Ding et al. 2021).
Ding_ExpressionData <- as.data.frame(read_xlsx("Data\\ACRs Ding et al., 2021.xlsx", sheet = 1))

# Filter for R-genes.
Ding_ExpressionData <- Ding_ExpressionData[which(Ding_ExpressionData$Gene %in% sampleGenes[["NLRs"]] $Gene),]

# Calculate the mean expression between the 3 replicates for the control and ETI-inducung treatments.
controlExpression <- c()
ETIexpression <- c()

for (row in 1:nrow(Ding_ExpressionData)) {
  controlExpression <- append(controlExpression, mean(as.numeric(Ding_ExpressionData[row, 2:4])))
  ETIexpression <- append(ETIexpression, mean(as.numeric(Ding_ExpressionData[row, 5:7])))
}

# Re-write the dataframe with only the average expression data.
Ding_ExpressionData <- data.frame(Gene = Ding_ExpressionData$Gene,
                                  Control = controlExpression,
                                  ETI = ETIexpression)

# Add information needed for the ReMap analysis.
Ding_ExpressionData <- cbind(Ding_ExpressionData, sampleGenes[["NLRs"]][,c(2:8)])

# Perform ReMap analysis to determine the chromatin modification enrichment in each R-gene.
# Determine which genes appear to be the most responsive to infection - are there any trend in their chromatin modifications.

# Import filtered ReMap2022 data.
ReMap <- as.data.frame(read_xlsx("Data\\Filtered ReMap data.xlsx"))

# Create list of chromatin modifications.
epiMods <- unique(ReMap$epiMod)

# Create a hash with the ReMap data in a particular tissue for the current set of genes. 
allModifications <- ReMapPerGene(Ding_ExpressionData, "leafGenes")

# For each gene in the current set of genes, create a new hash with the occurrences of each chromatin modification.
geneModifications <- modificationOccurrences(allModifications)

rm(allModifications)

# For each gene in the current set of genes, merge the overlapping occurrences of each modification.
allOverlaps <- mergeOverlappingModifications(geneModifications)

print(length(names(allOverlaps)))

rm(geneModifications)

# Determine the % R-genes with a chromatin mark in each gene region (frequency)
# and the proportion of each gene region with that mark.
geneRegions <- getGeneCoordinates(Ding_ExpressionData)

modFrequencyPerRegion <- modFrequenciesFunction(geneRegions, allOverlaps, epiMods)
modProportionPerRegion <- modProportionsFunction(geneRegions, allOverlaps, epiMods)

# Add a column to modFrequencyPerRegion and modProportionPerRegion with the numbers for 
# each gene region that will correspond with their position on the x axis.
modFrequencyPerRegion <- geneRegionAxisLocations(modFrequencyPerRegion, geneRegions)
modProportionPerRegion <- geneRegionAxisLocations(modProportionPerRegion, geneRegions)

rm(geneRegions)

# Which R-genes overlap with ACRs?

# Which TFs are the R-genes associated with - is there a correlation between the enrichment of particlar chromatin 
# modifications and particular TFs?

# Are there similarities in chromatin modification and TF enrichment between co-expressed genes?

# Which R-genes are asssociated with the nuclear envelope?

