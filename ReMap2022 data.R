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

NLRgenes <- dataToUse[which(dataToUse$Gene %in% ArabidopsisNLRs$Gene),]

# Create a ranges column by merging the start and end columns.
dataToUse <- NLRgenes

NLRgenes$ranges <- mergeCoordinates(dataToUse)


# Add R-genes to sampleGenes.
sampleGenes[["NLRs"]] <- NLRgenes

rm(ArabidopsisNLRs, NLRgenes, Atgenes)


# Get filtered expression data for each set of sample genes in each tissue. 
# Add dataframes to sampleGenes for gene sets with particular expression levels.
source("Functions\\PlantExp.R")

sampleGenes <- PlantExp(sampleGenes)


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

for (test in names(sampleGenes)[c(42:44)]) {

  dataToUse <- sampleGenes[[test]]
  
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
  modFrequencyPerRegion <- expressionColumn(modFrequencyPerRegion)
  modProportionPerRegion <- expressionColumn(modProportionPerRegion)
  
  # Store final results on the appropriate hash.
  sampleGenesFrequencies[[test]] <- modFrequencyPerRegion
  sampleGenesProportions[[test]] <- modProportionPerRegion
}


# Merge all data from all sample gene sets into one big dataframe.
allResultsFrequencies <- data.frame()
allResultsProportions <- data.frame()

for (test in names(sampleGenesFrequencies)) {

  df1 <- sampleGenesFrequencies[[test]]
  df1 <- cbind(df1, data.frame(Test = rep(test, times = nrow(df1))))
  
  allResultsFrequencies <- rbind(allResultsFrequencies, df1)
  
  df2 <- sampleGenesProportions[[test]]
  df2 <- cbind(df2, data.frame(Test = rep(test, times = nrow(df2))))
  
  allResultsProportions <- rbind(allResultsProportions, df2)
}

rm(test, df1, df2, tissueForAnalysis, allOverlaps, modFrequencyPerRegion, modProportionPerRegion, dataToUse, ReMap)


# Plot the the results.
axisText <- c("Intergenic", "Promotor \n(1kb)", "Promotor \n(500bp)", "TSS",
              "20%", "40%", "60%", "80%", "100%", 
              "Downstream \n(200bp)", "Intergenic")

# Frequencies plot for R-genes & controls.
dataToUse <- allResultsFrequencies[c(which(allResultsFrequencies$Test %in% c(names(sampleGenes)[c(1:10,33)]))),]

for (mod in epiMods) {
  df <- dataToUse[dataToUse$Modification==mod,]
  
  plot <- ggplot(df, aes(x = axisGroup, y = Frequency, color = Test)) + 
    scale_x_continuous(limits = c(-60, 140), breaks = seq(-60, 140, 20), labels = axisText) +
    geom_line(aes(group = Test),linewidth = 1.3) +
    geom_point(aes(group = Test), size = 2) + theme_minimal() + 
    scale_colour_manual(limits = c("control1", "NLRs"), 
                        values=c("grey43", "black"), labels = c("Controls", "R-genes")) +
    labs(x = "", y = "Frequency of occurrence (%)") +
    geom_vline(xintercept=0, color="grey", size=1) +
    coord_cartesian(ylim= c(0,100), clip = "off") + theme(plot.margin = unit(c(1,1,2,1), "lines")) +
    annotation_custom(textGrob("% of gene length from TSS", gp=gpar(fontsize=14, col = "grey33")),xmin=0,xmax=100,ymin=-22,ymax=-22) + 
    annotation_custom(textGrob("Gene region", gp=gpar(fontsize=16)),xmin=0,xmax=100,ymin=-30,ymax=-30) +
    theme(axis.text.x = element_text(size = 13, colour = "black", angle = 45, vjust = 1, hjust = 1), axis.text.y = element_text(size = 14,colour = "black"), 
          axis.title.y = element_text(size = 16, vjust = 2)) 
  
  ggsave(paste(mod, "_", ".LeavesFrequencies.pdf", sep = ""), plot = plot, width = 12, height = 6)
}

  
# Frequencies plot for R-genes only.
dataToUse <- allResultsFrequencies[allResultsFrequencies$Tissue=="leafExpression",]

for (mod in epiMods) {
  df <- dataToUse[dataToUse$Modification==mod,]
  
  plot <- ggplot(df, aes(x = axisGroup, y = Frequency, color = factor(Expression, levels = expressionLevel))) + 
    scale_x_continuous(limits = c(-60, 140), breaks = seq(-60, 140, 20), labels = axisText) +
    geom_line(aes(group = Expression),linewidth = 1.3) +
    geom_point(aes(group = Expression), size = 2) + theme_minimal() + 
    scale_colour_brewer(name = "Expression Level", palette = "Reds", direction=-1) +
    labs(x = "", y = "Frequency of occurrence (%)") +
    geom_vline(xintercept=0, color="grey", linewidth=1) +
    coord_cartesian(ylim= c(0,100), clip = "off") + theme(plot.margin = unit(c(1,1,2,1), "lines")) +
    annotation_custom(textGrob("% of gene length from TSS", gp=gpar(fontsize=12, col = "grey33")),xmin=0,xmax=100,ymin=-14,ymax=-14) + 
    annotation_custom(textGrob("Gene region", gp=gpar(fontsize=14)),xmin=0,xmax=100,ymin=-20,ymax=-20) +
    theme(axis.text.x = element_text(size = 11, colour = "black"), axis.text.y = element_text(size = 12,colour = "black"), 
          axis.title.y = element_text(size = 14, vjust = 2)) 

  ggsave(paste(mod, "_", ".LeavesFrequenciesPerExprression.pdf", sep = ""), plot = plot, width = 12, height = 6)
}


# Calculate the mean proportion of overlap and add as a new column to the dataframe.
allResultsAverageProportions <- data.frame()

for (test in unique(allResultsProportions$Test)) {
  df <- allResultsProportions[allResultsProportions$Test==test,]
  
    for (level in unique(allResultsProportions$Expression)) { 
      df1 <- df[df$Expression==level,]
      
      if (nrow(df1) >= 1) {
        for (mod in unique(allResultsProportions$Modification)) {
          df2 <- df1[df1$Modification==mod,]
          
          for (r in unique(allResultsProportions$Region)) {
            df3 <- df2[df2$Region==r,]
            
            allResultsAverageProportions <- rbind(allResultsAverageProportions, data.frame(Region = r,
                                                                                           Modification = mod,
                                                                                           Proportion = mean(df3$Proportion),
                                                                                           Tissue = tissue,
                                                                                           axisGroup = df3$axisGroup[1],
                                                                                           Expression = level,
                                                                                           n = nrow(df3)))
          }
        }
      } else allResultsAverageProportions <- allResultsAverageProportions
    }
}

dataToUse <- allResultsAverageProportions[allResultsAverageProportions$Tissue=="leafExpression",]

# Proportions plot for R-genes only.
for (mod in epiMods) {
  df <- dataToUse[dataToUse$Modification==mod,]
  
    plot <- ggplot(df, aes(x = axisGroup, y = Proportion, color = factor(Expression, levels = expressionLevel))) + 
    scale_x_continuous(limits = c(-60, 140), breaks = seq(-60, 140, 20), labels = axisText) +
    geom_line(aes(x = axisGroup, y = Proportion, group = Expression),size = 1.3) + 
    geom_point(aes(x = axisGroup, y = Proportion),size = 2) +
    scale_colour_brewer(name = "Expression Level", palette = "Reds", direction=-1) +
    theme_minimal() + 
    labs(x = "", y = "Average proportion of gene region") +
    geom_vline(xintercept=0, color="grey", size=1) +
    coord_cartesian(ylim= c(0,1), clip = "off") + theme(plot.margin = unit(c(1,1,2,1), "lines")) +
    annotation_custom(textGrob("% of gene length from TSS", gp=gpar(fontsize=12, col = "grey33")),xmin=0,xmax=100,ymin=-0.15,ymax=-0.15) + 
    annotation_custom(textGrob("Gene region", gp=gpar(fontsize=14)),xmin=0,xmax=100,ymin=-0.2,ymax=-0.2) +
    theme(axis.text.x = element_text(size = 11, colour = "black"), axis.text.y = element_text(size = 12,colour = "black"), 
          axis.title.y = element_text(size = 14, vjust = 2))
  
  ggsave(paste(mod, "_", ".LeafProportionsMeanPerExpression.pdf", sep = ""), plot = plot, width = 12, height = 6)
}


