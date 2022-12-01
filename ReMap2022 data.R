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
for (tissue in c("Leaf", "Root", "Seedling")) {
  df <- allResultsFrequencies[grepl(tissue, allResultsFrequencies$SampleGenes),]
  
  hypergeometricTest <- data.frame()
  
  for (mod in unique(allResultsFrequencies$Modification)) {
    df1 <- df[df$Modification==mod,]
    
    for (r in unique(allResultsFrequencies$Region)) {
      df2 <- df1[df1$Region==r,]
      
      for (level in unique(allResultsFrequencies$Expression)) {
        df3 <- df2[df2$Expression==level,]
        df3 <- df3[!grepl("luster", df3$SampleGenes),]
        
        # Create a list of genes from each each sample set.
        geneList <- c()
        
        for (row in 1:nrow(df3)) {
          geneList <- append(geneList, rep(df3[row,"SampleGenes"], times = df3[row, "n"]))
        }
        
        meanNLRs <- c()
        for (n in 1:10) {
          geneSample <- sample(geneList, length(geneList)*0.1)
          meanNLRs <- append(meanNLRs, length(geneSample[grepl("NLRs", geneSample)]))
        }
        
        meanNLRs <- mean(meanNLRs)
        
        if (length(geneSample) > 1) {
          
          statTest <- fisher.test(matrix(c(meanNLRs, length(geneList[grepl("NLRs", geneList)])-meanNLRs,
                                           length(geneSample) - meanNLRs, length(geneList[grepl("control", geneList)])-length(geneSample) - meanNLRs), 2,2), 
                                  alternative = "less")
          
          hypergeometricTest <- rbind(hypergeometricTest, data.frame(Expression = level,
                                                                     Modification = mod,
                                                                     Region = r,
                                                                     Sample.R.genes = meanNLRs,
                                                                     Sample.all.genes = length(geneSample),
                                                                     All.R.genes = length(geneList[grepl("NLRs", geneList)]),
                                                                     All.genes = length(geneList),
                                                                     p.value = statTest$p.value))
        }
      }
    }
  }
write.csv(hypergeometricTest, file = paste("Tests\\", tissue, "_Fisher.Test_frequencies.csv", sep = ""))
}

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

for (test in unique(allResultsProportions$SampleGenes)) {
  df <- allResultsProportions[allResultsProportions$SampleGenes==test,]
  
    for (level in unique(allResultsProportions$Expression)) { 
      df1 <- df[df$Expression==level,]
      
      if (nrow(df1) >= 1) {
        for (mod in unique(allResultsProportions$Modification)) {
          df2 <- df1[df1$Modification==mod,]
          
          for (r in unique(allResultsProportions$Region)) {
            df3 <- df2[df2$Region==r,]
            
            if (nrow(df3) >= 10) {
              allResultsAverageProportions <- rbind(allResultsAverageProportions, data.frame(Region = r,
                                                                                             Modification = mod,
                                                                                             Proportion = mean(df3$Proportion),
                                                                                             Tissue = df3$SampleGenes[1],
                                                                                             axisGroup = df3$axisGroup[1],
                                                                                             Expression = level,
                                                                                             SampleSize = paste(level, 
                                                                                                                paste("(n = ", nrow(df3), ")", sep = ""), sep = " ")))
            }
            else allResultsAverageProportions <- allResultsAverageProportions
          }
        }
      } else allResultsAverageProportions <- allResultsAverageProportions
    }
}


# Is there a significant difference in the average proportion of coverage of each gene region by a 
# particular modification between R-genes and controls?
controlSets <- c("control1","control2","control3","control4","control5",
                 "control6","control7","control8","control9","control10")

wilcoxGeneSets <- hash()
ksGeneSets <- hash()

df <- allResultsProportions[grepl("Seedling", allResultsProportions$SampleGenes) & 
                              !grepl("luster", allResultsProportions$SampleGenes),]

for (set in controlSets) {
  wilcoxDF <- data.frame(Modification = character(),
                         Region = character(),
                         W.statistic = numeric(),
                         p.value = numeric())
  
  ksDF <- data.frame(Modification = character(),
                     Region = character(),
                     W.statistic = numeric(),
                     p.value = numeric())
  
  for (mod in unique(allResultsProportions$Modification)) {
    df1 <- df[df$Modification==mod,]
    
    for (r in unique(allResultsProportions$Region)) {
      df2 <- df1[df1$Region==r,]
      
   
      wilcoxTest <- wilcox.test(Proportion~SampleGenes, df2[c(which(df2$SampleGenes == paste(set, "_Seedling", sep = "")), which(df2$SampleGenes == "NLRs_Seedling")),])
      wilcoxDF <- rbind(wilcoxDF, data.frame(Modification = mod,
                                             Region = r,
                                             W.statistic = wilcoxTest$statistic,
                                             p.value = wilcoxTest$p.value))
      
      ksTest <- ks.test(df2[c(which(df2$SampleGenes == paste(set, "_Seedling", sep = ""))),"Proportion"], df2[c(which(df2$SampleGenes == "NLRs_Seedling")),"Proportion"])
      ksDF <- rbind(ksDF, data.frame(Modification = mod,
                                     Region = r,
                                     W.statistic = ksTest$statistic,
                                     p.value = ksTest$p.value))
      
      
    }
  }
  wilcoxGeneSets[[set]] <- wilcoxDF
  ksGeneSets[[set]] <- ksDF
}



dataToUse <- allResultsAverageProportions[c(which(allResultsAverageProportions$Tissue == "control1_Seedling")),]

# Proportions plot for R-genes only.
for (mod in epiMods) {
  df <- dataToUse[dataToUse$Modification==mod,]
  
    plot <- ggplot(df, aes(x = axisGroup, y = Proportion, color = factor(Expression, levels = exLevel))) + 
    scale_x_continuous(limits = c(-60, 140), breaks = seq(-60, 140, 20), labels = axisText) +
    geom_line(aes(x = axisGroup, y = Proportion, group = Expression),size = 1.3) + 
    geom_point(aes(x = axisGroup, y = Proportion),size = 2) +
    scale_colour_brewer(name = "Expression Level", palette = "Reds", direction=1, labels = c(unique(df$SampleSize))) +
    theme_minimal() + 
    labs(x = "", y = "Average proportion of gene region") +
    geom_vline(xintercept=0, color="grey", size=1) +
    coord_cartesian(ylim= c(0,1), clip = "off") + theme(plot.margin = unit(c(1,1,2,1), "lines")) +
    annotation_custom(textGrob("% of gene length from TSS", gp=gpar(fontsize=12, col = "grey33")),xmin=0,xmax=100,ymin=-0.15,ymax=-0.15) + 
    annotation_custom(textGrob("Gene region", gp=gpar(fontsize=14)),xmin=0,xmax=100,ymin=-0.2,ymax=-0.2) +
    theme(axis.text.x = element_text(size = 11, colour = "black"), axis.text.y = element_text(size = 12,colour = "black"), 
          axis.title.y = element_text(size = 14, vjust = 2))
  
  ggsave(paste("NLRs_Seedling", "_",mod , ".pdf", sep = ""), plot = plot, width = 12, height = 6)
}

# Proportions plot for R-genes & controls.
dataToUse <- allResultsAverageProportions[grepl("Seedling", allResultsAverageProportions$SampleGenes) & 
                                     !grepl("luster", allResultsAverageProportions$SampleGenes),]

for (mod in epiMods) {
  df <- dataToUse[dataToUse$Modification==mod,]
  
  for (level in exLevel) {
    df1 <- df[df$Expression==level,]
    
    plot <- ggplot(df1, aes(x = axisGroup, y = Proportion, color = SampleGenes)) + 
      scale_x_continuous(limits = c(-60, 140), breaks = seq(-60, 140, 20), labels = axisText) +
      scale_colour_manual(limits = c("control1_Seedling", "NLRs_Seedling"), 
                          values=c("grey43", "black"), labels = c("Controls", "R-genes")) +
      geom_line(aes(group = SampleGenes),linewidth = 1.3) +
      geom_point(aes(group = SampleGenes), size = 2) + theme_minimal() +
      labs(x = "", y = "Average proportion of gene region") +
      geom_vline(xintercept=0, color="grey", linewidth=1) +
      coord_cartesian(ylim= c(0,1), clip = "off") + theme(plot.margin = unit(c(1,1,2,1), "lines")) +
      annotation_custom(textGrob("% of gene length from TSS", gp=gpar(fontsize=14, col = "grey33")),xmin=0,xmax=100,ymin=-22,ymax=-22) + 
      annotation_custom(textGrob("Gene region", gp=gpar(fontsize=16)),xmin=0,xmax=100,ymin=-30,ymax=-30) +
      theme(axis.text.x = element_text(size = 13, colour = "black", angle = 45, vjust = 1, hjust = 1), axis.text.y = element_text(size = 14,colour = "black"), 
            axis.title.y = element_text(size = 16, vjust = 2)) 
    
    ggsave(paste(mod, "_", ".LeavesFrequencies.pdf", sep = ""), plot = plot, width = 12, height = 6) 
  }
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

# Determine whether the expression is more similar between genes within a cluster than between all genes.
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

# Create a hash for storing the proportion of each R-gene region modified with each mark.
sampleGenesProportions <- hash()

# Choose ecotype and tissue for analysis.
# Options: ColLeaf, ColRoot
tissueForAnalysis <- "ColLeaf"

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
    
    # Collect all hashes for modProportionPerRegion into single dataframes.
    modProportionPerRegion <- mergeResults(modProportionPerRegion)
    
    # Add a column to modProportionPerRegion with the numbers for 
    # each gene region that will correspond with their position on the x axis.
    modProportionPerRegion <- geneRegionAxisLocations(modProportionPerRegion, geneRegions)
    
    # Add a column to modProportionPerRegion with the current expression level.
    modProportionPerRegion <- expressionColumn(modProportionPerRegion, cluster)
    
    # Store final results on the appropriate hash.
    sampleGenesProportions[[test]][[cluster]] <- modProportionPerRegion
  }
}


# Merge all data from all sample gene sets into one big dataframe.
allResultsProportions <- data.frame()

for (test in names(sampleGenesProportions)) {
  for (cluster in names(sampleGenesProportions[[test]])) {
    
    df2 <- sampleGenesProportions[[test]][[cluster]]
    df2 <- cbind(df2, data.frame(SampleGenes = rep(test, times = nrow(df2))))
    
    allResultsProportions <- rbind(allResultsProportions, df2)
  }
}

# Calculate the difference in the proportion of each gene region covered by each modification between R-genes of the same cluster.
modificationDifference_Clusters <- data.frame()

for (mod in unique(allResultsProportions$Modification)) {
  df <- allResultsProportions[allResultsProportions$Modification == mod,]
  
  for (r in unique(df$Region)) {
    df1 <- df[df$Region == r,]
    
    for (cluster in unique(df1[grepl("cluster", df1$Expression),"Expression"])) {
      df2 <- df1[df1$Expression == cluster,]
      
      # For all pairwise comparisons,
      for (i in 1:nrow(df2)) {
        for (j in 1:nrow(df2)) {
          if (i != j & j > i) {
            # calculate the % difference in the proportion of each gene region covered by each modification,
            # using the abs() function to make all values positive.
            modificationDifference_Clusters <- rbind(modificationDifference_Clusters, 
                                                   data.frame(Region = r,
                                                              Modification = mod,
                                                              axisGroup = df2$axisGroup,
                                                              ExpressionDifference = (abs(df2[i, "Proportion"]-df2[j, "Proportion"])/df2[i, "Proportion"])*100,
                                                              Cluster = cluster))
          }
        }
      }
    }
  }
}

write.csv(modificationDifference_Clusters, file = "Cluster modification comparisons.csv")

# Calculate the difference in the proportion of each gene region covered by each modification between all R-genes.
modificationDifference_all <- data.frame()


for (mod in c("H3K4me3","H3K36me3","H3K9ac","H3K27ac","H3K27me1","H2AK121ub","H3K27me3","H3K9me2")) {
  df <- allResultsProportions[allResultsProportions$Modification == mod,]
  
  for (r in unique(df$Region)) {
    df1 <- df[df$Region == r,]
    
    for (i in 1:nrow(df1)) {
      for (j in 1:nrow(df1)) {
        if (i != j & j > i) {
      modificationDifference_all <- rbind(modificationDifference_all, 
                                          data.frame(Region = r,
                                                     Modification = mod,
                                                     axisGroup = df1$axisGroup[1],
                                                     ExpressionDifference = (abs(df1[i, "Proportion"]-df1[j, "Proportion"])/df1[i, "Proportion"])*100))
        }
      }
    }
    print(r)
    
  }
  print(mod)
}

write.csv(modificationDifference_all, file = "All R-gene modification comparisons.csv")


expressionDifference <- data.frame(Comparison = rep("Between clustered R-genes", times = length(modificationDifference_Clusters)),
																	 ExpressionDifference = modificationDifference_Clusters[,c(1:4)])

expressionDifference <- rbind(expressionDifference, data.frame(Comparison = rep("Between all R-genes", times = length(modificationDifference_all)),
																															 ExpressionDifference = modificationDifference_all))

colnames(expressionDifference) <- c("Comparison", "Region", "Modification", "axisGroup", "Difference")

statTestDF <- data.frame(Modification = character(),
                         Region = character(),
                         W.statistic = numeric(),
                         p.value = numeric())

for (mod in c("H3K4me3","H3K36me3","H3K9ac","H3K27ac","H3K27me1","H2AK121ub","H3K27me3","H3K9me2")) {
	df <- expressionDifference[expressionDifference$Modification==mod,]
	
	for (r in unique(df$Region)) {
	  df1 <- df[df$Region==r,]
	  statTest <- wilcox.test(Difference~Comparison, df1)
	  
	  statTestDF <- rbind(statTestDF, data.frame(Modification= mod,
	                                             Region= r,
	                                             W.statistic = statTest$statistic,
	                                             p.value = statTest$p.value))
	}

	plot <- ggplot(df, aes(x = axisGroup, y = Difference)) + 
		scale_x_continuous(limits = c(-70, 150), breaks = seq(-60, 140, 20), labels = axisText) +
		geom_boxplot(aes(group = Region)) + theme_minimal() + coord_cartesian(ylim= c(0,1), clip = "off") +
		labs(x = "", y = "Difference in proportion of gene region modified") +
		geom_vline(xintercept=0, color="grey", linewidth=1) + theme(plot.margin = unit(c(1,1,1,1), "lines")) +
		annotation_custom(textGrob("% of gene length from TSS", gp=gpar(fontsize=12, col = "grey33")),xmin=-10,xmax=100,ymin=-.23,ymax=-.23) + 
		annotation_custom(textGrob("Gene region", gp=gpar(fontsize=14)),xmin=-10,xmax=100,ymin=-.3,ymax=-.3) +
		theme(axis.text.x = element_text(size = 11, colour = "black", angle = 45, vjust = 1, hjust = 1), axis.text.y = element_text(size = 12,colour = "black"), 
					axis.title.y = element_text(size = 14, vjust = 2), strip.text = element_text(size = 16)) + 
	  facet_wrap(~Comparison) + 
	  stat_summary(fun="mean", geom="point", color="black", size=2) +
	  stat_summary(fun="mean", geom="line", color="black", size=1)
	
	
	ggsave(paste(mod, "_", "Clustered vs Unclustered.pdf", sep = ""), plot = plot, width = 16, height = 6)
}


write.csv(statTestDF, file = "Difference in proportion of gene region modified.csv")