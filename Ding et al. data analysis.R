source("Functions\\Overlaps functions.R")
source("Functions\\Modifications per gene.R")
source("Functions\\Coordinates per gene region.R")
source("Functions\\Modification frequencies & proportions.R")
source("Functions\\Get range - merge gene coordinates.R")


# Import the expression data from the ACRs paper (Ding et al. 2021).
Ding_ExpressionData <- as.data.frame(read_xlsx("Data\\ACRs data Ding et al., 2021.xlsx", sheet = 1))

# Filter for R-genes.
Ding_ExpressionData <- Ding_ExpressionData[which(Ding_ExpressionData$Gene %in% sampleGenes[["NLRs"]]$Gene),]

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

write.csv(Ding_ExpressionData, file = "Data\\Ding_ExpressionData.csv")

# Add information needed for the ReMap analysis.
Ding_ExpressionData <- cbind(Ding_ExpressionData, sampleGenes[["NLRs"]][,c(2:8)])

# Perform ReMap analysis to determine the chromatin modification enrichment in each R-gene.
# Determine which genes appear to be the most responsive to infection - are there any trend in their chromatin modifications.

# Import filtered ReMap2022 data.
ReMap <- as.data.frame(read_xlsx("Data\\Filtered ReMap data.xlsx"))

# Create list of chromatin modifications.
epiMods <- unique(ReMap$epiMod)

# Create a hash with the ReMap data in a particular tissue for the current set of genes. 
allModifications <- ReMapPerGene(Ding_ExpressionData, "seedlingGenes")

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


# Create a scatter plot to determine whether there is a correlation between expression level and the enrichment of each
# chromatin modification in each gene region.
for (n in 5:12) {
  df <- as.data.frame(read_xlsx("Data\\ACRs data Ding et al., 2021.xlsx", sheet = n))
  df <- df[-which(df$Control_Expression > 150),]
  
  df1 <- data.frame()
  
    for (row in 1:nrow(df)) {
      for (r in colnames(df)[c(4:13)]) {
        
      df1 <- rbind(df1, data.frame(Gene = df[row, "Gene"],
                                   Expression = df[row, "Control_Expression"],
                                   Region = r,
                                   Proportion = df[row, r]))
    }
  }
  plot <- ggplot(df1, aes(x = Expression, y = Proportion)) +
    geom_point() + geom_smooth() + facet_wrap(~Region)
}

# Create a scatter plot to determine whether there is a correlation between the change in expression level between control 
# and ETI conditions and the enrichment of each chromatin modification in each gene region - is there evidence of particular 
# chromatin modifications being associated with priming?
for (n in 5:12) {
  df <- as.data.frame(read_xlsx("Data\\ACRs data Ding et al., 2021.xlsx", sheet = n))
  
  df1 <- data.frame()
  
  for (row in 1:nrow(df)) {
    for (r in colnames(df)[c(4:13)]) {
      
      df1 <- rbind(df1, data.frame(Gene = df[row, "Gene"],
                                   Expression = df[row, "ETI_Expression"] - df[row, "Control_Expression"],
                                   Region = r,
                                   Proportion = df[row, r]))
    }
  }
  plot <- ggplot(df1, aes(x = Expression, y = Proportion)) +
    geom_point() + geom_smooth() + facet_wrap(~Region)
}


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

# Create bed files to visualise in IGV.
Ding_Control_ACR_Bed <- GRanges(
  seqnames=Rle(Ding_Control_ACR$seqnames),
  ranges=IRanges(Ding_Control_ACR$ranges),
  name=Ding_Control_ACR$Gene)

rtracklayer::export.bed(Ding_Control_ACR_Bed, "Data\\Ding_Control_ACR_Bed.bed")

Ding_ETI_ACR_Bed <- GRanges(
  seqnames=Rle(Ding_ETI_ACR$seqnames),
  ranges=IRanges(Ding_ETI_ACR$ranges),
  name=Ding_ETI_ACR$Gene)

rtracklayer::export.bed(Ding_ETI_ACR_Bed, "Data\\Ding_ETI_ACR_Bed.bed")

# Determine which ACRs overlap with the R-genes.
Ding_Control_ACR <- cbind(Ding_Control_ACR, data.frame(Condition = rep("Control", times = nrow(Ding_Control_ACR))))
Ding_ETI_ACR <- cbind(Ding_ETI_ACR, data.frame(Condition = rep("ETI", times = nrow(Ding_ETI_ACR))))

ACR_data <- Ding_Control_ACR
ACR_data <- rbind(ACR_data, Ding_ETI_ACR)

rm(Ding_Control_ACR, Ding_Control_ACR_Bed, Ding_ETI_ACR, Ding_ETI_ACR_Bed)

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

Leaf_NE_associations_Bed <- GRanges(
  seqnames=Rle(Leaf_NE_data$seqnames),
  ranges=IRanges(Leaf_NE_data$ranges),
  name=Leaf_NE_data$Gene)

rtracklayer::export.bed(Leaf_NE_associations_Bed, "Data\\Leaf_NE_associations_Bed.bed") 

Root_NE_associations_Bed <- GRanges(
  seqnames=Rle(Root_NE_data$seqnames),
  ranges=IRanges(Root_NE_data$ranges),
  name=Root_NE_data$Gene)

rtracklayer::export.bed(Root_NE_associations_Bed, "Data\\Root_NE_associations_Bed.bed") 


# Which TFs are the R-genes associated with - is there a correlation between the enrichment of particular chromatin 
# modifications and particular TFs?

# Import the TF data from the ACRs paper (Ding et al. 2021).
Ding_TFs <- as.data.frame(read_xlsx("Data\\ACRs data Ding et al., 2021.xlsx", sheet = 4))

# Filter for R-genes.
Ding_TFs <- Ding_TFs[which(Ding_TFs$target %in% sampleGenes[["NLRs"]]$Gene),]


# Add sheets to 'ACRs Ding et al., 2021' spreadsheet giving a summary of the enrichment of chromatin modifications
# and the presence of ACRs in each gene region.
wb <- loadWorkbook("Data\\ACRs data Ding et al., 2021.xlsx")

for (mod in c("H3K9me2","H3K27me3","H2A-Z","H2AK121ub","H3K4me3","H3K36me3","H3K27ac","H3K9ac")) {
  
  df <- modProportionPerRegion[grepl(mod, modProportionPerRegion$Modification),]

  Leaf_NEAs <- c()
  Root_NEAs <- c()
  
  df1 <- df[df$Region=="Downstream",]
  
  for (row in 1:nrow(df1)) {
   
    overlappingLeafNEAs <- Leaf_NE_data[which(Leaf_NE_data$Gene==df1[row,"Gene"] & Leaf_NE_data$Region==df1[row,"Region"]),]
    
    if (nrow(overlappingLeafNEAs) != 0) {
      Leaf_NEAs <- append(Leaf_NEAs, "Yes")
    } else Leaf_NEAs <- append(Leaf_NEAs, "No")
    
    overlappingRootNEAs <- Root_NE_data[which(Root_NE_data$Gene==df1[row,"Gene"] & Root_NE_data$Region==df1[row,"Region"]),]
    
    if (nrow(overlappingRootNEAs) != 0) {
      Root_NEAs <- append(Root_NEAs, "Yes")
    } else Root_NEAs <- append(Root_NEAs, "No")
  }
  
  associatedTFs <- data.frame()
  for (gene in unique(df1$Gene)) {
    if (gene %in% Ding_TFs$target) {
      associatedTFs <- rbind(associatedTFs, data.frame(TFs = paste(Ding_TFs[which(Ding_TFs$target==gene), "TF_alias"], collapse = ", ")))
    } else associatedTFs <- rbind(associatedTFs, data.frame(TFs = " "))
  }
  
  DingDataResults <- data.frame(Gene = Ding_ExpressionData$Gene,
                                Control_Expression = Ding_ExpressionData$Control,
                                ETI_Expression = Ding_ExpressionData$ETI,
                                Leaf_NE_association = Leaf_NEAs,
                                Root_NE_association = Root_NEAs)
  
  names <- c()
  for (r in unique(df$Region)) {
    df1 <- df[df$Region==r,]
    
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
    
    DingDataResults <- cbind(DingDataResults, df1$Proportion)
    DingDataResults <- cbind(DingDataResults, control_ACR)
    DingDataResults <- cbind(DingDataResults, ETI_ACR)
    names <- append(names, c(paste("Enrichment_", r, sep = ""), paste("Control_ACR_", r, sep = ""), paste("ETI_ACR_", r, sep = "")))
  }
  
  DingDataResults <- cbind(DingDataResults, data.frame(TFs = associatedTFs$TFs))
  
  colnames(DingDataResults)[6:35] <- names
  
  addWorksheet(wb,mod)
  writeData(wb,mod,DingDataResults)
  saveWorkbook(wb,"Data\\ACRs data Ding et al., 2021.xlsx",overwrite = TRUE)
}


# Are any of the R-gene ACRs associated with super-enhancers?
SE_ACRs <- as.data.frame(read_xlsx("Data\\Super-enhancer ACRs.xlsx"))

SE_ACR_overlaps <- data.frame()

for (i in 1:nrow(ACR_data)) {
  for (j in 1:nrow(SE_ACRs)) {
    if (ACR_data[i,"seqnames"]==SE_ACRs[j,"seqnames"] & overlapsFunction(ACR_data[i, "start"], ACR_data[i, "end"],
                         SE_ACRs[j,"start"], SE_ACRs[j, "end"])==TRUE) {
      
      SE_ACR_overlaps <- rbind(SE_ACR_overlaps, data.frame(seqnames = SE_ACRs[j, "seqnames"],
                                                           Gene = ACR_data[i, "Gene"],
                                                           Region = ACR_data[i, "Region"],
                                                           ACR_start = ACR_data[i, "start"],
                                                           ACR_end = ACR_data[i, "end"],
                                                           start = SE_ACRs[j, "start"],
                                                           end = SE_ACRs[j, "end"],
                                                           Tissues = SE_ACRs[j, "tissues"]))
    }
    else SE_ACR_overlaps <- SE_ACR_overlaps
  }
}

SE_ACR_overlaps$SE_ranges <- mergeCoordinates(SE_ACR_overlaps)

# One SE_ACR in the intergenic region between AT4G19520 and AT4G19530.
SE_ACR_Bed <- GRanges(
  seqnames=Rle(SE_ACR_overlaps$seqnames),
  ranges=IRanges(SE_ACR_overlaps$SE_ranges),
  name=SE_ACR_overlaps$Gene)

rtracklayer::export.bed(SE_ACR_Bed, "Data\\SE_ACR_Bed.bed") 


# Are there similarities in chromatin modification and TF enrichment between co-expressed genes?
coExpressionData <- as.data.frame(read_xlsx("Data\\Arabidopsis NLRs.xlsx"))
coExpressionData <- coExpressionData[,-c(5:8,15)]

coExpressionMatrix <- matrix(rep(0, times = nrow(coExpressionData)^2), nrow = nrow(coExpressionData), ncol = nrow(coExpressionData))
colnames(coExpressionMatrix) <- coExpressionData$Gene
rownames(coExpressionMatrix) <- coExpressionData$Gene

for (row in 1:nrow(coExpressionData)) {
  coExpressedGenes <- unlist(str_split(coExpressionData[row, "Co-expression"], ", "))
  
  if (length(coExpressedGenes) != 0 & NA %!in% coExpressedGenes) {
    for (gene in coExpressedGenes) {
      coExpressionMatrix[coExpressionData[row,"Gene"], gene] <- 1
      coExpressionMatrix[gene, coExpressionData[row,"Gene"]] <- 1
    }
  }
  else coExpressionMatrix <- coExpressionMatrix
}