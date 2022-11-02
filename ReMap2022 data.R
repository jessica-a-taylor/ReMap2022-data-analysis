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

# Import ReMap2022 data.
ReMap <- rtracklayer::import.bed("C:\\Users\\jexy2\\OneDrive\\Documents\\PhD\\remap2022_histone_all_macs2_TAIR10_v1_0.bed.gz")

# Convert to a dataframe and define column names.
ReMap <- as.data.frame(ReMap, colnames = c("seqnames", "start", "end", "width",
                                           "strand", "name", "score", "itemRgb",
                                           "thick.start", "thick.end", "thick.width"))

# Remove unwanted columns.
ReMap <- ReMap[,-c(9:11)]

# Import all Arabidopsis genes.
Atgenes <- as.data.frame(transcriptsBy(TxDb.Athaliana.BioMart.plantsmart28, by="gene"))

# Remove duplicate genes (different versions).
Atgenes <- Atgenes[-c(which(Atgenes$tx_name == str_match(Atgenes$tx_name, "^([0-9a-zA-Z]+)([.])([2-9])$")[,1])),]

pericentromericRegions <- data.frame(Chromosome = c(1:5),
                                     Start = c("11500000", "1100000", "10300000", "1500000", "9000000"),
                                     End = c("17700000", "7200000", "17300000", "6300000", "16000000"))

euchromaticRegions <- data.frame()

for (row in 1:nrow(pericentromericRegions)) {
  df <- Atgenes[c(which(Atgenes$seqnames==row & Atgenes$start < as.numeric(pericentromericRegions[row, "Start"]))),]
  df <- rbind(df, Atgenes[c(which(Atgenes$seqnames==row & Atgenes$end < as.numeric(pericentromericRegions[row, "End"]))),])
  
  euchromaticRegions <- rbind(euchromaticRegions, df)
}

rm(pericentromericRegions)

euchromaticRegions <- euchromaticRegions[,-c(1,8,9)]
colnames(euchromaticRegions)[1] <- "Gene"
euchromaticRegions$ranges <- paste(euchromaticRegions$start,"-",euchromaticRegions$end, sep = "")

# Remove duplicate genes.
newEuchromaticRegions <- euchromaticRegions
euchromaticRegions <- data.frame()

for (gene in unique(newEuchromaticRegions$Gene)) {
  euchromaticRegions <- rbind(euchromaticRegions, newEuchromaticRegions[newEuchromaticRegions$Gene==gene,][1,])
}

rm(newEuchromaticRegions)

# Remove TEs from the euchromaticRegions dataframe.
transposableElements <- as.data.frame(read_xlsx("C:\\Users\\jexy2\\OneDrive\\Documents\\PhD\\Arabidopsis TE genes.xlsx"))

withoutTEs <- euchromaticRegions[-c(which(euchromaticRegions$Gene %in% transposableElements$Locus)),]

rm(transposableElements)


# Create function that determines whether value a is between values b and c.
betweenFunction <- function(a,b,c) {
  return(b<a & a<c)
}

# Create a function that determines whether two ranges overlap using the between function.
overlapsFunction <- function(S1, E1, S2, E2) {
  if (betweenFunction(S1, S2, E2)) {
    return (TRUE)
  }
  if (betweenFunction(E1, S2, E2)) {
    return (TRUE)
  }
  if (betweenFunction(S2, S1, E1)) {
    return (TRUE)
  }
  if (betweenFunction(E2, S1, E1)) {
    return (TRUE)
  }
  return(FALSE)
}

# Function to return the index of the set containing "item" in "overlapSets"
findItem <- function(item, overlapSets) {
  
  # For each set in overlapSets
  for (setIndex in 1:length(overlapSets)) {
    
    # If item is contained within that set, return the index of that set
    if (item %in% overlapSets[[setIndex]]) {
      return(setIndex)
    }
  }
}



# Get 10 sets of random genes and store in a hash, using either euchromaticRegions or withoutTEs.
dataToUse <- withoutTEs
  
testData <- hash()

for (n in c(1:10)) {
  testData[[paste("control", n, sep = "")]] <- dataToUse[c(sample(nrow(dataToUse), 200)),]
}


# Get the coordinates for the gene bodies of each NLR.
NLRgenes <- as.data.frame(read_xlsx("C:\\Users\\jexy2\\OneDrive\\Documents\\PhD\\Arabidopsis NLRs.xlsx", sheet = 1))

NLRgenebody <- Atgenes[which(Atgenes$group_name %in% NLRgenes$Gene),]

# Create a ranges column by merging the start and end columns.
NLRgenebody$ranges <- paste(NLRgenebody$start,"-",NLRgenebody$end, sep = "")
colnames(NLRgenebody)[2] <- "Gene"

genebodyBed <- GRanges(
  seqnames=Rle(NLRgenebody$seqnames),
  ranges=IRanges(NLRgenebody$ranges),
  name=NLRgenebody$Gene)

#rtracklayer::export.bed(genebodyBed, "NLRgenebody.bed")


# Add NLR genes to testData.
testData[["NLRs"]] <- NLRgenebody



# Perform all analyses on the testData.
for (test in names(testData)) {
  dataToUse <- testData[[test]]
  
  # Create new dataframes for chunks of the gene body (20% intervals of the gene length).
  geneChunks <- hash(width20 = dataToUse[,c(which(colnames(dataToUse) != "start" & colnames(dataToUse) != "end" &
                                                    colnames(dataToUse) != "width" & colnames(dataToUse) != "ranges"))], 
                     width40 = dataToUse[,c(which(colnames(dataToUse) != "start" & colnames(dataToUse) != "end" &
                                                    colnames(dataToUse) != "width" & colnames(dataToUse) != "ranges"))], 
                     width60 = dataToUse[,c(which(colnames(dataToUse) != "start" & colnames(dataToUse) != "end" &
                                                    colnames(dataToUse) != "width" & colnames(dataToUse) != "ranges"))], 
                     width80 = dataToUse[,c(which(colnames(dataToUse) != "start" & colnames(dataToUse) != "end" &
                                                    colnames(dataToUse) != "width" & colnames(dataToUse) != "ranges"))], 
                     width100 = dataToUse[,c(which(colnames(dataToUse) != "start" & colnames(dataToUse) != "end" &
                                                     colnames(dataToUse) != "width" & colnames(dataToUse) != "ranges"))])
  
  geneWidth <- hash(width20 = c(), 
                    width40 = c(), 
                    width60 = c(), 
                    width80 = c(), 
                    width100 = c())
  
  for (row in 1:nrow(dataToUse)) {
    if (dataToUse[row, "strand"]=="+") {
      geneWidth[["width20"]] <- append(geneWidth[["width20"]], paste(dataToUse[row,"start"],"-", dataToUse[row,"start"] + dataToUse[row,"width"]*0.2, sep = ""))
      geneWidth[["width40"]] <- append(geneWidth[["width40"]], paste(dataToUse[row,"start"] + dataToUse[row,"width"]*0.2+1,"-", dataToUse[row,"start"] + dataToUse[row,"width"]*0.4, sep = ""))
      geneWidth[["width60"]] <- append(geneWidth[["width60"]], paste(dataToUse[row,"start"] + dataToUse[row,"width"]*0.4+1,"-", dataToUse[row,"start"] + dataToUse[row,"width"]*0.6, sep = ""))
      geneWidth[["width80"]] <- append(geneWidth[["width80"]], paste(dataToUse[row,"start"] + dataToUse[row,"width"]*0.6+1,"-", dataToUse[row,"start"] + dataToUse[row,"width"]*0.8, sep = ""))
      geneWidth[["width100"]] <- append(geneWidth[["width100"]], paste(dataToUse[row,"start"] + dataToUse[row,"width"]*0.8+1,"-",dataToUse[row,"end"], sep = ""))
    }
    else if (dataToUse[row, "strand"]=="-"){
      geneWidth[["width20"]] <- append(geneWidth[["width20"]], paste(dataToUse[row,"end"] - dataToUse[row,"width"]*0.2,"-", dataToUse[row,"end"], sep = ""))
      geneWidth[["width40"]] <- append(geneWidth[["width40"]], paste(dataToUse[row,"end"] - dataToUse[row,"width"]*0.4,"-", dataToUse[row,"end"] - dataToUse[row,"width"]*0.2-1, sep = ""))
      geneWidth[["width60"]] <- append(geneWidth[["width60"]], paste(dataToUse[row,"end"] - dataToUse[row,"width"]*0.6,"-", dataToUse[row,"end"] - dataToUse[row,"width"]*0.4-1, sep = ""))
      geneWidth[["width80"]] <- append(geneWidth[["width80"]], paste(dataToUse[row,"end"] - dataToUse[row,"width"]*0.8,"-", dataToUse[row,"end"] - dataToUse[row,"width"]*0.6-1, sep = ""))
      geneWidth[["width100"]] <- append(geneWidth[["width100"]], paste(dataToUse[row,"start"],"-",dataToUse[row,"end"] - dataToUse[row,"width"]*0.8-1, sep = ""))
    }
    }
  
  for (n in names(geneChunks)) {
    geneChunks[[n]]$ranges <- geneWidth[[n]]
  }
  
  rm(geneWidth)
  
  for (n in names(geneChunks)) {
    start <- c()
    end <- c()
    
    for (row in 1:nrow(geneChunks[[n]])) {
      start <- append(start, str_match(geneChunks[[n]][row,"ranges"], "^(\\d*\\.?\\d+)(-)(\\d*\\.?\\d+)$")[,2])
      end <- append(end, str_match(geneChunks[[n]][row,"ranges"], "^(\\d*\\.?\\d+)(-)(\\d*\\.?\\d+)$")[,4])
    }
    geneChunks[[n]]$start <- start
    geneChunks[[n]]$end <- end
  }
  
  for (n in names(geneChunks)) {
    geneBed <- GRanges(
      seqnames=Rle(geneChunks[[n]]$seqnames),
      ranges=IRanges(geneChunks[[n]]$ranges),
      name=geneChunks[[n]]$Gene)
    
    #rtracklayer::export.bed(geneBed, paste("NLR", n, ".bed", sep = ""))
  }
  
  # Create new dataframe for the coordinates of the regions 200bp downstream of the TTS.
  
  downstreamRegion <- c()
  for (row in 1:nrow(dataToUse)) {
    if (dataToUse[row, "strand"]=="+") {
      downstreamRegion <- append(downstreamRegion, paste(dataToUse[row,"end"],"-",dataToUse[row,"end"]+200, sep = ""))
    }
    else if (dataToUse[row, "strand"]=="-") {
      downstreamRegion <- append(downstreamRegion, paste(dataToUse[row,"start"]-200,"-", dataToUse[row,"start"], sep = ""))
    }
  }
  
  downstream <- dataToUse[,c(which(colnames(dataToUse) != "start" & colnames(dataToUse) != "end" &
                                     colnames(dataToUse) != "width" & colnames(dataToUse) != "ranges"))]
  
  downstream$ranges <- downstreamRegion
  rm(downstreamRegion)
  
  downstream$start <- str_match(downstream$ranges, "^([0-9]+)(-)([0-9]+)$")[,2]
  downstream$end <- str_match(downstream$ranges, "^([0-9]+)(-)([0-9]+)$")[,4]
  
  downstreamBed <- GRanges(
    seqnames=Rle(downstream$seqnames),
    ranges=IRanges(downstream$ranges),
    name=downstream$Gene)
  
  #rtracklayer::export.bed(downstreamBed, "downstream.bed")
  
  # Get the coordinates for the promotors of each gene.
  ATpromotors500 <- promoters(TxDb.Athaliana.BioMart.plantsmart28, upstream=500, downstream=0, use.names = TRUE)
  ATpromotors1000 <- promoters(TxDb.Athaliana.BioMart.plantsmart28, upstream=1000, downstream=0, use.names = TRUE)
  
  # Remove duplicate genes (different versions).
  ATpromotors500 <- ATpromotors500[-c(which(ATpromotors500$tx_name == str_match(ATpromotors500$tx_name, "^([0-9a-zA-Z]+)([.])([2-9])$")[,1])),]
  ATpromotors1000 <- ATpromotors1000[-c(which(ATpromotors1000$tx_name == str_match(ATpromotors1000$tx_name, "^([0-9a-zA-Z]+)([.])([2-9])$")[,1])),]
  
  promotor500 <- data.frame(seqnames = numeric(),
                            start = numeric(),
                            end = numeric(),
                            width = numeric(),
                            strand = factor(),
                            tx_id = numeric(),
                            tx_name = character())
  
  promotor1000 <- promotor500
  
  for (gene in dataToUse$Gene) {
    promotor500 <- rbind(promotor500, as.data.frame(ATpromotors500[grepl(gene,ATpromotors500$tx_name),]))
    promotor1000 <- rbind(promotor1000, as.data.frame(ATpromotors1000[grepl(gene,ATpromotors1000$tx_name),]))
  }
  
  # Alter coordinated of promotor1000 to be only 500bp upstream of promototr500.
  for (row in 1:nrow(promotor1000)) {
    if (promotor1000[row, "strand"]=="+") {
      promotor1000[row, "end"] <- promotor1000[row, "start"]+500
    }
    else if (promotor1000[row, "strand"]=="-") {
      promotor1000[row, "start"] <- promotor1000[row, "end"]-500
    }
  }
  
  # Create a ranges column by merging the start and end columns.
  promotor500$ranges <- paste(promotor500$start,"-",promotor500$end, sep = "")
  promotor1000$ranges <- paste(promotor1000$start,"-",promotor1000$end, sep = "")
  
  
  # Add a new column for the gene name, removing ".1" from the end.
  promotor500$group_name <- str_match(promotor500$tx_name, "^([0-9a-zA-Z]+)([.])([1])$")[,2]
  
  promotor500Bed <- GRanges(
    seqnames=Rle(promotor500$seqnames),
    ranges=IRanges(promotor500$ranges),
    name=promotor500$tx_name)
  
  #rtracklayer::export.bed(promotor500Bed, "promotor500.bed")
  
  promotor1000$group_name <- str_match(promotor1000$tx_name, "^([0-9a-zA-Z]+)([.])([1])$")[,2]
  
  promotor1000Bed <- GRanges(
    seqnames=Rle(promotor1000$seqnames),
    ranges=IRanges(promotor1000$ranges),
    name=promotor1000$tx_name)
  
  #rtracklayer::export.bed(promotor1000Bed, "promotor1000.bed")
  
  rm(ATpromotors500, ATpromotors1000)
  
  # Get the coordinates for the upstream intergenic regions of each gene.
  usCoordinates <- c()
  
  for (gene in dataToUse$Gene) {
    currentGene <- which(Atgenes$group_name==gene)
    
    if (Atgenes[currentGene, "strand"]=="+") {
      previousGene <- currentGene - 1
      
      if (previousGene > 0 & as.numeric(Atgenes[currentGene, "seqnames"])==as.numeric(Atgenes[previousGene, "seqnames"])) {
        if (previousGene > 0 & Atgenes[previousGene, "strand"]=="+") {
          distance <- (Atgenes[currentGene, "start"] - 1001) - (Atgenes[previousGene, "end"] + 201)
          
          if (distance > 0) {
            usCoordinates <- append(usCoordinates, paste(Atgenes[previousGene, "end"] + 201, "-", Atgenes[previousGene, "end"] + 201 + distance, sep = "")) 
          } 
          else usCoordinates <- append(usCoordinates, NA)
        }
        
        else if (previousGene > 0 & Atgenes[previousGene, "strand"]=="-") {
          distance <- (Atgenes[currentGene, "start"] - 1001) - (Atgenes[previousGene, "end"] + 1001)
          
          if (distance > 0) {
            usCoordinates <- append(usCoordinates, paste(Atgenes[previousGene, "end"] + 1001, "-", Atgenes[previousGene, "end"] + 1001 + distance, sep = "")) 
          } 
          else usCoordinates <- append(usCoordinates, NA)
          
        } 
      }
      else usCoordinates <- append(usCoordinates, NA)
    }
    
    else if (Atgenes[currentGene, "strand"]=="-") {
      previousGene <- currentGene + 1
      
      if (previousGene > 0 & as.numeric(Atgenes[currentGene, "seqnames"])==as.numeric(Atgenes[previousGene, "seqnames"])) {
        if (previousGene > 0 & Atgenes[previousGene, "strand"]=="+") {
          distance <- (Atgenes[previousGene, "start"] - 1001) - (Atgenes[currentGene, "end"] + 1001)
          
          if (distance > 0) {
            usCoordinates <- append(usCoordinates, paste(Atgenes[previousGene, "start"] - 1001 - distance, "-", Atgenes[previousGene, "start"] - 1001, sep = "")) 
          } 
          else usCoordinates <- append(usCoordinates, NA)
          
        }
        
        else if (previousGene > 0 & Atgenes[previousGene, "strand"]=="-") {
          distance <- (Atgenes[previousGene, "start"] - 201) - (Atgenes[currentGene, "end"] + 1001)
          if (distance > 0) {
            usCoordinates <- append(usCoordinates, paste(Atgenes[previousGene, "start"] - 201 - distance, "-", Atgenes[previousGene, "start"] - 201, sep = "")) 
          } 
          else usCoordinates <- append(usCoordinates, NA)
          
        }
      }
      else usCoordinates <- append(usCoordinates, NA)
    } 
  }
  
  upstreamIntergenic <- Atgenes[which(Atgenes$group_name %in% dataToUse$Gene),]
  upstreamIntergenic$ranges <- usCoordinates 
  
  upstreamIntergenicBed <- upstreamIntergenic[-c(which(is.na(upstreamIntergenic$ranges))),]
  
  upstreamIntergenicBed <- GRanges(
    seqnames=Rle(as.numeric(upstreamIntergenicBed$seqnames)),
    ranges=IRanges(upstreamIntergenicBed$ranges),
    name=upstreamIntergenicBed$group_name)
  
  #rtracklayer::export.bed(upstreamIntergenicBed, "upstreamIntergenic.bed")
  
  
  # Get the coordinates for the downstream intergenic regions of each NLR.
  dsCoordinates <- c()
  
  for (gene in dataToUse$Gene) {
    currentGene <- which(Atgenes$group_name==gene)
    
    if (Atgenes[currentGene, "strand"]=="-") {
      nextGene <- currentGene - 1
      
      if (as.numeric(Atgenes[currentGene, "seqnames"])==as.numeric(Atgenes[nextGene, "seqnames"])) {
        if (Atgenes[nextGene, "strand"]=="-") {
          distance <- (Atgenes[currentGene, "start"] - 201) - (Atgenes[nextGene, "end"] + 1001)
          
          if (distance > 0) {
            dsCoordinates <- append(dsCoordinates, paste(Atgenes[nextGene, "end"] + 1001, "-", Atgenes[nextGene, "end"] + 1001 + distance, sep = "")) 
          } 
          else dsCoordinates <- append(dsCoordinates, NA) 
        }
        
        else if (Atgenes[nextGene, "strand"]=="+") {
          distance <- (Atgenes[currentGene, "start"] - 201) - (Atgenes[nextGene, "end"] + 201)
          
          if (distance > 0) {
            dsCoordinates <- append(dsCoordinates, paste(Atgenes[nextGene, "end"] + 201, "-", Atgenes[nextGene, "end"] + 201 + distance, sep = "")) 
          } 
          else dsCoordinates <- append(dsCoordinates, NA)
          
        }
      }
      else dsCoordinates <- append(dsCoordinates, NA)
    }
    
    else if (Atgenes[currentGene, "strand"]=="+") {
      nextGene <- currentGene + 1
      
      if (as.numeric(Atgenes[currentGene, "seqnames"])==as.numeric(Atgenes[nextGene, "seqnames"])) {
        if (Atgenes[nextGene, "strand"]=="+") {
          distance <- (Atgenes[nextGene, "start"] - 1001) - (Atgenes[currentGene, "end"] + 201)
          
          if (distance > 0) {
            dsCoordinates <- append(dsCoordinates, paste(Atgenes[nextGene, "start"] - 1001 - distance, "-", Atgenes[nextGene, "start"] - 1001, sep = "")) 
          } 
          else dsCoordinates <- append(dsCoordinates, NA) 
          
        }
        
        else if (Atgenes[nextGene, "strand"]=="-") {
          distance <- (Atgenes[nextGene, "start"] - 201) - (Atgenes[currentGene, "end"] + 201)
          if (distance > 0) {
            dsCoordinates <- append(dsCoordinates, paste(Atgenes[nextGene, "start"] - 201 - distance, "-", Atgenes[nextGene, "start"] - 201, sep = "")) 
          } 
          else dsCoordinates <- append(dsCoordinates, NA) 
          
        }
      }
      else dsCoordinates <- append(dsCoordinates, NA)
    } 
  }
  
  downstreamIntergenic <- Atgenes[which(Atgenes$group_name %in% dataToUse$Gene),]
  downstreamIntergenic$ranges <- dsCoordinates 
  
  downstreamIntergenicBed <- downstreamIntergenic[-c(which(is.na(downstreamIntergenic$ranges))),]
  
  downstreamIntergenicBed <- GRanges(
    seqnames=Rle(as.numeric(downstreamIntergenicBed$seqnames)),
    ranges=IRanges(downstreamIntergenicBed$ranges),
    name=downstreamIntergenicBed$group_name)
  
  #rtracklayer::export.bed(downstreamIntergenicBed, "downstreamIntergenic.bed")
  
  
  # Create a hash to store ReMap data for each NLR/cluster.
  colnames(NLRgenebody)[2] <- "Gene"
  
  gene_hash <- hash()
  
  for (row in 1:nrow(dataToUse)) {
    # Select rows that are within the range of each gene and on the same chromosome.
    ReMapRows <- c(which(ReMap[,"start"] > dataToUse[row, "start"]-5000 & ReMap[,"end"] < dataToUse[row, "end"]+5000 & ReMap[,"seqnames"] == as.numeric(dataToUse[row, "seqnames"])))
    gene_hash[[dataToUse[row,"Gene"]]] <- ReMap[ReMapRows,]
  }
  
  rm(ReMapRows, dataToUse)
  
  for(n in names(gene_hash)) {
    if (nrow(gene_hash[[n]])>=1) {
      # Run regex on name column, extracting each section
      # (experiment, epigenetic modification, ecotype, other info)
      gene_hash[[n]][c("exp.", "epiMod", "ecotype", "info")] <- str_match(gene_hash[[n]][,"name"], "^([0-9a-zA-Z]+)\\.([0-9a-zA-Z-]+)\\.([0-9a-zA-Z-]+)[_\\.](.*)$")[,-1]
      
      # Filter epiMod column, excluding unwanted modifications
      gene_hash[[n]] <- gene_hash[[n]][!gene_hash[[n]]$epiMod %in% c("H3", "HTR12", "H2A", "H2B", "H3T3ph", "H1", "H2A-X",
                                                                     "H2AV", "HTA6", "H3-1") & 
                                         !gene_hash[[n]]$ecotype %in% c("C24", "undef", "Col-x-Ler", "Ler-x-Col", "Col-x-C24"),]
      
      gene_hash[[n]] <- gene_hash[[n]][,-6]
      
      
      # Filter info column, excluding unwanted conditions (too old, too young, wrong part of plant, etc)
      gene_hash[[n]] <-
        gene_hash[[n]][!grepl("mutant", gene_hash[[n]]$info) &
                         !grepl("mature", gene_hash[[n]]$info) &
                         !grepl("senescent", gene_hash[[n]]$info) &
                         !grepl("inflorescence", gene_hash[[n]]$info) &
                         !grepl("drought", gene_hash[[n]]$info) &
                         !grepl("old", gene_hash[[n]]$info) &
                         !grepl("min", gene_hash[[n]]$info) &
                         !grepl("endosperm", gene_hash[[n]]$info) &
                         !grepl("-se-", gene_hash[[n]]$info) &
                         !grepl("-TSA-", gene_hash[[n]]$info) &
                         !grepl("-GSNO-", gene_hash[[n]]$info) &
                         !grepl("flg22", gene_hash[[n]]$info) &
                         !grepl("transgenic", gene_hash[[n]]$info) &
                         !grepl("GSH", gene_hash[[n]]$info) &
                         !grepl("-acc1", gene_hash[[n]]$info) &
                         !grepl("-ethylene", gene_hash[[n]]$info) &
                         !grepl("-C2H4", gene_hash[[n]]$info) &
                         !grepl("leaves_3w-K36M-homoz", gene_hash[[n]]$info) &
                         !grepl("undef_seedling_10d-h3-1kd-1", gene_hash[[n]]$info) &
                         !grepl("-air", gene_hash[[n]]$info) &
                         !grepl("-ehylene", gene_hash[[n]]$info) &
                         !grepl("-swap", gene_hash[[n]]$info) &
                         !grepl("-K36M", gene_hash[[n]]$info) &
                         !grepl("-H3-KD", gene_hash[[n]]$info) &
                         !grepl("-water", gene_hash[[n]]$info) &
                         !grepl("undef_seedling_10d-h3-1kd-2", gene_hash[[n]]$info) &
                         !grepl("seedling_3d-wt-ehylene", gene_hash[[n]]$info) &
                         !grepl("GSE67322", gene_hash[[n]]$exp.) &
                         !grepl("GSE42695", gene_hash[[n]]$exp.) &
                         !grepl("GSE75071", gene_hash[[n]]$exp.) &
                         !grepl("GSE62615", gene_hash[[n]]$exp.) &
                         !grepl("GSE103361", gene_hash[[n]]$exp.) &
                         !grepl("GSE50636", gene_hash[[n]]$exp.) &
                         !grepl("GSE93223", gene_hash[[n]]$exp.) &
                         !grepl("GSE37644", gene_hash[[n]]$exp.) &
                         !grepl("GSE108414", gene_hash[[n]]$exp.) &
                         !grepl("GSE22276", gene_hash[[n]]$exp.) &
                         !grepl("GSE89768", gene_hash[[n]]$exp.) &
                         !grepl("GSE117391", gene_hash[[n]]$exp.) &
                         !grepl("brm", gene_hash[[n]]$info) &
                         !grepl("lhp1", gene_hash[[n]]$info) &
                         !grepl("atbmi", gene_hash[[n]]$info) &
                         !grepl("swn", gene_hash[[n]]$info) &
                         !grepl("ref6", gene_hash[[n]]$info) &
                         !grepl("arp6", gene_hash[[n]]$info) &
                         !grepl("clf", gene_hash[[n]]$info) &
                         !grepl("caa39", gene_hash[[n]]$info) &
                         !grepl("sdg8", gene_hash[[n]]$info) &
                         !grepl("atxr", gene_hash[[n]]$info) &
                         !grepl("hag1", gene_hash[[n]]$info) &
                         !grepl("OTU5", gene_hash[[n]]$info), ]
      
      # Filter info column, checking written plant ages and removing BAD AGES
      
      # For each row
      for (row in nrow(gene_hash[[n]]):1) {
        # Find the age if it exists (in the format 1w / 8d / 10h)
        matches <- str_match(gene_hash[[n]][row, "info"], "_([0-9]+)([dwh])")
        
        # matches will be of format ["_30h", "30", "h"] (or [NA, NA, NA])
        
        # If we found an age (ie, matches[1] is not NA)
        if (!is.na(matches[1])) {
          # Convert string number into integer
          timeValue <- as.numeric(matches[2])
          
          # Maximum allowed age in weeks
          maxWeeks <- 3
          
          badAge <- FALSE
          
          if (matches[3]=="h") {
            # If the age is measured in hours, it's too young. BAD AGE.
            badAge <- TRUE
          }
          
          else if (matches[3]=="d" & timeValue > maxWeeks*7) {
            # If the age is measured in days and it's longer than maxWeeks (converting maxWeeks to days), it's too old. BAD AGE.
            badAge <- TRUE
          }
          
          else if (matches[3]=="w" & timeValue > maxWeeks) {
            # If the age is measured in weeks and it's longer than maxWeeks, it's too old. BAD AGE.
            badAge <- TRUE
          }
          
          # If we had a BAD AGE, delete the corresponding row. (Otherwise, move on to the next row without deleting.)
          if (badAge) {
            gene_hash[[n]] <- gene_hash[[n]][-row,]
          }
        }
      }
      
      # Tidy up 🧹
      rm(matches, badAge, maxWeeks, row, timeValue)
    }
  }
  
  
  # Create hash for root data.
  RootNLRs <- hash()
  LeafNLRs <- hash()
  
  
  for(n in names(gene_hash)) {
    # Extract root data from gene_hash and store in RootNLRs.
    RootNLRs[[n]] <- gene_hash[[n]][grepl("roots",gene_hash[[n]]$info),]
    # Remove root data from gene_hash.
    LeafNLRs[[n]] <- gene_hash[[n]][!grepl("roots",gene_hash[[n]]$info),]
  }
  
  # Delete gene_hash.
  rm(gene_hash)
  
  
  # Create hashes for leaf data in Col and Ler ecotypes.
  ColLeafNLRs <- hash()
  LerLeafNLRs <- hash()
  
  for (n in names(LeafNLRs)) {
    ColLeafNLRs[[n]] <- LeafNLRs[[n]][grepl("Col-0",LeafNLRs[[n]]$ecotype),]
    LerLeafNLRs[[n]] <- LeafNLRs[[n]][grepl("Ler",LeafNLRs[[n]]$ecotype),]
  }
  
  rm(LeafNLRs)
  
  
  # Select tissue type for analysis: ColLeafNLRs, LerLeafNLRs or RootNLRs.
  dataToUse <- ColLeafNLRs

  # Create list of chromatin modifications.
  epiMods <- c()
  for (n in names(dataToUse)) {
    epiMods <- append(epiMods, unique(dataToUse[[n]]$epiMod))
  }
  
  epiMods <- unique(epiMods)
  
  # Create a dictionary (hash containing hashes) with dataframes for each epiMod for each NLR.
  ColWTdata <- hash()
  
  for (n in names(dataToUse)) {
    modHash <- hash()
    
    for (mod in epiMods) {
      modHash[[mod]] <- dataToUse[[n]][dataToUse[[n]]$epiMod==mod,]
    }
    ColWTdata[[n]] <- modHash
  }
  
  if (test = "NLRs") {
    forIGV <- ColWTdata
  }
  else next
  
  rm(dataToUse, modHash)
  
  # Merge start and end coordinates columns to create a ranges column.
  for (n in names(ColWTdata)) {
    for (mod in epiMods) {
      if (nrow(ColWTdata[[n]][[mod]]) >= 1) {
        ColWTdata[[n]][[mod]]$ranges <- paste(ColWTdata[[n]][[mod]]$start,"-",ColWTdata[[n]][[mod]]$end, sep = "")
      }
      else next
    }
  }
  
  
  # Create a hash with the data on the coordinates of each gene region.
  # First rename the columns in each dataset so they can be indexed the same way.
  colnames(promotor500)[9] ="Gene"
  colnames(promotor1000)[9] ="Gene"
  colnames(downstream)[2] ="Gene"
  colnames(upstreamIntergenic)[2] ="Gene"
  colnames(downstreamIntergenic)[2] ="Gene"
  
  for (n in names(geneChunks)) {
    colnames(geneChunks[[n]])[2]="Gene"
  }
  
  
  regions <- hash(UpstreamIntergenic = upstreamIntergenic, Promotor1000 = promotor1000, Promotor500 = promotor500,
                  Gene20 = geneChunks[["width20"]],  Gene40 = geneChunks[["width40"]], Gene60 = geneChunks[["width60"]], 
                  Gene80 = geneChunks[["width80"]], Gene100 = geneChunks[["width100"]], Downstream = downstream,
                  DownstreamIntergenic = downstreamIntergenic)
  
  # Create a dictionary containing the frequency of each chromatin modification occurring in each region of each NLR.
  modsPerRegion <- hash()
  
  for (r in names(regions)) {
    # Create a hash containing a list of chromatin modifications overlapping with the gene region.
    GBmod <- hash()
    
    for (n in names(ColWTdata)) {
      modList <- c()
      
      for (mod in epiMods) {
        modPresent <- FALSE
        
        if (nrow(ColWTdata[[n]][[mod]]) >= 1 & !is.na(regions[[r]][regions[[r]]$Gene==n,]$ranges)) {
          for (row in 1:nrow(ColWTdata[[n]][[mod]])) {
            if (overlapsFunction(ColWTdata[[n]][[mod]][row, "start"], ColWTdata[[n]][[mod]][row, "end"],
                                 regions[[r]][regions[[r]]$Gene==n,]$start, regions[[r]][regions[[r]]$Gene==n,]$end)==TRUE) {
              modPresent <- TRUE
            }
            else modPresent <- modPresent
          }
          if (modPresent == TRUE) {
            modList <- append(modList, mod)
          }
          else modPresent <- modPresent
        }
        else next
      }
      GBmod[[n]] <- modList
    }
    
    rm(modList, modPresent)
    
    # Calculate the percentage of NLRs with each chromatin modification within the gene body.
    modFrequencies <- hash()
    
    for (mod in epiMods) {
      modPresent <- 0
      
      for (n in names(GBmod)) {
        if (mod %in% GBmod[[n]]) {
          modPresent <- modPresent+1
        }
        else modPresent <- modPresent
      }
      modFrequencies[[mod]] <- modPresent/length(names(GBmod))*100
    }
    
    rm(modPresent)
    
    # Create a dataframe containing the percentage of NLRs with each chromatin modification within the gene body.
    modFrequenciesDF <- data.frame(Region = character(),
                                   Modification = character(),
                                   Frequency = numeric())
    
    for (n in names(modFrequencies)) {
      df <- data.frame(Region = rep(r, times = length(modFrequencies[[n]])),
                       Modification = n,
                       Frequency = modFrequencies[[n]])
      
      modFrequenciesDF <- rbind(modFrequenciesDF, df)
    }
    
    rm(df, modFrequencies)
    
    modsPerRegion[[r]] <- modFrequenciesDF
  }
  
  rm(modFrequenciesDF, regions)
  
  
  # Make a big dataframe with the modification frequencies for each gene region.
  level = c("UpstreamIntergenic", "Promotor1000", "Promotor500",
            "Gene20", "Gene40", "Gene60", "Gene80", "Gene100", 
            "Downstream", "DownstreamIntergenic")
  
  modFrequenciesDF <- data.frame(Region = character(),
                                 Modification = character(),
                                 Frequency = numeric())
  
  for (r in level) {
    modFrequenciesDF <- rbind(modFrequenciesDF, modsPerRegion[[r]])
  }
  
  # Add a column to the dataframe with the numbers on the x axis that will correspond with each gene region.
  grouping <- c(seq(from = -60, to = -20, by = 20),seq(from = 20, to = 140, by = 20))
  regions <- unique(modFrequenciesDF$Region)
  
  axisGroup <- c()
  for (c in 1:length(regions)) {
    axisGroup <- append(axisGroup, rep(grouping[c], times = nrow(modFrequenciesDF[modFrequenciesDF$Region==regions[c],])))
  }
  
  modFrequenciesDF <- cbind(modFrequenciesDF, axisGroup)
  
  rm(grouping, regions, axisGroup)
  
  testData[[test]] <- modFrequenciesDF
}
  
rm(ReMap)

# Merge all data into one big dataframe.
allResults <- data.frame()

for (test in names(testData)) {
  df <- testData[[test]]
  df <- cbind(df, data.frame(Test = rep(test, times = nrow(df))))
  
  allResults <- rbind(allResults, df)
}

rm(modFrequenciesDF)

# Plot the percentage of NLRs with each chromatin modification within the gene body.
axisText <- c("Intergenic", "Promotor \n(1kb)", "Promotor \n(500bp)", "TSS",
              "20%", "40%", "60%", "80%", "100%", 
              "Downstream \n(200bp)", "Intergenic")

for (mod in epiMods) {
  df <- allResults[allResults$Modification==mod,]

  modFrequenciesPlot <- ggplot(df, aes(x = axisGroup, y = Frequency, color = Test)) + 
    scale_x_continuous(limits = c(-60, 140), breaks = seq(-60, 140, 20), labels = axisText) +
    geom_line(aes(group = Test),size = 1.3) +
    geom_point(aes(group = Test), size = 2) + theme_minimal() + 
    scale_colour_manual(limits = c("control1", "NLRs"), values=c("grey43", "black"), labels = c("Controls", "R-genes")) +
    labs(x = "", y = "Frequency of occurrence (%)") +
    geom_vline(xintercept=0, color="grey", size=1) +
    coord_cartesian(ylim= c(0,100), clip = "off") + theme(plot.margin = unit(c(1,1,2,1), "lines")) +
    annotation_custom(textGrob("% of gene length from TSS", gp=gpar(fontsize=12, col = "grey33")),xmin=0,xmax=100,ymin=-14,ymax=-14) + 
    annotation_custom(textGrob("Gene region", gp=gpar(fontsize=14)),xmin=0,xmax=100,ymin=-20,ymax=-20) +
    theme(axis.text.x = element_text(size = 11, colour = "black"), axis.text.y = element_text(size = 12,colour = "black"), 
          axis.title.y = element_text(size = 14, vjust = 2)) 
  
  ggsave(paste(mod, ".plotLeaves.pdf", sep = ""), plot = modFrequenciesPlot, width = 12, height = 6)
}
  



# Create a bed file with chromatin modifications in each NLR to view in IGV.

allOverlaps <- hash()

# For each epigenetic modification name
for (n in names(forIGV)) {
  modOverlaps <- hash()
  
  for (mod in epiMods) {
    
    # Generate overlapSets as a list of single-item sets
    # eg, [ {1}, {2}, {3}, {4}, {5}, {6} ]
    overlapSets <- list()
    if (nrow(forIGV[[n]][[mod]])>0) {
      
      for (r in 1:nrow(forIGV[[n]][[mod]])) {
        overlapSets <- append(overlapSets, list(sets::set(as.numeric(r))))
      }
      #For each gene co-ordinate comparison [k, l]
      for (k in 1:nrow(forIGV[[n]][[mod]])) {
        for (l in 1:k) {
          
          # If the co-ordinate ranges overlap
          if (overlapsFunction(forIGV[[n]][[mod]][k, "start"], forIGV[[n]][[mod]][k, "end"], 
                               forIGV[[n]][[mod]][l, "start"], forIGV[[n]][[mod]][l, "end"])==TRUE) {
            
            # Find the indexes of the sets containing each range
            kIndex <- findItem(k, overlapSets)
            lIndex <- findItem(l, overlapSets)
            
            # No need to merge if the co-ordinate ranges are already in the same sets
            if (kIndex!=lIndex) {
              
              # If they are in different sets, merge the two sets, replacing the old ones
              newSet <- set_union(overlapSets[[kIndex]], overlapSets[[lIndex]])
              overlapSets <- overlapSets[-c(kIndex, lIndex)]
              overlapSets <- append(overlapSets, list(newSet))
            }
          }
        }
      }
    } 
    else next
    modOverlaps[[mod]] <- overlapSets
  }
  allOverlaps[[n]] <- modOverlaps
}

rm(modOverlaps, overlapSets, kIndex, lIndex, newSet)


# Find the maximum range for the overlapping epigenetic modifications.
for (n in names(allOverlaps)) {
  for (mod in epiMods) {
    if (length(allOverlaps[[n]][[mod]])>0) {
      
      for (l in 1:length(allOverlaps[[n]][[mod]])) {
        modStart <- c()
        modEnd <- c() 
        
        for (o in allOverlaps[[n]][[mod]][l]) {
          modStart <- append(modStart, forIGV[[n]][[mod]][as.numeric(o), "start"])
          modEnd <- append(modEnd, forIGV[[n]][[mod]][as.numeric(o), "end"])
          
          allOverlaps[[n]][[mod]][l] <- paste(min(modStart), max(modEnd), sep = "-")
        }
      }
    }
  }
}

rm(modStart, modEnd, o, l)


# Create dataframes with the information needed in the bed file.
for (n in names(forIGV)) {
  for (mod in epiMods) {
    df <- data.frame(seqname = numeric(),
                     ranges = character(),
                     strand = factor(),
                     epiMod = character(),
                     colour = character())
    
    if (length(allOverlaps[[n]][[mod]])>0) {
      for (l in 1:length(allOverlaps[[n]][[mod]])) {
        df <- rbind(df, data.frame(seqname = forIGV[[n]][[mod]][1,"seqnames"],
                                   ranges = allOverlaps[[n]][[mod]][[l]],
                                   strand = forIGV[[n]][[mod]][1,"strand"],
                                   epiMod = mod,
                                   colour = forIGV[[n]][[mod]][1,"itemRgb"]))
      }
    }
    allOverlaps[[n]][[mod]] <- df
  }
}

rm(df, l)


# Combine the dataframes for each epigenetic modification into one dataframe.
modBed <- data.frame(seqname = numeric(),
                     ranges = character(),
                     strand = factor(),
                     epiMod = character(),
                     colour = character())

for (n in names(allOverlaps)) {
  for (mod in epiMods) {
    modBed <- rbind(modBed, allOverlaps[[n]][[mod]])
    
  }
}

# Create bed file.
modBed <- GRanges(
  seqnames=Rle(modBed$seqname),
  ranges=IRanges(modBed$ranges),
  name=modBed$epiMod,
  itemRgb=modBed$colour)

# Export bed file.
rtracklayer::export.bed(modBed, "~/allNLRs.bed")
