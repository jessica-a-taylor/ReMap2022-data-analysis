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

# Import all Arabidopsis genes.
Atgenes <- as.data.frame(transcriptsBy(TxDb.Athaliana.BioMart.plantsmart28, by="gene"))

# Import ATAC-seq data (no treatment files). 
# This will be used to select eurchromatic genes and determine the average size of promoter regions.
sheets <- c(3,15,18)

openChromatin <- hash()
for (s in sheets) {
  openChromatin[[paste("ACR", s, sep = "")]] <-  as.data.frame(read_xlsx("C:\\Users\\jexy2\\OneDrive\\Documents\\PhD\\PhD reading\\Data\\ACRs paper.xlsx", sheet = s))
}

rm(sheets, s)

# From the openChromatin dataset, extract the rows corresponding to ATgene promotors.
openGenes <- openChromatin
for (s in names(openGenes)) {
  openGenes[[s]] <- openGenes[[s]][which(openGenes[[s]]$geneId %in% Atgenes$group_name &
                                                   grepl("Promoter", openGenes[[s]]$annotation)==TRUE),]
}

# Merge the data in to one big dataframe.
openGenesBed <- data.frame()

for (s in names(openGenes)) {
  openGenesBed <- rbind(openGenesBed, openGenes[[s]])
}

# Remove TEs from the openGenesBed dataframe.
transposableElements <- as.data.frame(read_xlsx("C:\\Users\\jexy2\\OneDrive\\Documents\\PhD\\Arabidopsis TE genes.xlsx"))

openGenesBed <- openGenesBed[-c(which(openGenesBed$geneId %in% transposableElements$Locus)),]

# Select a random sample of Atgenes from the openGenesBed dataframe.
controlGenes <- openGenesBed[c(sample(nrow(openGenesBed), 200)),]

# Get the coordinates for the gene bodies of each NLR.
NLRgenes <- as.data.frame(read_xlsx("C:\\Users\\jexy2\\OneDrive\\Documents\\PhD\\Arabidopsis NLRs.xlsx", sheet = 1))

# Remove duplicate genes (different versions).
Atgenes <- Atgenes[-c(which(Atgenes$tx_name == str_match(Atgenes$tx_name, "^([0-9a-zA-Z]+)([.])([2-9])$")[,1])),]


NLRgenebody <- data.frame(seqnames = numeric(),
                          start = numeric(),
                          end = numeric(),
                          width = numeric(),
                          strand = factor(),
                          tx_id = numeric(),
                          tx_name = character())

for (gene in NLRgenes$Gene) {
  NLRgenebody <- rbind(NLRgenebody, as.data.frame(Atgenes[grepl(gene, Atgenes$tx_name),]))
}

# Remove duplicate genes.
#NLRgenebody <- NLRgenebody[-c(which(NLRgenebody$tx_name == str_match(NLRgenebody$tx_name, "^([0-9a-zA-Z]+)([.])([2-9])$")[,1])),]

# Create a ranges column by merging the start and end columns.
NLRgenebody$ranges <- paste(NLRgenebody$start,"-",NLRgenebody$end, sep = "")

genebodyBed <- GRanges(
  seqnames=Rle(NLRgenebody$seqnames),
  ranges=IRanges(NLRgenebody$ranges),
  name=NLRgenebody$tx_name)

rtracklayer::export.bed(genebodyBed, "NLRgenebody.bed")


# Create new dataframes for chunks of the gene body (20% intervals of the gene length).
geneChunks <- hash(width20 = NLRgenebody[-c(4:6,10)], 
                  width40 = NLRgenebody[-c(4:6,10)], 
                  width60 = NLRgenebody[-c(4:6,10)], 
                  width80 = NLRgenebody[-c(4:6,10)], 
                  width100 = NLRgenebody[-c(4:6,10)])

geneWidth <- hash(width20 = c(), 
                   width40 = c(), 
                   width60 = c(), 
                   width80 = c(), 
                   width100 = c())

for (row in 1:nrow(NLRgenebody)) {
  geneWidth[["width20"]] <- append(geneWidth[["width20"]], paste(NLRgenebody[row,"start"],"-", NLRgenebody[row,"start"] + NLRgenebody[row,"width"]*0.2, sep = ""))
  geneWidth[["width40"]] <- append(geneWidth[["width40"]], paste(NLRgenebody[row,"start"] + NLRgenebody[row,"width"]*0.2+1,"-", NLRgenebody[row,"start"] + NLRgenebody[row,"width"]*0.4, sep = ""))
  geneWidth[["width60"]] <- append(geneWidth[["width60"]], paste(NLRgenebody[row,"start"] + NLRgenebody[row,"width"]*0.4+1,"-", NLRgenebody[row,"start"] + NLRgenebody[row,"width"]*0.6, sep = ""))
  geneWidth[["width80"]] <- append(geneWidth[["width80"]], paste(NLRgenebody[row,"start"] + NLRgenebody[row,"width"]*0.6+1,"-", NLRgenebody[row,"start"] + NLRgenebody[row,"width"]*0.8, sep = ""))
  geneWidth[["width100"]] <- append(geneWidth[["width100"]], paste(NLRgenebody[row,"start"] + NLRgenebody[row,"width"]*0.8+1,"-",NLRgenebody[row,"end"], sep = ""))
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
    name=geneChunks[[n]]$tx_name)
  
  rtracklayer::export.bed(geneBed, paste("NLR", n, ".bed", sep = ""))
}

# Create new dataframe for the coordinates of the regions 200bp downstream of the TTS.

downstreamRegion <- c()
for (row in 1:nrow(NLRgenebody)) {
  downstreamRegion <- append(downstreamRegion, paste(NLRgenebody[row,"end"],"-",NLRgenebody[row,"end"]+200, sep = ""))
}

NLRdownstream <- NLRgenebody[-c(4:6,10)]

NLRdownstream$ranges <- downstreamRegion
rm(downstreamRegion)

NLRdownstream$start <- str_match(NLRdownstream$ranges, "^([0-9]+)(-)([0-9]+)$")[,2]
NLRdownstream$end <- str_match(NLRdownstream$ranges, "^([0-9]+)(-)([0-9]+)$")[,4]

downstreamBed <- GRanges(
  seqnames=Rle(NLRdownstream$seqnames),
  ranges=IRanges(NLRdownstream$ranges),
  name=NLRdownstream$tx_name)

rtracklayer::export.bed(downstreamBed, "NLRdownstream.bed")


# From the openChromatin dataset, extract the rows corresponding to NLR promotors.
openPromotors <- openChromatin
for (s in names(openChromatin)) {
  openPromotors[[s]] <- openPromotors[[s]][which(openPromotors[[s]]$geneId %in% NLRgenes$Gene &
                                                   grepl("Promoter", openPromotors[[s]]$annotation)==TRUE),]
}

# Merge the data in to one big dataframe.
openPromotorsBed <- data.frame()

for (s in names(openPromotors)) {
  openPromotorsBed <- rbind(openPromotorsBed, openPromotors[[s]])
}

rm(openPromotors, openChromatin)

# Get the min and max coordinates for open chromatin at the promoters of each NLR.
openPromoters <- data.frame(seqnames = character(),
                            name = character(),
                            start = numeric(),
                            end = numeric(),
                            geneLength = numeric())

for (s in unique(openPromotorsBed$geneId)) {
  df <- openPromotorsBed[openPromotorsBed$geneId==s,]
  
  openPromoters <- rbind(openPromoters, data.frame(seqnames = df[1,"seqnames"],
                                                   name = s,
                                                   start = min(df$start),
                                                   end = max(df$end),
                                                   geneLength = df[1,"geneLength"]))
}


# Merge start and end coordinates columns to create a ranges column.
openPromoters$ranges <- paste(openPromoters$start,"-",openPromoters$end, sep = "")

promoterSize <- c()
# Plot a histogram for the promoter sizes.
for (row in 1:nrow(openPromoters)) {
  promoterSize <- append(promoterSize, openPromoters[row, "end"]-openPromoters[row, "start"])
}

hist(promoterSize) # most promotors are less than 1kb

# Create bed file.
openPromotersBed <- GRanges(
  seqnames=Rle(openPromoters$seqnames),
  ranges=IRanges(openPromoters$ranges),
  name=openPromoters$name)

# Export bed file.
rtracklayer::export.bed(openPromotersBed, "~/openPromoters.bed")


# Get the coordinates for the promotors of each NLR.
ATpromotors500 <- promoters(TxDb.Athaliana.BioMart.plantsmart28, upstream=500, downstream=0, use.names = TRUE)
ATpromotors1000 <- promoters(TxDb.Athaliana.BioMart.plantsmart28, upstream=1000, downstream=0, use.names = TRUE)

NLRpromotor500 <- data.frame(seqnames = numeric(),
                          start = numeric(),
                          end = numeric(),
                          width = numeric(),
                          strand = factor(),
                          tx_id = numeric(),
                          tx_name = character())

NLRpromotor1000 <- NLRpromotor500

for (gene in NLRgenes$Gene) {
  NLRpromotor500 <- rbind(NLRpromotor500, as.data.frame(ATpromotors500[grepl(gene,ATpromotors500$tx_name),]))
  NLRpromotor1000 <- rbind(NLRpromotor1000, as.data.frame(ATpromotors1000[grepl(gene,ATpromotors1000$tx_name),]))
}

# Create a ranges column by merging the start and end columns.
NLRpromotor500$ranges <- paste(NLRpromotor500$start,"-",NLRpromotor500$end, sep = "")
NLRpromotor1000$ranges <- paste(NLRpromotor1000$start,"-",NLRpromotor1000$end, sep = "")

# Remove duplicate genes.
NLRpromotor500 <- NLRpromotor500[-c(which(NLRpromotor500$tx_name == str_match(NLRpromotor500$tx_name, "^([0-9a-zA-Z]+)([.])([2-9])$")[,1])),]
NLRpromotor1000 <- NLRpromotor1000[-c(which(NLRpromotor1000$tx_name == str_match(NLRpromotor1000$tx_name, "^([0-9a-zA-Z]+)([.])([2-9])$")[,1])),]

# Add a new column for the gene name, removing ".1" from the end.
NLRpromotor500$group_name <- str_match(NLRpromotor500$tx_name, "^([0-9a-zA-Z]+)([.])([1])$")[,2]

promotor500Bed <- GRanges(
  seqnames=Rle(NLRpromotor500$seqnames),
  ranges=IRanges(NLRpromotor500$ranges),
  name=NLRpromotor500$tx_name)

rtracklayer::export.bed(promotor500Bed, "NLRpromotor500.bed")

NLRpromotor1000$group_name <- str_match(NLRpromotor1000$tx_name, "^([0-9a-zA-Z]+)([.])([1])$")[,2]

promotor1000Bed <- GRanges(
  seqnames=Rle(NLRpromotor1000$seqnames),
  ranges=IRanges(NLRpromotor1000$ranges),
  name=NLRpromotor1000$tx_name)

rtracklayer::export.bed(promotor1000Bed, "NLRpromotor1000.bed")

rm(ATpromotors500, ATpromotors1000)

# Get the coordinates for the upstream intergenic regions of each NLR.
usCoordinates <- c()

for (gene in NLRgenes$Gene) {
  currentGene <- which(Atgenes$group_name==gene)

  if (Atgenes[currentGene, "strand"]=="+") {
    previousGene <- currentGene - 1
    
    if (Atgenes[previousGene, "strand"]=="+") {
      distance <- (Atgenes[currentGene, "start"] - 1001) - (Atgenes[previousGene, "end"] + 201)
      
      if (distance > 0) {
        usCoordinates <- append(usCoordinates, paste(Atgenes[previousGene, "end"] + 201, "-", Atgenes[previousGene, "end"] + 201 + distance, sep = "")) 
      } 
      else usCoordinates <- append(usCoordinates, NA)
    }
    
    else if (Atgenes[previousGene, "strand"]=="-") {
      distance <- (Atgenes[currentGene, "start"] - 1001) - (Atgenes[previousGene, "end"] + 1001)
      
      if (distance > 0) {
        usCoordinates <- append(usCoordinates, paste(Atgenes[previousGene, "end"] + 1001, "-", Atgenes[previousGene, "end"] + 1001 + distance, sep = "")) 
      } 
      else usCoordinates <- append(usCoordinates, NA)
      
    }
  }
  
  else if (Atgenes[currentGene, "strand"]=="-") {
    previousGene <- currentGene + 1
    
    if (Atgenes[previousGene, "strand"]=="+") {
      distance <- (Atgenes[previousGene, "start"] - 1001) - (Atgenes[currentGene, "end"] + 1001)
      
      if (distance > 0) {
        usCoordinates <- append(usCoordinates, paste(Atgenes[previousGene, "start"] - 1001 - distance, "-", Atgenes[previousGene, "start"] - 1001, sep = "")) 
      } 
      else usCoordinates <- append(usCoordinates, NA)
      
    }
    
    else if (Atgenes[previousGene, "strand"]=="-") {
      distance <- (Atgenes[previousGene, "start"] - 201) - (Atgenes[currentGene, "end"] + 1001)
      if (distance > 0) {
        usCoordinates <- append(usCoordinates, paste(Atgenes[previousGene, "start"] - 201 - distance, "-", Atgenes[previousGene, "start"] - 201, sep = "")) 
      } 
      else usCoordinates <- append(usCoordinates, NA)
      
    }
  } 
}

upstreamIntergenic <- Atgenes[which(Atgenes$group_name %in% NLRgenes$Gene),]
upstreamIntergenic$ranges <- usCoordinates 

upstreamIntergenicBed <- upstreamIntergenic[-c(which(is.na(upstreamIntergenic$ranges))),]

upstreamIntergenicBed <- GRanges(
  seqnames=Rle(as.numeric(upstreamIntergenic$seqnames)),
  ranges=IRanges(upstreamIntergenic$ranges),
  name=upstreamIntergenic$group_name)

rtracklayer::export.bed(upstreamIntergenicBed, "upstreamIntergenic.bed")


# Get the coordinates for the downstream intergenic regions of each NLR.
dsCoordinates <- c()

for (gene in NLRgenes$Gene) {
  currentGene <- which(Atgenes$group_name==gene)
  
  if (Atgenes[currentGene, "strand"]=="-") {
    nextGene <- currentGene - 1
    
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
  
  else if (Atgenes[currentGene, "strand"]=="+") {
    nextGene <- currentGene + 1
    
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
}

downstreamIntergenic <- Atgenes[which(Atgenes$group_name %in% NLRgenes$Gene),]
downstreamIntergenic$ranges <- dsCoordinates 

downstreamIntergenicBed <- downstreamIntergenic[-c(which(is.na(downstreamIntergenic$ranges))),]

downstreamIntergenicBed <- GRanges(
  seqnames=Rle(as.numeric(downstreamIntergenic$seqnames)),
  ranges=IRanges(downstreamIntergenic$ranges),
  name=downstreamIntergenic$group_name)

rtracklayer::export.bed(downstreamIntergenicBed, "downstreamIntergenic.bed")


# Import ReMap2022 data.
ReMap <- rtracklayer::import.bed("C:\\Users\\jexy2\\OneDrive\\Documents\\PhD\\remap2022_histone_all_macs2_TAIR10_v1_0.bed.gz")

# Convert to a dataframe and define column names.
ReMap <- as.data.frame(ReMap, colnames = c("seqnames", "start", "end", "width",
                                           "strand", "name", "score", "itemRgb",
                                           "thick.start", "thick.end", "thick.width"))

# Remove unwanted columns.
ReMap <- ReMap[,-c(9:11)]

# Create a hash to store ReMap data for each NLR/cluster.
colnames(NLRgenebody)[2] <- "Gene"
NLR_hash <- hash()

for (row in 1:nrow(NLRgenebody)) {
  # Select rows that are within the range of each NLR/cluster and on the same chromosome.
  ReMapRows <- c(which(ReMap[,"start"] > NLRgenebody[row, "start"]-5000 & ReMap[,"end"] < NLRgenebody[row, "end"]+5000 & ReMap[,"seqnames"] == as.numeric(NLRgenebody[row, "seqnames"])))
  NLR_hash[[NLRgenebody[row,"Gene"]]] <- ReMap[ReMapRows,]
}

rm(ReMap, ReMapRows)

for(n in names(NLR_hash)) {
  if (nrow(NLR_hash[[n]])>=1) {
    # Run regex on name column, extracting each section
    # (experiment, epigenetic modification, ecotype, other info)
    NLR_hash[[n]][c("exp.", "epiMod", "ecotype", "info")] <- str_match(NLR_hash[[n]][,"name"], "^([0-9a-zA-Z]+)\\.([0-9a-zA-Z-]+)\\.([0-9a-zA-Z-]+)[_\\.](.*)$")[,-1]
    
    # Filter epiMod column, excluding unwanted modifications
    NLR_hash[[n]] <- NLR_hash[[n]][!NLR_hash[[n]]$epiMod %in% c("H3", "HTR12", "H2A", "H2B", "H3T3ph", "H1", "H4K16ac", "H2A-X",
                                                                "H2AV", "HTA6", "H3-1", "H4K12ac", "H4K8ac", "H3K5ac", "H4K5ac") & 
                                     !NLR_hash[[n]]$ecotype %in% c("C24", "undef", "Col-x-Ler", "Ler-x-Col", "Col-x-C24"),]
    
    NLR_hash[[n]] <- NLR_hash[[n]][,-6]
    
    # Filter info column, excluding unwanted conditions (too old, too young, wrong part of plant, etc)
    NLR_hash[[n]] <- NLR_hash[[n]][!grepl("mutant",NLR_hash[[n]]$info) & !grepl("mature",NLR_hash[[n]]$info) & !grepl("senescent",NLR_hash[[n]]$info) & 
                                     !grepl("inflorescence",NLR_hash[[n]]$info) & !grepl("drought",NLR_hash[[n]]$info) & !grepl("old",NLR_hash[[n]]$info) & 
                                     !grepl("min",NLR_hash[[n]]$info) & !grepl("endosperm",NLR_hash[[n]]$info) & !grepl("-se-",NLR_hash[[n]]$info) &
                                     !grepl("-TSA-",NLR_hash[[n]]$info) & !grepl("-GSNO-",NLR_hash[[n]]$info) & !grepl("flg22",NLR_hash[[n]]$info) &
                                     !grepl("transgenic",NLR_hash[[n]]$info) & !grepl("GSH",NLR_hash[[n]]$info) &
                                     !grepl("-acc1",NLR_hash[[n]]$info) & !grepl("-ethylene",NLR_hash[[n]]$info) & !grepl("-C2H4",NLR_hash[[n]]$info) &
                                     !grepl("leaves_3w-K36M-homoz",NLR_hash[[n]]$info) & !grepl("undef_seedling_10d-h3-1kd-1",NLR_hash[[n]]$info) &
                                     !grepl("-air",NLR_hash[[n]]$info) & !grepl("-ehylene",NLR_hash[[n]]$info) & !grepl("-swap",NLR_hash[[n]]$info) &
                                     !grepl("-K36M",NLR_hash[[n]]$info) & !grepl("-H3-KD",NLR_hash[[n]]$info) & !grepl("-water",NLR_hash[[n]]$info) &
                                     !grepl("undef_seedling_10d-h3-1kd-2",NLR_hash[[n]]$info) & !grepl("seedling_3d-wt-ehylene",NLR_hash[[n]]$info) &
                                     !grepl("GSE67322",NLR_hash[[n]]$exp.) & !grepl("GSE42695",NLR_hash[[n]]$exp.) & !grepl("GSE75071",NLR_hash[[n]]$exp.) &
                                     !grepl("GSE62615",NLR_hash[[n]]$exp.) & !grepl("GSE103361",NLR_hash[[n]]$exp.) & !grepl("GSE50636",NLR_hash[[n]]$exp.) &
                                     !grepl("GSE93223",NLR_hash[[n]]$exp.) & !grepl("GSE37644",NLR_hash[[n]]$exp.) & !grepl("GSE108414",NLR_hash[[n]]$exp.) &
                                     !grepl("GSE22276",NLR_hash[[n]]$exp.) & !grepl("GSE89768",NLR_hash[[n]]$exp.) & !grepl("GSE117391",NLR_hash[[n]]$exp.),] 
    
    # Filter info column, checking written plant ages and removing BAD AGES
    
    # For each row
    for (row in nrow(NLR_hash[[n]]):1) {
      # Find the age if it exists (in the format 1w / 8d / 10h)
      matches <- str_match(NLR_hash[[n]][row, "info"], "_([0-9]+)([dwh])")
      
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
          NLR_hash[[n]] <- NLR_hash[[n]][-row,]
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


for(n in names(NLR_hash)) {
  # Extract root data from NLR_hash and store in RootNLRs.
  RootNLRs[[n]] <- NLR_hash[[n]][grepl("roots",NLR_hash[[n]]$info),]
  # Remove root data from NLR_hash.
  LeafNLRs[[n]] <- NLR_hash[[n]][!grepl("roots",NLR_hash[[n]]$info),]
}

# Delete NLR_hash.
rm(NLR_hash)

# Create hashes for WT and mutant root data.
mutantRootNLRs <- hash()
WTRootNLRs <- hash()

for (n in names(RootNLRs)) {
  WTRootNLRs[[n]] <- RootNLRs[[n]][!grepl("-hag1",RootNLRs[[n]]$info) & !grepl("-OTU5",RootNLRs[[n]]$info),]
  mutantRootNLRs[[n]] <- RootNLRs[[n]][grepl("-hag1",RootNLRs[[n]]$info) & grepl("-OTU5",RootNLRs[[n]]$info),]
}

rm(RootNLRs)

# Create hashes for WT and mutant leaf data.
mutantLeafNLRs <- hash()
WTLeafNLRs <- hash()

allConditions <- c()
for (n in names(LeafNLRs)) {
  allConditions <- append(allConditions, unique(LeafNLRs[[n]]$info))
}

allConditions <- unique(allConditions)
WTConditions <- allConditions[c(1,2,3,5,6,7,9,12,13,14,15,17,19,20)]
mutantConditions <- allConditions[-c(1,2,3,5,6,7,9,12,13,14,15,17,19,20)]
rm(allConditions)

for (n in names(LeafNLRs)) {
  WTLeafNLRs[[n]] <- LeafNLRs[[n]][LeafNLRs[[n]]$info %in% c(WTConditions),]
  mutantLeafNLRs[[n]] <- LeafNLRs[[n]][LeafNLRs[[n]]$info %in% c(mutantConditions),]
}

rm(LeafNLRs, WTConditions, mutantConditions)

# Create hashes for leaf data in Col and Ler ecotypes.
ColWTLeafNLRs <- hash()
LerWTLeafNLRs <- hash()

for (n in names(WTLeafNLRs)) {
  ColWTLeafNLRs[[n]] <- WTLeafNLRs[[n]][grepl("Col-0",WTLeafNLRs[[n]]$ecotype),]
  LerWTLeafNLRs[[n]] <- WTLeafNLRs[[n]][grepl("Ler",WTLeafNLRs[[n]]$ecotype),]
}

rm(WTLeafNLRs)


# Select tissue type for analysis.
dataToUse <- WTRootNLRs

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


allOverlaps <- hash()

# For each epigenetic modification name
for (n in names(ColWTdata)) {
  modOverlaps <- hash()
  
  for (mod in epiMods) {
    
    # Generate overlapSets as a list of single-item sets
    # eg, [ {1}, {2}, {3}, {4}, {5}, {6} ]
    overlapSets <- list()
    if (nrow(ColWTdata[[n]][[mod]])>0) {
      
      for (r in 1:nrow(ColWTdata[[n]][[mod]])) {
        overlapSets <- append(overlapSets, list(sets::set(as.numeric(r))))
      }
      #For each gene co-ordinate comparison [k, l]
      for (k in 1:nrow(ColWTdata[[n]][[mod]])) {
        for (l in 1:k) {
          
          # If the co-ordinate ranges overlap
          if (overlapsFunction(ColWTdata[[n]][[mod]][k, "start"], ColWTdata[[n]][[mod]][k, "end"], 
                               ColWTdata[[n]][[mod]][l, "start"], ColWTdata[[n]][[mod]][l, "end"])==TRUE) {
            
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
          modStart <- append(modStart, ColWTdata[[n]][[mod]][as.numeric(o), "start"])
          modEnd <- append(modEnd, ColWTdata[[n]][[mod]][as.numeric(o), "end"])
          
          allOverlaps[[n]][[mod]][l] <- paste(min(modStart), max(modEnd), sep = "-")
        }
      }
    }
  }
}

rm(modStart, modEnd, o, l)


# Create dataframes with the information needed in the bed file.
for (n in names(ColWTdata)) {
  for (mod in epiMods) {
    df <- data.frame(seqname = numeric(),
                     ranges = character(),
                     strand = factor(),
                     epiMod = character(),
                     colour = character())
    
    if (length(allOverlaps[[n]][[mod]])>0) {
      for (l in 1:length(allOverlaps[[n]][[mod]])) {
        df <- rbind(df, data.frame(seqname = ColWTdata[[n]][[mod]][1,"seqnames"],
                                   ranges = allOverlaps[[n]][[mod]][[l]],
                                   strand = ColWTdata[[n]][[mod]][1,"strand"],
                                   epiMod = mod,
                                   colour = ColWTdata[[n]][[mod]][1,"itemRgb"]))
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



# Create a hash with the data on the coordinates of each gene region.
# First rename the columns in each dataset so they can be indexed the same way.
colnames(NLRpromotor500)[9] ="Gene"
colnames(NLRpromotor1000)[9] ="Gene"
colnames(NLRdownstream)[2] ="Gene"
colnames(upstreamIntergenic)[2] ="Gene"
colnames(downstreamIntergenic)[2] ="Gene"

for (n in names(geneChunks)) {
  colnames(geneChunks[[n]])[2]="Gene"
}


regions <- hash(UpstreamIntergenic = upstreamIntergenic, Promotor1000 = NLRpromotor1000, Promotor500 = NLRpromotor500,
                Gene20 = geneChunks[["width20"]],  Gene40 = geneChunks[["width40"]], Gene60 = geneChunks[["width60"]], 
                Gene80 = geneChunks[["width80"]], Gene100 = geneChunks[["width100"]],Downstream = NLRdownstream,
                DownstreamIntergenic = downstreamIntergenic)

# Create a dictionary containing the frequency of each chromatin modification occurring in each region of each NLR.
modsPerRegion <- hash()

for (r in names(regions)) {
  # Create a hash containing a list of chromatin modifications overlapping with the NLR region.
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

# Plot the percentage of NLRs with each chromatin modification within the gene body.
modFrequenciesBarPlot <- ggplot(modFrequenciesDF, aes(x=Region, y=Frequency)) + facet_wrap(~Modification) +
  geom_bar(stat = "identity", position = "dodge") + scale_fill_brewer(palette = "RdYlBu") +
  theme_minimal() + labs(x = "Gene Region", y = "Frequency of occurrence (%)")

modFrequenciesLinePlot <- ggplot(modFrequenciesDF, aes(x=factor(Region, level = level), y=Frequency)) + 
  geom_line(aes(group = Modification, color = Modification)) + theme_classic() + 
  labs(x = "Gene region", y = "Frequency of occurrence (%)")

axisText <- c("Intergenic", "Promotor \n(1kb)", "Promotor \n(500bp)", "TSS",
              "20%", "40%", "60%", "80%", "100%", 
              "Downstream \n(200bp)", "Intergenic")

modFrequenciesPlot <- ggplot(modFrequenciesDF, aes(x = axisGroup, y = Frequency)) + 
  scale_x_continuous(limits = c(-60, 140), breaks = seq(-60, 140, 20), labels = axisText) +
  geom_line(aes(group = Modification, color = Modification),size = 1.3) +
  geom_point(aes(group = Modification, color = Modification), size = 2) + theme_minimal() + 
  labs(x = "", y = "Frequency of occurrence (%)") +
  geom_vline(xintercept=0, color="grey", size=1) +
  coord_cartesian(ylim= c(0,100), clip = "off") + theme(plot.margin = unit(c(1,1,2,1), "lines")) +
  annotation_custom(textGrob("% of gene length from TSS", gp=gpar(fontsize=12, col = "grey33")),xmin=0,xmax=100,ymin=-14,ymax=-14) + 
  annotation_custom(textGrob("Gene region", gp=gpar(fontsize=14)),xmin=0,xmax=100,ymin=-20,ymax=-20) +
  theme(axis.text.x = element_text(size = 11, colour = "black"), axis.text.y = element_text(size = 12,colour = "black"), 
        axis.title.y = element_text(size = 14, vjust = 2))

modFrequenciesPlot