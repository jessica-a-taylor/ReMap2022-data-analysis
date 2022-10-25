library(readxl)
library(karyoploteR)
library(rtracklayer)
library(dplyr)
library(stringr)
library(hash)
library(sets)
library(TxDb.Athaliana.BioMart.plantsmart28)
library(ggplot2)

# Get the coordinates for the gene bodies of each NLR.
NLRgenes <- as.data.frame(read_xlsx("C:\\Users\\jexy2\\OneDrive\\Documents\\PhD\\Arabidopsis NLRs.xlsx", sheet = 1))

Atgenes <- as.data.frame(transcriptsBy(TxDb.Athaliana.BioMart.plantsmart28, by="gene"))

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

# Create a ranges column by merging the start and end columns.
NLRgenebody$ranges <- paste(NLRgenebody$start,"-",NLRgenebody$end, sep = "")

# Remove duplicate genes.
NLRgenebody <- NLRgenebody[-c(which(NLRgenebody$tx_name == str_match(NLRgenebody$tx_name, "^([0-9a-zA-Z]+)([.])([2-9])$")[,1])),]

genebodyBed <- GRanges(
  seqnames=Rle(NLRgenebody$seqnames),
  ranges=IRanges(NLRgenebody$ranges),
  name=NLRgenebody$tx_name)

rtracklayer::export.bed(genebodyBed, "NLRgenebody.bed")

rm(Atgenes)

# Get the coordinates for the promotors of each NLR.
ATpromotors <- promoters(TxDb.Athaliana.BioMart.plantsmart28, upstream=500, downstream=200, use.names = TRUE)

NLRpromotor <- data.frame(seqnames = numeric(),
                          start = numeric(),
                          end = numeric(),
                          width = numeric(),
                          strand = factor(),
                          tx_id = numeric(),
                          tx_name = character())

for (gene in NLRgenes$Gene) {
  NLRpromotor <- rbind(NLRpromotor, as.data.frame(ATpromotors[grepl(gene,ATpromotors$tx_name),]))
}

# Create a ranges column by merging the start and end columns.
NLRpromotor$ranges <- paste(NLRpromotor$start,"-",NLRpromotor$end, sep = "")

# Remove duplicate genes.
NLRpromotor <- NLRpromotor[-c(which(NLRpromotor$tx_name == str_match(NLRpromotor$tx_name, "^([0-9a-zA-Z]+)([.])([2-9])$")[,1])),]

# Add a new column for the gene name, removing ".1" from the end.
NLRpromotor$group_name <- str_match(NLRpromotor$tx_name, "^([0-9a-zA-Z]+)([.])([1])$")[,2]

promotorsBed <- GRanges(
  seqnames=Rle(NLRpromotor$seqnames),
  ranges=IRanges(NLRpromotor$ranges),
  name=NLRpromotor$tx_name)

rtracklayer::export.bed(promotorsBed, "NLRpromotor.bed")

rm(ATpromotors)

# Get the coordinates for the upstream and downstream intergenic regions of each NLR.
upstreamIntergenic <- data.frame(Chrom = NLRgenes$Chromosome,
                                 start = NLRgenes$`Upstream Intergenic Start`,
                                 end = NLRgenes$`Upstream Intergenic End`,
                                 Gene = NLRgenes$Gene)

upstreamIntergenicBed <- upstreamIntergenic[-c(which(is.na(upstreamIntergenic$start))),]

upstreamIntergenicBed$upstreamRanges <- paste(upstreamIntergenicBed$start,"-",upstreamIntergenicBed$end, sep = "")

upstreamIntergenicBed <- GRanges(
  seqnames=Rle(upstreamIntergenicBed$Chrom),
  ranges=IRanges(upstreamIntergenicBed$upstreamRanges),
  name=upstreamIntergenicBed$Gene)

rtracklayer::export.bed(upstreamIntergenicBed, "upstreamIntergenic.bed")

downstreamIntergenic <- data.frame(Chrom = NLRgenes$Chromosome,
                                   start = NLRgenes$`Downstream Intergenic Start`,
                                   end = NLRgenes$`Downstream Intergenic End`,
                                   Gene = NLRgenes$Gene)

downstreamIntergenicBed <- downstreamIntergenic[-c(which(is.na(downstreamIntergenic$start))),]

downstreamIntergenicBed$downstreamRanges <- paste(downstreamIntergenicBed$start,"-",downstreamIntergenicBed$end, sep = "")

downstreamIntergenicBed <- GRanges(
  seqnames=Rle(downstreamIntergenicBed$Chrom),
  ranges=IRanges(downstreamIntergenicBed$downstreamRanges),
  name=downstreamIntergenicBed$Gene)

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

# Create list of chromatin modifications.
epiMods <- c()
for (n in names(ColWTLeafNLRs)) {
  epiMods <- append(epiMods, unique(ColWTLeafNLRs[[n]]$epiMod))
}

epiMods <- unique(epiMods)


# Create a dictionary (hash containing hashes) with dataframes for each epiMod for each NLR.
ColWTLeafData <- hash()

for (n in names(ColWTLeafNLRs)) {
  modHash <- hash()
  
  for (mod in epiMods) {
    modHash[[mod]] <- ColWTLeafNLRs[[n]][ColWTLeafNLRs[[n]]$epiMod==mod,]
  }
  ColWTLeafData[[n]] <- modHash
}

rm(ColWTLeafNLRs, modHash)

# Merge start and end coordinates columns to create a ranges column.
for (n in names(ColWTLeafData)) {
  for (mod in epiMods) {
    if (nrow(ColWTLeafData[[n]][[mod]]) >= 1) {
      ColWTLeafData[[n]][[mod]]$ranges <- paste(ColWTLeafData[[n]][[mod]]$start,"-",ColWTLeafData[[n]][[mod]]$end, sep = "")
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
for (n in names(ColWTLeafData)) {
  modOverlaps <- hash()
  
  for (mod in epiMods) {
    
    # Generate overlapSets as a list of single-item sets
    # eg, [ {1}, {2}, {3}, {4}, {5}, {6} ]
    overlapSets <- list()
    if (nrow(ColWTLeafData[[n]][[mod]])>0) {
      
      for (r in 1:nrow(ColWTLeafData[[n]][[mod]])) {
        overlapSets <- append(overlapSets, list(set(as.numeric(r))))
      }
      #For each gene co-ordinate comparison [k, l]
      for (k in 1:nrow(ColWTLeafData[[n]][[mod]])) {
        for (l in 1:k) {
          
          # If the co-ordinate ranges overlap
          if (overlapsFunction(ColWTLeafData[[n]][[mod]][k, "start"], ColWTLeafData[[n]][[mod]][k, "end"], 
                               ColWTLeafData[[n]][[mod]][l, "start"], ColWTLeafData[[n]][[mod]][l, "end"])==TRUE) {
            
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
          modStart <- append(modStart, ColWTLeafData[[n]][[mod]][as.numeric(o), "start"])
          modEnd <- append(modEnd, ColWTLeafData[[n]][[mod]][as.numeric(o), "end"])
          
          allOverlaps[[n]][[mod]][l] <- paste(min(modStart), max(modEnd), sep = "-")
        }
      }
    }
  }
}

rm(modStart, modEnd, o, l)


# Create dataframes with the information needed in the bed file.
for (n in names(ColWTLeafData)) {
  for (mod in epiMods) {
    df <- data.frame(seqname = numeric(),
                     ranges = character(),
                     strand = factor(),
                     epiMod = character(),
                     colour = character())
    
    if (length(allOverlaps[[n]][[mod]])>0) {
      for (l in 1:length(allOverlaps[[n]][[mod]])) {
        df <- rbind(df, data.frame(seqname = ColWTLeafData[[n]][[mod]][1,"seqnames"],
                                   ranges = allOverlaps[[n]][[mod]][[l]],
                                   strand = ColWTLeafData[[n]][[mod]][1,"strand"],
                                   epiMod = mod,
                                   colour = ColWTLeafData[[n]][[mod]][1,"itemRgb"]))
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
colnames(NLRpromotor)[9] ="Gene"
colnames(NLRgenebody)[2] ="Gene"


regions <- hash(Promotor = NLRpromotor, GeneBody = NLRgenebody, UpstreamIntergenic = as.data.frame(upstreamIntergenic), 
                DownstreamIntergenic = as.data.frame(downstreamIntergenic))

# Create a dictionary containing the frequency of each chromatin modification occurring in each region of each NLR.
modsPerRegion <- hash()

for (r in names(regions)) {
  # Create a hash containing a list of chromatin modifications overlapping with the NLR region.
  GBmod <- hash()
  
  for (n in names(ColWTLeafData)) {
    modList <- c()
    
    for (mod in epiMods) {
      modPresent <- FALSE
      
      if (nrow(ColWTLeafData[[n]][[mod]]) >= 1 & !is.na(regions[[r]][regions[[r]]$Gene==n,]$start) & regions[[r]][regions[[r]]$Gene==n,]$end) {
        for (row in 1:nrow(ColWTLeafData[[n]][[mod]])) {
          if (overlapsFunction(ColWTLeafData[[n]][[mod]][row, "start"], ColWTLeafData[[n]][[mod]][row, "end"],
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
modFrequenciesDF <- data.frame(Region = character(),
                               Modification = character(),
                               Frequency = numeric())

for (r in names(modsPerRegion)) {
  modFrequenciesDF <- rbind(modFrequenciesDF, modsPerRegion[[r]])
}

# Plot the percentage of NLRs with each chromatin modification within the gene body.
modFrequenciesPlot <- ggplot(modFrequenciesDF, aes(x=Modification, y=Frequency, fill=Region)) + 
  geom_bar(stat = "identity", position = "dodge") + scale_fill_brewer(palette = "Spectral") +
  theme_minimal() + labs(x = "Chromatin Modification", y = "Frequency of occurance in NLR gene bodies (%)")




# Import chromatin states dataset.
chromStatesHash <- hash()
for(i in 1:9) {
  chromStatesHash[[paste("state", i, sep = "")]] <- read_xlsx("C:\\Users\\jexy2\\OneDrive\\Documents\\PhD\\Sequeira-Mendes ChromStates.xlsx", sheet = i)
}

chromStates <- data.frame(Chrom = numeric(),
                          Start = character(),
                          End = character(),
                          State = numeric())

for (n in names(chromStatesHash)) {
  df <- data.frame(Chrom = chromStatesHash[[n]]$Chrom,
                   Start = chromStatesHash[[n]]$start,
                   End = chromStatesHash[[n]]$end,
                   State = as.numeric(str_extract(n, "([0-9]+)")))
  
  chromStates <- rbind(chromStates, df)
}

# Merge start and end coordinates columns to create a ranges column.
chromStates$ranges <- paste(chromStates$Start,"-",chromStates$End, sep = "")

# Create a bed file for chromatin states dataset.
chromStatesBed <- GRanges(
  seqnames=Rle(chromStates$Chrom),
  ranges=IRanges(chromStates$ranges),
  name=chromStates$State)

# Export bed file.
rtracklayer::export.bed(chromStatesBed, "~/CS.bed")



# Import chromatin states dataset, sheet 2.
WTacr <- as.data.frame(read_xlsx("C:\\Users\\jexy2\\OneDrive\\Documents\\PhD\\RPP5 ChromStates.xlsx", sheet = 2))

# Merge start and end coordinates columns to create a ranges column.
WTacr$ranges = paste(WTacr$start,"-",WTacr$end, sep = "")

# Create a bed file for WT ACR dataset.
WTacrBed <- GRanges(
  seqnames=Rle("chr4",nrow(WTacr)),
  ranges=IRanges(WTacr$ranges))

# Export bed file.
rtracklayer::export.bed(WTacrBed, "~/WTacr.bed")

# Import chromatin states dataset, sheet 3.
Macr <- as.data.frame(read_xlsx("C:\\Users\\jexy2\\OneDrive\\Documents\\PhD\\RPP5 ChromStates.xlsx", sheet = 3))

# Merge start and end coordinates columns to create a ranges column.
Macr$ranges = paste(Macr$start,"-",Macr$end, sep = "")

# Create a bed file for WT ACR dataset.
MacrBed <- GRanges(
  seqnames=Rle("chr4",nrow(Macr)),
  ranges=IRanges(Macr$ranges))

# Export bed file.
rtracklayer::export.bed(MacrBed, "~/Mutant acr.bed")

# Import chromatin states dataset, sheet 4.
RootNEassociation <- as.data.frame(read_xlsx("C:\\Users\\jexy2\\OneDrive\\Documents\\PhD\\RPP5 ChromStates.xlsx", sheet = 4))

# Merge start and end coordinates columns to create a ranges column.
RootNEassociation$ranges = paste(RootNEassociation$start,"-",RootNEassociation$end, sep = "")

# Create a bed file for WT ACR dataset.
RootNEassociationBed <- GRanges(
  seqnames=Rle("chr4",nrow(RootNEassociation)),
  ranges=IRanges(RootNEassociation$ranges))

# Export bed file.
rtracklayer::export.bed(RootNEassociationBed, "~/RootNEassociation.bed")