install.packages(c("readxl", "karyoploteR", "rtracklayer"))
library(readxl)
library(karyoploteR)
library(rtracklayer)

# Import chromatin states dataset.
chromStates <- as.data.frame(read_xlsx("C:\\Users\\jexy2\\OneDrive\\Documents\\PhD\\RPP5 ChromStates.xlsx"))

# Merge start and end coordinates columns to create a ranges column.
chromStates$ranges <- paste(chromStates$start,"-",chromStates$end, sep = "")
chromStates

# Create a bed file for chromatin states dataset.
chromStatesBed <- GRanges(
  seqnames=Rle("chr4",nrow(chromStates)),
  ranges=IRanges(chromStates$ranges),
  name=chromStates$state)

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




# Import all ReMap2022 ChIP data (bed file).
ReMap <- rtracklayer::import.bed("C:\\Users\\jexy2\\OneDrive\\Documents\\PhD\\remap2022_histone_all_macs2_TAIR10_v1_0.bed.gz")

# Convert to a dataframe and define column names.
ReMap <- as.data.frame(ReMap, colnames = c("seqnames", "start", "end", "width",
                                           "strand", "name", "score", "itemRgb",
                                           "thick.start", "thick.end", "thick.width"))

# Remove unwanted columns.
ReMap <- ReMap[,-c(9:11)]

# Import coordainated of all NLRs/clusters and store in a dataframe.
allNLRs <- as.data.frame(read_xlsx("C:\\Users\\jexy2\\OneDrive\\Documents\\PhD\\Arabidopsis NLRs.xlsx", sheet = 2))

library(dplyr)
library(stringr)


# Create a hash to store ReMap data for each NLR/cluster.
NLR_hash <- hash()

for (row in 1:nrow(allNLRs)) {
  ReMapRows <- c(which(ReMap[,"start"] > allNLRs[row, "start"]-5000 & ReMap[,"end"] < allNLRs[row, "end"]+5000))
  NLR_hash[[allNLRs[row,"name"]]] <- ReMap[ReMapRows,]
}

rm(ReMapRows, allNLRs)

for(n in names(NLR_hash)) {
  # Run regex on name column, extracting each section
  # (experiment, epigenetic modification, ecotype, other info)
  NLR_hash[[n]][c("exp.", "epiMod", "ecotype", "info")] <- str_match(NLR_hash[[n]][,"name"], "^([0-9a-zA-Z]+)\\.([0-9a-zA-Z-]+)\\.([0-9a-zA-Z-]+)[_\\.](.*)$")[,-1]
  
  # Filter epiMod column, excluding unwanted modifications
  NLR_hash[[n]] <- NLR_hash[[n]][!NLR_hash[[n]]$epiMod %in% c("H3", "HTR12", "H2A", "H2B", "H3T3ph", "H1", "H2AV") & !NLR_hash[[n]]$ecotype %in% c("C24", "undef", "Col-x-Ler", "Ler-x-Col", "Col-x-C24"),]
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
                                   !grepl("undef_seedling_10d-h3-1kd-2",NLR_hash[[n]]$info) & !grepl("seedling_3d-wt-ehylene",NLR_hash[[n]]$info),] 
  
  # Filter info column, checking written plant ages and removing BAD AGES
  
  # For each row
  for (row in nrow(NLR_hash[[n]]):1) {
    # Find the age if it exists (in the format 1w / 8d / 10h)
    matches <- str_match(NLR_hash[[n]][row, "info"], "_([0-9]+)([dwh])")
    
    if (row %% 1000 == 0) {
      print(row)
    }
    
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


# Create hash for root data.
RootNLRs <- hash()

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
WTConditions <- allConditions[c(1,2,4:7,13,14,16:18,21,23:25,28,31,32,36,44,48,49)]
mutantConditions <- allConditions[-c(1,2,4:7,13,14,16:18,21,23:25,28,31,32,36,44,48,49)]
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

# Create a dictionary (hash containing hashes) with dataframes for each epiMod for each NLR.
ColWTLeafData <- hash()

for (n in names(ColWTLeafNLRs)) {
  epiMods <- unique(ColWTLeafNLRs[[n]]$epiMod)
  
  modHash <- hash()
  for (mod in epiMods) {
    modHash[[mod]] <- ColWTLeafNLRs[[n]][ColWTLeafNLRs[[n]]$epiMod==mod,]
  }
  ColWTLeafData[[n]] <- modHash
}



# Merge start and end coordinates columns to create a ranges column.
WTonly_Col$ranges = paste(WTonly_Col$start,"-",WTonly_Col$end, sep = "")


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

library(hash)
library(sets)

modData <- hash()
modOverlaps <- hash()

# For each epigenetic modification name
for (mod in unique(WTonly_Col$epiMod)) {
  
  # Grab all entries for that modification and store in the hash.
  modDF <- WTonly_Col[WTonly_Col$epiMod==mod,]
  modData[[mod]] <- modDF
 
  # Generate overlapSets as a list of single-item sets
  # eg, [ {1}, {2}, {3}, {4}, {5}, {6} ]
  overlapSets <- list()
  for (l in 1:nrow(modDF)) {
    overlapSets <- append(overlapSets, list(set(as.numeric(l))))
  }
  
  # For each gene co-ordinate comparison [k, l] in modDF
  for (k in 1:nrow(modDF)) {
    for (l in 1:k) {
      
      # If the co-ordinate ranges overlap
      if (overlapsFunction(modDF[k, "start"], modDF[k, "end"], 
                           modDF[l, "start"], modDF[l, "end"])==TRUE) {
        
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
  modOverlaps[[mod]] <- overlapSets
}

# Find the maximum range for the overlapping epigenetic modifications.
for (m in unique(WTonly_Col$epiMod)) {
  for (n in 1:length(modOverlaps[[m]])) {
    modStart <- c()
    modEnd <- c()
    
    for (o in modOverlaps[[m]][[n]]) {
      modStart <- append(modStart, modData[[m]][as.numeric(o), "start"])
      modEnd <- append(modEnd, modData[[m]][as.numeric(o), "end"])
    }
    modOverlaps[[m]][n] <- paste(min(modStart), max(modEnd), sep = "-")
  }
}

# Create dataframes with the information needed in the bed file.
for (p in unique(WTonly_Col$epiMod)) {
  df <- data.frame(seqname = numeric(),
                   range = character(),
                   strand = character(),
                   name = character(),
                   colour = character())
  
  for (q in 1:length(modOverlaps[[p]])) {
    df <- rbind(df, data.frame(seqname = 4,
                               range = modOverlaps[[p]][[q]],
                               strand = "*",
                               name = p,
                               colour = modData[[p]][1,"itemRgb"]))
  }
  modData[[p]] <- df
}

# Combine the dataframes for each epigenetic modification into one dataframe.
modBed <- data.frame(seqname = numeric(),
                     range = character(),
                     strand = character(),
                     name = character(),
                     colour = character())

for (g in unique(WTonly_Col$epiMod)) {
  modBed <- rbind(modBed, modData[[g]])
}

# Create bed file.
modBed <- GRanges(
  seqnames=Rle("chr5",nrow(modBed)),
  ranges=IRanges(modBed$range),
  name=modBed$name,
  itemRgb=modBed$colour)

# Export bed file.
rtracklayer::export.bed(modBed, "~/WTonlychr5.bed")


mutantsOnlyConditions <- allConditions[c(1,3,6,9,13,14,16,22,23,24,27,30,34,35,36,38,40,43,44,45,46,47,48,49,51,52,53,54,55,56,57)]

mutantData <- hash()

for (m in mutantsOnlyConditions) {
  mutantDF <- ReMapRPP5_Col[ReMapRPP5_Col$info %in% m,]
  
  # Merge start and end coordinates columns to create a ranges column.
  mutantDF$ranges = paste(mutantDF$start,"-",mutantDF$end, sep = "")
  
  # Store each mutant dataframe in the mutantData hash.
  mutantData[[m]] <- mutantDF
}

for (g in mutantsOnlyConditions) {
    
    # Create bed file.
    mutantBed <- GRanges(
      seqnames=Rle("chr4",nrow(mutantData[[g]])),
      ranges=IRanges(mutantData[[g]]$ranges),
      name=mutantData[[g]]$epiMod,
      itemRgb=mutantData[[g]]$itemRgb,
      info=mutantData[[g]]$info)
    
    # Export bed file.
    rtracklayer::export.bed(mutantBed, paste("~/",g, ".bed", sep = ""))
}


methylationBed <- modBed[modBed$name %in% c("H3K27me3", "H3K4me3", "H3K36me3", "H3K4me1", "H3K9me2"),]
rtracklayer::export.bed(methylationBed, "~/methylationChr5.bed")

acetylationBed <- modBed[modBed$name %in% c("H3K9ac", "H3K14ac", "H3K36ac", "H3K56ac"),]
rtracklayer::export.bed(acetylationBed, "~/acetylationChr5.bed")



