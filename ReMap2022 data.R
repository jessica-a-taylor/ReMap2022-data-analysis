install.packages(c("readxl", "karyoploteR", "rtracklayer"))
library(readxl)
library(karyoploteR)
library(rtracklayer)

# Import chromatin states dataset.
chromStates <- as.data.frame(read_xlsx("C:\\Users\\jexy2\\OneDrive\\Documents\\PhD\\RPP5 ChromStates.xlsx"))

# Merge start and end coordinates columns to create a ranges column.
chromStates$ranges = paste(chromStates$start,"-",chromStates$end, sep = "")
chromStates

# Create a bed file for chromatin states dataset.
chromStatesBed <- GRanges(
  seqnames=Rle("chr4",nrow(chromStates)),
  ranges=IRanges(chromStates$ranges),
  name=chromStates$state)

# Export bed file.
rtracklayer::export.bed(chromStatesBed, "~/CS.bed")

# Import RPP5 gene coordinates dataset.
RPP5 <- as.data.frame(read_xlsx("C:\\Users\\jexy2\\OneDrive\\Documents\\PhD\\RPP5.xlsx"))

# Merge start and end coordinates columns to create a ranges column.
RPP5$ranges = paste(RPP5$start,"-",RPP5$end, sep = "")
RPP5

# Create a bed file for the RPP5 gene coordinates dataset.
RPP5Bed <- GRanges(
  seqnames=Rle("chr4",nrow(RPP5)),
  ranges=IRanges(RPP5$ranges),
  name=RPP5$name)

# Export bed file.
rtracklayer::export.bed(RPP5Bed, "~/RPP5.bed")

# Import ReMap2022 ChIP dataset for SNC1, RPP4 and RPP5 only.
epiMods <- as.data.frame(read_xlsx("C:\\Users\\jexy2\\OneDrive\\Documents\\PhD\\ReMap2022 data.xlsx"))

# Merge start and end coordinates columns to create a ranges column.
epiMods$ranges = paste(epiMods$Start,"-",epiMods$End, sep = "")
epiMods

# Create a bed file for the SNC1, RPP4 and RPP5 data.
epiModsBed <- GRanges(
  seqnames=Rle("chr4",nrow(epiMods)),
  ranges=IRanges(epiMods$ranges),
  name=epiMods$`Epigenetic modifications (ReMap2022)`)

# Export bed file.
rtracklayer::export.bed(epiMods_foo, "~/epiMods.bed")


# Import all ReMap2022 ChIP data (bed file).
ReMap <- rtracklayer::import.bed("C:\\Users\\jexy2\\OneDrive\\Documents\\PhD\\remap2022_histone_all_macs2_TAIR10_v1_0.bed.gz")

# Convert to a dataframe and define column names.
ReMap <- as.data.frame(ReMap, colnames = c("seqnames", "start", "end", "width",
                                           "strand", "name", "score", "itemRgb",
                                           "thick.start", "thick.end", "thick.width"))

# Filter for chromosome 4.
ReMapRPP5 <- ReMap[ReMap$seqnames==4,]

# Remove original dataset.
rm(ReMap)

# Filter for the coordinates of the RPP5 gene cluster.
ReMapRPP5 <- ReMapRPP5[ReMapRPP5$start>9480000,]
ReMapRPP5 <- ReMapRPP5[ReMapRPP5$end<9570000,]

install.packages(c("dplyr", "stringr"))
library(dplyr)
library(stringr)

# Run regex on name column, extracting each section
# (experiment, epigenetic modification, ecotype, other info)
ReMapRPP5[c("exp.", "epiMod", "ecotype", "info")] <- str_match(ReMapRPP5$name, "^([0-9a-zA-Z]+)\\.([0-9a-zA-Z-]+)\\.([0-9a-zA-Z-]+)[_\\.](.*)$")[,-1]

# Filter epiMod column, excluding unwanted modifications
ReMapRPP5 <- ReMapRPP5[!ReMapRPP5$epiMod %in% c("H3", "HTR12", "H2A", "H2B", "H3T3ph", "H1") & !ReMapRPP5$ecotype %in% c("C24", "undef", "Col-x-Ler", "Ler-x-Col", "Col-x-C24"),]

# Filter info column, excluding unwanted conditions (too old, too young, wrong part of plant, etc)
ReMapRPP5 <- ReMapRPP5[!grepl("mutant",ReMapRPP5$info) & !grepl("mature",ReMapRPP5$info) & !grepl("senescent",ReMapRPP5$info) & 
                         !grepl("inflorescence",ReMapRPP5$info) & !grepl("drought",ReMapRPP5$info) & !grepl("old",ReMapRPP5$info) & 
                         !grepl("min",ReMapRPP5$info) & !grepl("endosperm",ReMapRPP5$info) & !grepl("-se-",ReMapRPP5$info) &
                         !grepl("-TSA-",ReMapRPP5$info) & !grepl("-GSNO-",ReMapRPP5$info) & !grepl("flg22",ReMapRPP5$info) &
                         !grepl("GSE79354",ReMapRPP5$exp.) & !grepl("transgenic",ReMapRPP5$info) & !grepl("GSH",ReMapRPP5$info) &
                         !grepl("GSE93223",ReMapRPP5$exp.) & !grepl("-acc1",ReMapRPP5$info) & !grepl("GSE93875",ReMapRPP5$exp.) &
                         !grepl("leaves_3w-K36M-homoz",ReMapRPP5$info) & !grepl("undef_seedling_10d-h3-1kd-1",ReMapRPP5$info) &
                         !grepl("undef_seedling_10d-h3-1kd-2",ReMapRPP5$info) & !grepl("seedling_3d-wt-ehylene",ReMapRPP5$info) &
                         !grepl("GSE77394",ReMapRPP5$exp.) & !grepl("GSE100965",ReMapRPP5$exp.),] 

# Filter info column, checking written plant ages and removing BAD AGES

# For each row
for (row in nrow(ReMapRPP5):1) {
  # Find the age if it exists (in the format 1w / 8d / 10h)
  matches <- str_match(ReMapRPP5[row, "info"], "_([0-9]+)([dwh])")
  
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
      ReMapRPP5 <- ReMapRPP5[-row,]
    }
  }
}

# Tidy up 🧹
rm(matches, badAge, maxWeeks, row, timeValue)

# Create separate dataframes for each ecotype.
ecotypes <- unique(ReMapRPP5$ecotype)
ReMapRPP5_Col <- ReMapRPP5[ReMapRPP5$ecotype=="Col-0",]
ReMapRPP5_Ler <- ReMapRPP5[ReMapRPP5$ecotype=="Ler",]

# Create lists of WT and mutant conditions from the info column.
allConditions <- unique(ReMapRPP5$info)
WTonlyConditions <- allConditions[-c(1,3,6,9,13,14,16,22,23,24,27,30,34,35,36,38,40,43,44,45,46,47,48,49,51,52,53,54,55,56,57)]

# Create a ReMapRPP5 dataset for Col-0 WT only.
WTonly_Col <- ReMapRPP5_Col[ReMapRPP5_Col$info %in% c(WTonlyConditions),]

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

install.packages(c("hash", "sets"))
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
  seqnames=Rle("chr4",nrow(modBed)),
  ranges=IRanges(modBed$range),
  name=modBed$name,
  itemRgb=modBed$colour)

# Export bed file.
rtracklayer::export.bed(modBed, "~/WTonlyRPP5.bed")


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