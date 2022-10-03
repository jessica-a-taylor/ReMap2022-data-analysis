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
ReMapRPP5 <- ReMapRPP5[ReMapRPP5$start>9488000,]
ReMapRPP5 <- ReMapRPP5[ReMapRPP5$end<9566000,]

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