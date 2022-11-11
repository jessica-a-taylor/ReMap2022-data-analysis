# Import ReMap2022 data.
ReMap <- rtracklayer::import.bed("C:\\Users\\jexy2\\OneDrive\\Documents\\PhD\\remap2022_histone_all_macs2_TAIR10_v1_0.bed.gz")

# Convert to a dataframe and define column names.
ReMap <- as.data.frame(ReMap, colnames = c("seqnames", "start", "end", "width",
                                           "strand", "name", "score", "itemRgb",
                                           "thick.start", "thick.end", "thick.width"))

# Remove unwanted columns.
ReMap <- ReMap[,-c(9:11)]

# Run regex on name column, extracting each section
# (experiment, epigenetic modification, ecotype, other info)
ReMap[c("exp.", "epiMod", "ecotype", "info")] <- str_match(ReMap[,"name"], "^([0-9a-zA-Z]+)\\.([0-9a-zA-Z-]+)\\.([0-9a-zA-Z-]+)[_\\.](.*)$")[,-1]


# Filter for conditions of interest.
conditions <- c("seedling_14d-wt","aerial-part_6d-wt","undef_seedling_10d-wt","rosette-leaves_10d-wt","seedling_12d-wt",
                "seedling_7d-wt","roots_14d-wt","seedling_14d-mock","aerial-part_12-13d-wt","roots_5d-wt","leaves_3w",
                "seedling_10d","roots_7d-wt","leaves_3w-wt","rosette-leaves","whole-plant_2w-wt","leaves_20d")

ReMap <- ReMap[c(which(ReMap$info %in% conditions)),]

ReMap <- ReMap[,-c(6,7,9)]

# Save filtered ReMap data.
write.csv(ReMap, file = "Filtered ReMap data.csv")
