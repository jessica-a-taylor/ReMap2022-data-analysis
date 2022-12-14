library(TxDb.Athaliana.BioMart.plantsmart28)
library(readr)
library(stringr)
library(hash)
library(readxl)

source("Functions\\Overlaps functions.R")

# Import all Arabidopsis genes.
Atgenes <- as.data.frame(transcriptsBy(TxDb.Athaliana.BioMart.plantsmart28, by="gene"))
colnames(Atgenes)[2] <- "Gene"

# Remove duplicate genes (different versions).
Atgenes <- Atgenes[-c(which(Atgenes$tx_name == str_match(Atgenes$tx_name, "^([0-9a-zA-Z]+)([.])([2-9])$")[,1])),]


# Create separate datasets for the euchromatic and heterochromatic regions.
pericentromericgeneRegions <- data.frame(Chromosome = c(1:5),
                                         Start = c("11500000", "1100000", "10300000", "1500000", "9000000"),
                                         End = c("17700000", "7200000", "17300000", "6300000", "16000000"))

euchromaticgeneRegions <- data.frame()
heterochromaticRegions <- data.frame()

for (row in 1:nrow(pericentromericgeneRegions)) {
  euchromaticgeneRegions <- rbind(euchromaticgeneRegions,
                                  Atgenes[c(which(Atgenes$seqnames==row & 
                                                    betweenFunction(Atgenes$start, as.numeric(pericentromericgeneRegions[row, "Start"]), as.numeric(pericentromericgeneRegions[row, "End"]))==FALSE &
                                                    betweenFunction(Atgenes$end, as.numeric(pericentromericgeneRegions[row, "Start"]), as.numeric(pericentromericgeneRegions[row, "End"]))==FALSE)),])
  
  heterochromaticRegions <- rbind(heterochromaticRegions,
                                  Atgenes[c(which(Atgenes$seqnames==row & 
                                                    betweenFunction(Atgenes$start, as.numeric(pericentromericgeneRegions[row, "Start"]), as.numeric(pericentromericgeneRegions[row, "End"]))==TRUE &
                                                    betweenFunction(Atgenes$end, as.numeric(pericentromericgeneRegions[row, "Start"]), as.numeric(pericentromericgeneRegions[row, "End"]))==TRUE)),])
}

rm(pericentromericgeneRegions)

genomicRegions <- hash(euchromatic = euchromaticgeneRegions, 
                       heterochromatic = heterochromaticRegions)

for (region in names(genomicRegions)) {
  # Remove duplicate genes.
  genomicRegions[[region]] <- genomicRegions[[region]][c(which(is.na(str_match(genomicRegions[[region]]$tx_name, "^[0-9a-zA-Z]+[.][2-9]+|^[0-9a-zA-Z]+[.][1][0]$")[,1]))),]
  
  genomicRegions[[region]] <- genomicRegions[[region]][,-c(1,8,9)]
  genomicRegions[[region]]$ranges <- paste(genomicRegions[[region]]$start,"-",genomicRegions[[region]]$end, sep = "")
  
  write.csv(genomicRegions[[region]], file = paste("Data\\", region, ".csv", sep=""))
}

# Remove TEs from the euchromaticgeneRegions dataframe.
transposableElements <- as.data.frame(read_xlsx("Data\\Arabidopsis TE genes.xlsx"))
write.csv(genomicRegions[["euchromatic"]][-c(which(genomicRegions[["euchromatic"]]$Gene %in% transposableElements$Locus)),], 
          file = "Data\\euchromaticWithoutTEs.csv")