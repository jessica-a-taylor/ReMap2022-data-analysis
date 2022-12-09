source("Functions\\Overlaps functions.R")

region = c("UpstreamIntergenic", "Promotor1000", "Promotor500",
          "Gene20", "Gene40", "Gene60", "Gene80", "Gene100", 
          "Downstream", "DownstreamIntergenic")

# Functions for determining the % R-genes with a chromatin mark in each gene region (frequency)
# and the proportion of each gene region with that mark.

modFrequenciesFunction <- function (geneRegions, allOverlaps, epiMods) {
  modFrequencyPerRegion <- hash()
 
  for (r in names(geneRegions)) {
    modFrequenciesDF <- data.frame(Region= character(), 
                                   Modification = character(),
                                   Frequency = numeric())
    
    if (length(names(allOverlaps)) >= 1) {
      for (mod in epiMods) {
        geneList <- c()
        
        for (n in names(allOverlaps)) {
          modPresent <- FALSE
          
          if (nrow(allOverlaps[[n]][[mod]]) >= 1 & n %in% geneRegions[[r]]$Gene == TRUE) {
            
            for (row in 1:nrow(allOverlaps[[n]][[mod]])) {
              if (overlapsFunction(as.numeric(allOverlaps[[n]][[mod]][row, "start"]), as.numeric(allOverlaps[[n]][[mod]][row, "end"]),
                                   as.numeric(geneRegions[[r]][geneRegions[[r]]$Gene==n,]$start), as.numeric(geneRegions[[r]][geneRegions[[r]]$Gene==n,]$end))==TRUE) {
                modPresent <- TRUE
              }
              else modPresent <- modPresent
            }
            if (modPresent == TRUE) {
              geneList <- append(geneList, n)
            }
            else geneList <- geneList
          }
        }
        modFrequenciesDF <- rbind(modFrequenciesDF, data.frame(Region = r,
                                                               Modification = mod,
                                                               Frequency = length(geneList)/length(names(allOverlaps))*100))
      }
    } else modFrequenciesDF <- modFrequenciesDF
    
    modFrequencyPerRegion[[r]] <- modFrequenciesDF
  }
  
  # Collect all hashes into a single dataframe.
  DF <- data.frame(Region = character(),
                   Modification = character(),
                   Measure = numeric(),
                   n = numeric())
  
  for (r in region) {
    DF <- rbind(DF, modFrequencyPerRegion[[r]])
  }
  return(DF)
}


modProportionsFunction <- function (geneRegions, allOverlaps, epiMods) {
  modProportionPerRegion <- hash()
  
  for (r in names(geneRegions)) {
    modProportionDF <- data.frame(Region= character(), 
                                  Modification = character(),
                                  Proportion = numeric())
    
    if (length(names(allOverlaps)) >= 1) {
      for (mod in epiMods) {

        for (n in names(allOverlaps)) {
          modPresent <- FALSE
          modOverlaps <- c()
          
          if (nrow(allOverlaps[[n]][[mod]]) >= 1 & n %in% geneRegions[[r]]$Gene == TRUE) {
            
            for (row in 1:nrow(allOverlaps[[n]][[mod]])) {
              modOverlaps <- append(modOverlaps, newOverlapsFunction(as.numeric(allOverlaps[[n]][[mod]][row, "start"]), as.numeric(allOverlaps[[n]][[mod]][row, "end"]),
                                                                     as.numeric(geneRegions[[r]][geneRegions[[r]]$Gene==n,]$start), as.numeric(geneRegions[[r]][geneRegions[[r]]$Gene==n,]$end)))
            }
            modProportionDF <- rbind(modProportionDF, data.frame(Region = r,
                                                                 Modification = mod,
                                                                 Proportion = sum(modOverlaps)/geneRegions[[r]][geneRegions[[r]]$Gene==n,]$width))
          }
          else modProportionDF <- rbind(modProportionDF, data.frame(Region = r,
                                                                    Modification = mod,
                                                                    Proportion = 0))
        }
      }
    } else modProportionDF <- modProportionDF
    
    modProportionPerRegion[[r]] <- modProportionDF
  }
  # Collect all hashes into a single dataframe.
  DF <- data.frame(Region = character(),
                   Modification = character(),
                   Measure = numeric())
  
  for (r in region) {
    DF <- rbind(DF, modProportionPerRegion[[r]])
  }
  return(DF)
}


# Function to add a column to the dataframe with the numbers on the x axis that will correspond with each gene region.
geneRegionAxisLocations <- function(dataToUse, geneRegions) {
  grouping <- c(seq(from = -60, to = -20, by = 20),seq(from = 20, to = 140, by = 20))
  axisGroup <- c()
  
  for (c in 1:length(names(geneRegions))) {
    axisGroup <- append(axisGroup, rep(grouping[c], times = nrow(dataToUse[dataToUse$Region == names(geneRegions)[c],])))
  }
  dataToUse <- cbind(dataToUse, axisGroup)
 
  return(dataToUse) 
}


# Function to add a column for expression level.
expressionColumn <- function(dataToUse, level) {
  if (nrow(dataToUse) >= 1) {
    dataToUse <- cbind(dataToUse, data.frame(Expression = rep(level, times = nrow(dataToUse))))
  }
  else dataToUse <- dataToUse
  return(dataToUse)
}