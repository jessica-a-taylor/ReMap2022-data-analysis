source("Functions\\Overlaps functions.R")

# Function for creating a hash with the ReMap data in a particular tissue for the current set of genes. 

ReMapPerGene <- function(dataToUse, tissueForAnalysis) {
  allModifications <- hash()
  
  if (nrow(dataToUse) >= 1) {
    for (row in 1:nrow(dataToUse)) {
      # Select rows that are within the range of each gene and on the same chromosome.
      ReMapRows <- c(which(ReMap[,"start"] > dataToUse[row, "start"]-5000 & ReMap[,"end"] < dataToUse[row, "end"]+5000 & ReMap[,"seqnames"] == as.numeric(dataToUse[row, "seqnames"])))
      allModifications[[dataToUse[row,"Gene"]]] <- ReMap[ReMapRows,]
    }
  }
  
  # Create hash for leaf and root data.
  rootGenes <- hash()
  leafGenes <- hash()
  seedlingGenes <- hash()
  
  
  for(n in names(allModifications)) {
    rootGenes[[n]] <- allModifications[[n]][c(which(allModifications[[n]]$info %in% c("roots_14d-wt", "roots_5d-wt","roots_7d-wt"))),]
    
    leafGenes[[n]] <- allModifications[[n]][c(which(allModifications[[n]]$info %in% c("leaves_3w", "leaves_3w-wt", "rosette-leaves", "leaves_20d"))),]
    
    seedlingGenes[[n]] <- allModifications[[n]][c(which(allModifications[[n]]$info %in% c("seedling_14d-wt", "aerial-part_6d-wt", "undef_seedling_10d-wt",
                                                                                          "rosette-leaves_10d-wt", "seedling_12d-wt", "seedling_7d-wt",
                                                                                          "seedling_14d-mock", "aerial-part_12-13d-wt", "seedling_10d",
                                                                                          "whole-plant_2w-wt"))),]
  }
  
  # Delete allModifications.
  rm(allModifications)
  
  if (tissueForAnalysis == "leafGenes") {
    ColGenes <- hash()
    LerGenes <- hash()
  
    for (n in names(leafGenes)) {
      ColGenes[[n]] <- leafGenes[[n]][grepl("Col-0",leafGenes[[n]]$ecotype),]
      LerGenes[[n]] <- leafGenes[[n]][grepl("Ler",leafGenes[[n]]$ecotype),]
    }
  return(ColGenes)
  } else if (tissueForAnalysis == "rootGenes") {
    ColGenes <- hash()
    LerGenes <- hash()
    
    for (n in names(rootGenes)) {
      ColGenes[[n]] <- rootGenes[[n]][grepl("Col-0",rootGenes[[n]]$ecotype),]
      LerGenes[[n]] <- rootGenes[[n]][grepl("Ler",rootGenes[[n]]$ecotype),]
    }
    return(ColGenes)
  }else if (tissueForAnalysis == "seedlingGenes") {
    ColGenes <- hash()
    LerGenes <- hash()
    
    for (n in names(seedlingGenes)) {
      ColGenes[[n]] <- seedlingGenes[[n]][grepl("Col-0",seedlingGenes[[n]]$ecotype),]
      LerGenes[[n]] <- seedlingGenes[[n]][grepl("Ler",seedlingGenes[[n]]$ecotype),]
    }
    return(ColGenes)
  }
  rm(leafGenes)
}


# Function for creating a hash with the occurrences of a chromatin modification for each gene.

modificationOccurrences <- function(allModifications) {
  geneModifications <- hash()
  
  for (n in names(allModifications)) {
    modHash <- hash()
    
    for (mod in epiMods) {
      modHash[[mod]] <- allModifications[[n]][allModifications[[n]]$epiMod==mod,]
    }
    geneModifications[[n]] <- modHash
  }

  # Merge start and end coordinates columns to create a ranges column.
  source("Functions\\Get range - merge gene coordinates.R")
  
  for (n in names(geneModifications)) {
    for (mod in epiMods) {
      
      if (nrow(geneModifications[[n]][[mod]]) >= 1) {
        dataToUse <- geneModifications[[n]][[mod]]
        
        geneModifications[[n]][[mod]]$ranges <- mergeCoordinates(dataToUse)
      }
      else next
    }
  }
  return(geneModifications)
}


# Function for creating a hash with the overlapping occurrences of each modification merged together.

mergeOverlappingModifications <- function(geneModifications) {
  allOverlaps <- hash()
  
  # For each epigenetic modification name
  for (n in names(geneModifications)) {
    modOverlaps <- hash()
    
    for (mod in epiMods) {
      
      # Generate overlapSets as a list of single-item sets
      # eg, [ {1}, {2}, {3}, {4}, {5}, {6} ]
      overlapSets <- list()
      if (nrow(geneModifications[[n]][[mod]])>0) {
        
        for (r in 1:nrow(geneModifications[[n]][[mod]])) {
          overlapSets <- append(overlapSets, list(sets::set(as.numeric(r))))
        }
        #For each gene co-ordinate comparison [k, l]
        for (k in 1:nrow(geneModifications[[n]][[mod]])) {
          for (l in 1:k) {
            
            # If the co-ordinate ranges overlap
            if (overlapsFunction(geneModifications[[n]][[mod]][k, "start"], geneModifications[[n]][[mod]][k, "end"], 
                                 geneModifications[[n]][[mod]][l, "start"], geneModifications[[n]][[mod]][l, "end"])==TRUE) {
              
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

  
  # Find the maximum range for the overlapping epigenetic modifications.
  for (n in names(allOverlaps)) {
    for (mod in epiMods) {
      if (length(allOverlaps[[n]][[mod]])>0) {
        
        for (l in 1:length(allOverlaps[[n]][[mod]])) {
          modStart <- c()
          modEnd <- c() 
          
          for (o in allOverlaps[[n]][[mod]][l]) {
            modStart <- append(modStart, geneModifications[[n]][[mod]][as.numeric(o), "start"])
            modEnd <- append(modEnd, geneModifications[[n]][[mod]][as.numeric(o), "end"])
            
            allOverlaps[[n]][[mod]][l] <- paste(min(modStart), max(modEnd), sep = "-")
          }
        }
      }
    }
  }
  

  # Create dataframes with the information needed in the bed file.
  for (n in names(geneModifications)) {
    for (mod in epiMods) {
      df <- data.frame(seqnames = numeric(),
                       start = numeric(),
                       end = numeric(),
                       width = numeric(),
                       ranges = character(),
                       strand = factor(),
                       epiMod = character(),
                       colour = character())
      
      if (length(allOverlaps[[n]][[mod]])>0) {
        for (l in 1:length(allOverlaps[[n]][[mod]])) {
          df <- rbind(df, data.frame(seqnames = geneModifications[[n]][[mod]][1,"seqnames"],
                                     start = str_match(allOverlaps[[n]][[mod]][[l]], "^([0-9]+)-([0-9]+)$")[,2],
                                     end = str_match(allOverlaps[[n]][[mod]][[l]], "^([0-9]+)-([0-9]+)$")[,3],
                                     width = as.numeric(str_match(allOverlaps[[n]][[mod]][[l]], "^([0-9]+)-([0-9]+)$")[,3]) - as.numeric(str_match(allOverlaps[[n]][[mod]][[l]], "^([0-9]+)-([0-9]+)$")[,2]),
                                     ranges = allOverlaps[[n]][[mod]][[l]],
                                     strand = geneModifications[[n]][[mod]][1,"strand"],
                                     epiMod = mod,
                                     colour = geneModifications[[n]][[mod]][1,"itemRgb"]))
        }
      }
      allOverlaps[[n]][[mod]] <- df
    }
  }
  
  return(allOverlaps)
}

