# Function to get the coordinates of each gene region.

getGeneCoordinates <- function(dataToUse) {
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
  
  if (nrow(dataToUse) >= 1) {
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
  }
  
  for (n in names(geneChunks)) {
    geneChunks[[n]]$ranges <- geneWidth[[n]]
  }
  
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
    geneChunks[[n]]$width <- as.numeric(geneChunks[[n]]$end) - as.numeric(geneChunks[[n]]$start)
  }
  
  # Create new dataframe for the coordinates of the regions 200bp downstream of the TTS.
  
  downstreamRegion <- c()
  if (nrow(dataToUse) >= 1) {
    for (row in 1:nrow(dataToUse)) {
      if (dataToUse[row, "strand"]=="+") {
        downstreamRegion <- append(downstreamRegion, paste(dataToUse[row,"end"],"-",dataToUse[row,"end"]+200, sep = ""))
      }
      else if (dataToUse[row, "strand"]=="-") {
        downstreamRegion <- append(downstreamRegion, paste(dataToUse[row,"start"]-200,"-", dataToUse[row,"start"], sep = ""))
      }
    }
  }
  
  
  downstream <- dataToUse[,c(which(colnames(dataToUse) != "start" & colnames(dataToUse) != "end" &
                                     colnames(dataToUse) != "width" & colnames(dataToUse) != "ranges"))]
  
  downstream$ranges <- downstreamRegion

  downstream$start <- as.numeric(str_match(downstream$ranges, "^([0-9]+)(-)([0-9]+)$")[,2])
  downstream$end <- as.numeric(str_match(downstream$ranges, "^([0-9]+)(-)([0-9]+)$")[,4])
  downstream$width <- downstream$end - downstream$start
  
  
  # Create new dataframe for the coordinates of the 500bp promotor region.
  
  PromotorRegion <- c()
  if (nrow(dataToUse) >= 1) {
    for (row in 1:nrow(dataToUse)) {
      if (dataToUse[row, "strand"]=="+") {
        PromotorRegion <- append(PromotorRegion, paste(dataToUse[row,"start"]-500,"-",dataToUse[row,"start"], sep = ""))
      }
      else if (dataToUse[row, "strand"]=="-") {
        PromotorRegion <- append(PromotorRegion, paste(dataToUse[row,"end"],"-", dataToUse[row,"end"]+500, sep = ""))
      }
    }
  }
  
  
  Promotor500 <- dataToUse[,c(which(colnames(dataToUse) != "start" & colnames(dataToUse) != "end" &
                                     colnames(dataToUse) != "width" & colnames(dataToUse) != "ranges"))]
  
  Promotor500$ranges <- PromotorRegion
  
  Promotor500$start <- as.numeric(str_match(Promotor500$ranges, "^([0-9]+)(-)([0-9]+)$")[,2])
  Promotor500$end <- as.numeric(str_match(Promotor500$ranges, "^([0-9]+)(-)([0-9]+)$")[,4])
  Promotor500$width <- Promotor500$end - Promotor500$start
  
  
  # Create new dataframe for the coordinates of the 500bp-1000bp promotor region.
  
  PromotorRegion <- c()
  if (nrow(dataToUse) >= 1) {
    for (row in 1:nrow(dataToUse)) {
      if (dataToUse[row, "strand"]=="+") {
        PromotorRegion <- append(PromotorRegion, paste(dataToUse[row,"start"]-1000,"-",dataToUse[row,"start"], sep = ""))
      }
      else if (dataToUse[row, "strand"]=="-") {
        PromotorRegion <- append(PromotorRegion, paste(dataToUse[row,"end"],"-", dataToUse[row,"end"]+1000, sep = ""))
      }
    }
  }
  
  
  Promotor1000 <- dataToUse[,c(which(colnames(dataToUse) != "start" & colnames(dataToUse) != "end" &
                                      colnames(dataToUse) != "width" & colnames(dataToUse) != "ranges"))]
  
  Promotor1000$ranges <- PromotorRegion
  
  Promotor1000$start <- as.numeric(str_match(Promotor1000$ranges, "^([0-9]+)(-)([0-9]+)$")[,2])
  Promotor1000$end <- as.numeric(str_match(Promotor1000$ranges, "^([0-9]+)(-)([0-9]+)$")[,4])
  Promotor1000$width <- Promotor1000$end - Promotor1000$start
  
 
  # Get the coordinates for the upstream intergenic regions.
  usCoordinates <- c()
  
  if (nrow(dataToUse) >= 1) {
    for (gene in dataToUse$Gene) {
      currentGene <- which(genomicData$Gene==gene)
      
      if (genomicData[currentGene, "strand"]=="+") {
          previousGene <- currentGene - 1
          
          if (previousGene > 0 & previousGene < nrow(genomicData)) {
            if (as.numeric(genomicData[currentGene, "seqnames"])==as.numeric(genomicData[previousGene, "seqnames"]) & 
                genomicData[previousGene, "strand"]=="+") {
              
              distance <- (genomicData[currentGene, "start"] - 1001) - (genomicData[previousGene, "end"] + 201)
              
              if (distance > 0) {
                usCoordinates <- append(usCoordinates, paste(genomicData[previousGene, "end"] + 201, "-", genomicData[previousGene, "end"] + 201 + distance, sep = "")) 
              } else usCoordinates <- append(usCoordinates, NA)
            }
            
            else if (as.numeric(genomicData[currentGene, "seqnames"])==as.numeric(genomicData[previousGene, "seqnames"]) & 
                     genomicData[previousGene, "strand"]=="-") {
              
              distance <- (genomicData[currentGene, "start"] - 1001) - (genomicData[previousGene, "end"] + 1001)
              
              if (distance > 0) {
                usCoordinates <- append(usCoordinates, paste(genomicData[previousGene, "end"] + 1001, "-", genomicData[previousGene, "end"] + 1001 + distance, sep = "")) 
              } else usCoordinates <- append(usCoordinates, NA)
              
            } else usCoordinates <- append(usCoordinates, NA)
          } else usCoordinates <- append(usCoordinates, NA)
      }
      else if (genomicData[currentGene, "strand"]=="-") {
        previousGene <- currentGene + 1
        
        if (previousGene > 0 & previousGene < nrow(genomicData)) {
          if (as.numeric(genomicData[currentGene, "seqnames"])==as.numeric(genomicData[previousGene, "seqnames"]) &
              genomicData[previousGene, "strand"]=="+") {
            
            distance <- (genomicData[previousGene, "start"] - 1001) - (genomicData[currentGene, "end"] + 1001)
            
            if (distance > 0) {
              usCoordinates <- append(usCoordinates, paste(genomicData[previousGene, "start"] - 1001 - distance, "-", genomicData[previousGene, "start"] - 1001, sep = "")) 
            } else usCoordinates <- append(usCoordinates, NA)
            
          }
          
          else if (as.numeric(genomicData[currentGene, "seqnames"])==as.numeric(genomicData[previousGene, "seqnames"]) & 
                   genomicData[previousGene, "strand"]=="-") {
            
            distance <- (genomicData[previousGene, "start"] - 201) - (genomicData[currentGene, "end"] + 1001)
            if (distance > 0) {
              usCoordinates <- append(usCoordinates, paste(genomicData[previousGene, "start"] - 201 - distance, "-", genomicData[previousGene, "start"] - 201, sep = "")) 
            } else usCoordinates <- append(usCoordinates, NA)
            
          } else usCoordinates <- append(usCoordinates, NA)
        } else usCoordinates <- append(usCoordinates, NA)
      } 
    }
  }
  
  
  upstreamIntergenic <- genomicData[which(genomicData$Gene %in% dataToUse$Gene),]
  upstreamIntergenic$ranges <- usCoordinates 
  
  if (nrow(upstreamIntergenic) >= 1) {
    for (row in 1:nrow(upstreamIntergenic)) {
      upstreamIntergenic[row, "start"] <- str_match(upstreamIntergenic[row, "ranges"], "^([0-9]+)-([0-9]+)$")[,2]
      upstreamIntergenic[row, "end"] <- str_match(upstreamIntergenic[row, "ranges"], "^([0-9]+)-([0-9]+)$")[,3]
      
      if (!is.na(upstreamIntergenic[row, "ranges"])) {
        upstreamIntergenic[row, "width"] <- as.numeric(upstreamIntergenic[row, "end"]) - as.numeric(upstreamIntergenic[row, "start"])
      }
      else upstreamIntergenic[row, "width"] <- NA
    }
  }
  
  upstreamIntergenic <- upstreamIntergenic[-c(which(is.na(upstreamIntergenic$ranges))),]
  

  # Get the coordinates for the downstream intergenic regions.
  dsCoordinates <- c()
  
  if (nrow(dataToUse) >= 1) {
    for (gene in dataToUse$Gene) {
      currentGene <- which(genomicData$Gene==gene)
      
      if (genomicData[currentGene, "strand"]=="-") {
          nextGene <- currentGene - 1
        
        if (as.numeric(genomicData[currentGene, "seqnames"])==as.numeric(genomicData[nextGene, "seqnames"])) {
          if (genomicData[nextGene, "strand"]=="-") {
            distance <- (genomicData[currentGene, "start"] - 201) - (genomicData[nextGene, "end"] + 1001)
            
            if (distance > 0) {
              dsCoordinates <- append(dsCoordinates, paste(genomicData[nextGene, "end"] + 1001, "-", genomicData[nextGene, "end"] + 1001 + distance, sep = "")) 
            } else dsCoordinates <- append(dsCoordinates, NA) 
          }
          
          else if (genomicData[nextGene, "strand"]=="+") {
            distance <- (genomicData[currentGene, "start"] - 201) - (genomicData[nextGene, "end"] + 201)
            
            if (distance > 0) {
              dsCoordinates <- append(dsCoordinates, paste(genomicData[nextGene, "end"] + 201, "-", genomicData[nextGene, "end"] + 201 + distance, sep = "")) 
            } else dsCoordinates <- append(dsCoordinates, NA)
            
          } else dsCoordinates <- append(dsCoordinates, NA)

        } else dsCoordinates <- append(dsCoordinates, NA)
      }
      else if (genomicData[currentGene, "strand"]=="+") {
        nextGene <- currentGene + 1
        
        if (as.numeric(genomicData[currentGene, "seqnames"])==as.numeric(genomicData[nextGene, "seqnames"])) {
          if (genomicData[nextGene, "strand"]=="+") {
            distance <- (genomicData[nextGene, "start"] - 1001) - (genomicData[currentGene, "end"] + 201)
            
            if (distance > 0) {
              dsCoordinates <- append(dsCoordinates, paste(genomicData[nextGene, "start"] - 1001 - distance, "-", genomicData[nextGene, "start"] - 1001, sep = "")) 
            } else dsCoordinates <- append(dsCoordinates, NA) 
          }
          
          else if (genomicData[nextGene, "strand"]=="-") {
            distance <- (genomicData[nextGene, "start"] - 201) - (genomicData[currentGene, "end"] + 201)
            
            if (distance > 0) {
              dsCoordinates <- append(dsCoordinates, paste(genomicData[nextGene, "start"] - 201 - distance, "-", genomicData[nextGene, "start"] - 201, sep = "")) 
            } else dsCoordinates <- append(dsCoordinates, NA) 
            
          } else dsCoordinates <- append(dsCoordinates, NA)

        } else dsCoordinates <- append(dsCoordinates, NA)
      } 
    }
  }
  
  downstreamIntergenic <- genomicData[which(genomicData$Gene %in% dataToUse$Gene),]
  downstreamIntergenic$ranges <- dsCoordinates 
  
  if (nrow(downstreamIntergenic) >= 1) {
    for (row in 1:nrow(downstreamIntergenic)) {
      downstreamIntergenic[row, "start"] <- str_match(downstreamIntergenic[row, "ranges"], "^([0-9]+)-([0-9]+)$")[,2]
      downstreamIntergenic[row, "end"] <- str_match(downstreamIntergenic[row, "ranges"], "^([0-9]+)-([0-9]+)$")[,3]
      
      if (!is.na(downstreamIntergenic[row, "ranges"])) {
        downstreamIntergenic[row, "width"] <- as.numeric(downstreamIntergenic[row, "end"]) - as.numeric(downstreamIntergenic[row, "start"])
      }
      else downstreamIntergenic[row, "width"] <- NA
    }
  }
  
  downstreamIntergenic <- downstreamIntergenic[-c(which(is.na(downstreamIntergenic$ranges))),]

  geneRegions <- hash(UpstreamIntergenic = upstreamIntergenic, Promotor1000 = promotor1000, Promotor500 = promotor500,
                  Gene20 = geneChunks[["width20"]],  Gene40 = geneChunks[["width40"]], Gene60 = geneChunks[["width60"]], 
                  Gene80 = geneChunks[["width80"]], Gene100 = geneChunks[["width100"]], Downstream = downstream,
                  DownstreamIntergenic = downstreamIntergenic)
  
  return(geneRegions)
}