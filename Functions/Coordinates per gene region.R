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
  
  
  # Get the coordinates for the promotor regions.
  ATpromotors500 <- promoters(TxDb.Athaliana.BioMart.plantsmart28, upstream=500, downstream=0, use.names = TRUE)
  ATpromotors1000 <- promoters(TxDb.Athaliana.BioMart.plantsmart28, upstream=1000, downstream=0, use.names = TRUE)
  
  # Remove duplicate genes (different versions).
  ATpromotors500 <- ATpromotors500[-c(which(ATpromotors500$tx_name == str_match(ATpromotors500$tx_name, "^([0-9a-zA-Z]+)[.][2-9]+|^([0-9a-zA-Z]+)[.][1][0]$")[,1])),]
  ATpromotors1000 <- ATpromotors1000[-c(which(ATpromotors1000$tx_name == str_match(ATpromotors1000$tx_name, "^([0-9a-zA-Z]+)[.][2-9]+|^([0-9a-zA-Z]+)[.][1][0]$")[,1])),]
  
  promotor500 <- data.frame(seqnames = numeric(),
                            start = numeric(),
                            end = numeric(),
                            width = numeric(),
                            strand = factor(),
                            tx_id = numeric(),
                            tx_name = character())
  
  promotor1000 <- promotor500
  
  if (nrow(dataToUse) >= 1) {
    for (gene in dataToUse$Gene) {
      promotor500 <- rbind(promotor500, as.data.frame(ATpromotors500[grepl(gene,ATpromotors500$tx_name),]))
      promotor1000 <- rbind(promotor1000, as.data.frame(ATpromotors1000[grepl(gene,ATpromotors1000$tx_name),]))
    }
    # Alter coordinates of promotor1000 to be only 500bp upstream of promotor500.
    for (row in 1:nrow(promotor1000)) {
      if (promotor1000[row, "strand"]=="+") {
        promotor1000[row, "end"] <- promotor1000[row, "start"]+500
      }
      else if (promotor1000[row, "strand"]=="-") {
        promotor1000[row, "start"] <- promotor1000[row, "end"]-500
      }
    }
    
    # Create a ranges column by merging the start and end columns.
    promotor500$ranges <- paste(promotor500$start,"-",promotor500$end, sep = "")
    promotor1000$ranges <- paste(promotor1000$start,"-",promotor1000$end, sep = "")
    
    
    # Add a new column for the gene name, removing ".1" from the end.
    promotor500$Gene <- str_match(promotor500$tx_name, "^([0-9a-zA-Z]+)([.])([1])$")[,2]
    promotor1000$Gene <- str_match(promotor1000$tx_name, "^([0-9a-zA-Z]+)([.])([1])$")[,2]
  }
  
  
 
  # Get the coordinates for the upstream intergenic regions.
  usCoordinates <- c()
  
  if (nrow(dataToUse) >= 1) {
    for (gene in dataToUse$Gene) {
      currentGene <- which(dataToUse$Gene==gene)
      
      if (dataToUse[currentGene, "strand"]=="+") {
        previousGene <- currentGene - 1
        
        if (previousGene > 0 & previousGene < nrow(dataToUse)) {
          if (as.numeric(dataToUse[currentGene, "seqnames"])==as.numeric(dataToUse[previousGene, "seqnames"]) & 
              dataToUse[previousGene, "strand"]=="+") {
            
            distance <- (dataToUse[currentGene, "start"] - 1001) - (dataToUse[previousGene, "end"] + 201)
            
            if (distance > 0) {
              usCoordinates <- append(usCoordinates, paste(dataToUse[previousGene, "end"] + 201, "-", dataToUse[previousGene, "end"] + 201 + distance, sep = "")) 
            } 
            else usCoordinates <- append(usCoordinates, NA)
          }
          
          else if (as.numeric(dataToUse[currentGene, "seqnames"])==as.numeric(dataToUse[previousGene, "seqnames"]) & 
                   dataToUse[previousGene, "strand"]=="-") {
            
            distance <- (dataToUse[currentGene, "start"] - 1001) - (dataToUse[previousGene, "end"] + 1001)
            
            if (distance > 0) {
              usCoordinates <- append(usCoordinates, paste(dataToUse[previousGene, "end"] + 1001, "-", dataToUse[previousGene, "end"] + 1001 + distance, sep = "")) 
            } 
            else usCoordinates <- append(usCoordinates, NA)
            
          } 
        }
        else usCoordinates <- append(usCoordinates, NA)
      }
      
      else if (dataToUse[currentGene, "strand"]=="-") {
        previousGene <- currentGene + 1
        
        if (previousGene > 0 & previousGene < nrow(dataToUse)) {
          if (as.numeric(dataToUse[currentGene, "seqnames"])==as.numeric(dataToUse[previousGene, "seqnames"]) &
              dataToUse[previousGene, "strand"]=="+") {
            
            distance <- (dataToUse[previousGene, "start"] - 1001) - (dataToUse[currentGene, "end"] + 1001)
            
            if (distance > 0) {
              usCoordinates <- append(usCoordinates, paste(dataToUse[previousGene, "start"] - 1001 - distance, "-", dataToUse[previousGene, "start"] - 1001, sep = "")) 
            } 
            else usCoordinates <- append(usCoordinates, NA)
            
          }
          
          else if (as.numeric(dataToUse[currentGene, "seqnames"])==as.numeric(dataToUse[previousGene, "seqnames"]) & 
                   dataToUse[previousGene, "strand"]=="-") {
            
            distance <- (dataToUse[previousGene, "start"] - 201) - (dataToUse[currentGene, "end"] + 1001)
            if (distance > 0) {
              usCoordinates <- append(usCoordinates, paste(dataToUse[previousGene, "start"] - 201 - distance, "-", dataToUse[previousGene, "start"] - 201, sep = "")) 
            } 
            else usCoordinates <- append(usCoordinates, NA)
            
          }
        }
        else usCoordinates <- append(usCoordinates, NA)
      } 
    }
  }
  
  
  upstreamIntergenic <- dataToUse[which(dataToUse$Gene %in% dataToUse$Gene),]
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
      currentGene <- which(dataToUse$Gene==gene)
      
      if (dataToUse[currentGene, "strand"]=="-") {
        nextGene <- currentGene - 1
        
        if (as.numeric(dataToUse[currentGene, "seqnames"])==as.numeric(dataToUse[nextGene, "seqnames"])) {
          if (dataToUse[nextGene, "strand"]=="-") {
            distance <- (dataToUse[currentGene, "start"] - 201) - (dataToUse[nextGene, "end"] + 1001)
            
            if (distance > 0) {
              dsCoordinates <- append(dsCoordinates, paste(dataToUse[nextGene, "end"] + 1001, "-", dataToUse[nextGene, "end"] + 1001 + distance, sep = "")) 
            } 
            else dsCoordinates <- append(dsCoordinates, NA) 
          }
          
          else if (dataToUse[nextGene, "strand"]=="+") {
            distance <- (dataToUse[currentGene, "start"] - 201) - (dataToUse[nextGene, "end"] + 201)
            
            if (distance > 0) {
              dsCoordinates <- append(dsCoordinates, paste(dataToUse[nextGene, "end"] + 201, "-", dataToUse[nextGene, "end"] + 201 + distance, sep = "")) 
            } 
            else dsCoordinates <- append(dsCoordinates, NA)
            
          }
        }
        else dsCoordinates <- append(dsCoordinates, NA)
      }
      
      else if (dataToUse[currentGene, "strand"]=="+") {
        nextGene <- currentGene + 1
        
        if (as.numeric(dataToUse[currentGene, "seqnames"])==as.numeric(dataToUse[nextGene, "seqnames"])) {
          if (dataToUse[nextGene, "strand"]=="+") {
            distance <- (dataToUse[nextGene, "start"] - 1001) - (dataToUse[currentGene, "end"] + 201)
            
            if (distance > 0) {
              dsCoordinates <- append(dsCoordinates, paste(dataToUse[nextGene, "start"] - 1001 - distance, "-", dataToUse[nextGene, "start"] - 1001, sep = "")) 
            } 
            else dsCoordinates <- append(dsCoordinates, NA) 
            
          }
          
          else if (dataToUse[nextGene, "strand"]=="-") {
            distance <- (dataToUse[nextGene, "start"] - 201) - (dataToUse[currentGene, "end"] + 201)
            if (distance > 0) {
              dsCoordinates <- append(dsCoordinates, paste(dataToUse[nextGene, "start"] - 201 - distance, "-", dataToUse[nextGene, "start"] - 201, sep = "")) 
            } 
            else dsCoordinates <- append(dsCoordinates, NA) 
            
          }
        }
        else dsCoordinates <- append(dsCoordinates, NA)
      } 
    }
  }
  
  downstreamIntergenic <- dataToUse[which(dataToUse$Gene %in% dataToUse$Gene),]
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