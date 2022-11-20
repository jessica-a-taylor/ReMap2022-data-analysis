# Function to get filtered expression data for each set of sample genes in each tissue. 
PlantExp <- function(dataToUse) {
  
  # Get a list of tissue types for which there is expression data.
  tissue <- list.files(path = "Data\\PlantExp data")
  
  # For each set of sample genes...
  for (test in names(dataToUse)) {

    # For each tissue type...
    for (t in tissue) {
      # Create a hash.
      sampleGenes[[paste(test, "_", t, sep = "")]] <- hash()
      
      # Create a list of files in the folder for a particular tissue type.
      filenames <- list.files(path = paste("Data\\PlantExp data\\", t, sep = ""),pattern="*.tsv")
      data <- data.frame()
      PlantExpData <- data.frame()
      
      # Merge the data from all files.
      for (file in filenames) {
        data <- rbind(data, as.data.frame(read_tsv(paste("Data\\PlantExp data\\", t, "\\", file, sep = ""), show_col_types = FALSE)))
      }
      
      # Filter for the genes in the sample gene set.
      data <- data[c(which(data$geneId %in% dataToUse[[test]]$Gene)),]
      data <- data[,-c(2,3,5)]
      
      # For each gene, calculate the ,ean expression across experiments.
      for (gene in unique(data$geneId)) {
        df <- data[data$geneId==gene,]
        
        
        PlantExpData <- rbind(PlantExpData, data.frame(geneId = gene,
                                                       FPKM = mean(df$FPKM),
                                                       tissue = t))
      }
      # Add a column to PlantExpData with the expression level of each gene.
      expressionLevel <- c()
      
      for (row in 1:nrow(PlantExpData)) {
        if (0 <= PlantExpData[row, "FPKM"] & PlantExpData[row, "FPKM"] <= 5) {
          expressionLevel <- append(expressionLevel, "No Expression")
        }
        else if (5 < PlantExpData[row, "FPKM"] & PlantExpData[row, "FPKM"] <= 25) {
          expressionLevel <- append(expressionLevel, "V.Low Expression")
        }
        else if (25 < PlantExpData[row, "FPKM"] & PlantExpData[row, "FPKM"] <= 55) {
          expressionLevel <- append(expressionLevel, "Low Expression")
        }
        else if (55 < PlantExpData[row, "FPKM"] & PlantExpData[row, "FPKM"] <= 100) {
          expressionLevel <- append(expressionLevel, "Intermediate Expression")
        }
        else if (100 < PlantExpData[row, "FPKM"] & PlantExpData[row, "FPKM"]<= 200) {
          expressionLevel <- append(expressionLevel, "High Expression")
        }
        else if (PlantExpData[row, "FPKM"] > 200) {
          expressionLevel <- append(expressionLevel, "V.High Expression")
        }
      }
      
      PlantExpData <- cbind(PlantExpData, data.frame(ExpressionLevel = expressionLevel))
      
      # Sort genes to hashes based on expression level.
      for (level in unique(PlantExpData$ExpressionLevel)) {
        sampleGenes[[paste(test, "_", t, sep = "")]][[level]] <- PlantExpData[PlantExpData$ExpressionLevel==level,]
      } 
    }
  }
  return(sampleGenes)
}
