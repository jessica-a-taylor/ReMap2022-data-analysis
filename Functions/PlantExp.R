# Function to get filtered expression data for each set of sample genes in each tissue. 
PlantExp <- function(dataToUse, exLevel) {

  
  sampleGeneSets <- names(dataToUse)
 
  # For each tissue type...
  for (t in c("leaves","root","seedlings")) {
    
    # Create a list of files in the folder for a particular tissue type.
    filenames <- list.files(path = paste("Data\\PlantExp data\\", t, sep = ""),pattern="*.tsv")
    expressionData <- data.frame()
    
    # Merge the data from all files.
    for (file in filenames) {
      expressionData <- rbind(expressionData, as.data.frame(read_tsv(paste("Data\\PlantExp data\\", t, "\\", file, sep = ""), show_col_types = FALSE)))
    }
    
    # For each set of sample genes...
    for (test in sampleGeneSets) {
      PlantExpData <- data.frame()

      # Create a hash.
      dataToUse[[paste(test, "_", t, sep = "")]] <- hash()
      
      # Filter for the genes in the sample gene set.
      sampleData <- expressionData[c(which(expressionData$geneId %in% dataToUse[[test]]$Gene)),]
      sampleData <- sampleData[,-c(2,3,5)]
      sampleData <- sampleData[order(sampleData$geneId),]
      
      # For each gene, calculate the mean expression across experiments.
      for (gene in unique(sampleData$geneId)) {
        df <- sampleData[sampleData$geneId==gene,]
        
        
        PlantExpData <- rbind(PlantExpData, data.frame(geneId = gene,
                                                       FPKM = mean(df$FPKM),
                                                       tissue = t))
      }
      # Add a column to PlantExpData with the expression level of each gene.
      expressionLevel <- c()
      
      for (row in 1:nrow(PlantExpData)) {
        if (0 <= PlantExpData[row, "FPKM"] & PlantExpData[row, "FPKM"] <= 1) {
          expressionLevel <- append(expressionLevel, "No Expression")
        }
        else if (1 < PlantExpData[row, "FPKM"] & PlantExpData[row, "FPKM"] <= 50) {
          expressionLevel <- append(expressionLevel, "Low Expression")
        }
        else if (50 < PlantExpData[row, "FPKM"] & PlantExpData[row, "FPKM"] <= 100) {
          expressionLevel <- append(expressionLevel, "Intermediate Expression")
        }
        #else if (100 < PlantExpData[row, "FPKM"] & PlantExpData[row, "FPKM"]<= 200) {
        #  expressionLevel <- append(expressionLevel, "High Expression")
        #}
        else if (PlantExpData[row, "FPKM"] > 100) {
          expressionLevel <- append(expressionLevel, "High Expression")
        }
      }
      
      PlantExpData <- cbind(PlantExpData, data.frame(ExpressionLevel = expressionLevel))
      
      # Sort genes to hashes based on expression level.
      for (level in exLevel) {
        df <- PlantExpData[PlantExpData$ExpressionLevel==level,]

        dataToUse[[paste(test, "_", t, sep = "")]][[level]] <- dataToUse[[test]][c(which(dataToUse[[test]]$Gene %in% df$geneId)),]
        
        dataToUse[[paste(test, "_", t, sep = "")]][[level]] <- cbind(dataToUse[[paste(test, "_", t, sep = "")]][[level]], df[,c(2:4)])
      } 
    }
    write.csv(PlantExpData, file = paste("Data\\PlantExp data\\", t,"\\Expression data.csv", sep = ""))
  }
  return(dataToUse)
}
