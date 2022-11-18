# Function to get filtered expression data for each set of sample genes in each tissue. 
PlantExp <- function(dataToUse) {
  
  # Get a list of tissue types for which there is expression data.
  tissue <- list.files(path = "Data\\PlantExp data")
  
  for (test in names(dataToUse)) {

    for (t in tissue) {
      filenames <- list.files(path = paste("Data\\PlantExp data\\", t, sep = ""),pattern="*.tsv")
      PlantExpData <- data.frame()
      
      for (file in filenames) {
        data <- as.data.frame(read_tsv(paste("Data\\PlantExp data\\", t, "\\", file, sep = ""), show_col_types = FALSE))
        data <- data[c(which(data$geneId %in% dataToUse[[test]]$Gene)),]
        data <- data[,-c(2,3,5)]
        
        dataDF <- data.frame()
        
        for (gene in unique(data$geneId)) {
          df <- data[data$geneId==gene,]
          
          
          dataDF <- rbind(dataDF, data.frame(geneId = df$geneId[1],
                                             FPKM = mean(df$FPKM),
                                             tissue = t))
        }
      } 
      PlantExpData <- rbind(PlantExpData, dataDF)
      
      # Calculate the mean expression of each gene from each experiment.
      dataDF <- data.frame()
      
      for (gene in unique(PlantExpData$geneId)) {
        df <- PlantExpData[PlantExpData$geneId==gene,]
        
        dataDF <- rbind(dataDF, data.frame(geneId = df$geneId[1],
                                           FPKM = mean(df$FPKM),
                                           tissue = t))
      }
      sampleGenes[[paste(test, "_", t, sep = "")]] <- dataDF
    }
  }
  # Set thresholds for expression levels.
  tissueExpression <- hash(NoExpression = data.frame(),
                           V.LowExpression = data.frame(),
                           LowExpression = data.frame(),
                           MedExpression = data.frame(),
                           HighExpression = data.frame(),
                           V.HighExpression = data.frame())
  
  for (test in names(sampleGenes)[-c(1,5,9,13,17,21,25,29,33,37,41)]) {
    
    expressionLevel <- c()
    for (row in 1:nrow(sampleGenes[[test]])) {
      if (0 <= sampleGenes[[test]][row, "FPKM"] & sampleGenes[[test]][row, "FPKM"] <= 5) {
        expressionLevel <- append(expressionLevel, "No Expression")
      }
      else if (5 < sampleGenes[[test]][row, "FPKM"] & sampleGenes[[test]][row, "FPKM"] <= 25) {
        expressionLevel <- append(expressionLevel, "V.Low Expression")
      }
      else if (25 < sampleGenes[[test]][row, "FPKM"] & sampleGenes[[test]][row, "FPKM"] <= 55) {
        expressionLevel <- append(expressionLevel, "Low Expression")
      }
      else if (55 < sampleGenes[[test]][row, "FPKM"] & sampleGenes[[test]][row, "FPKM"] <= 100) {
        expressionLevel <- append(expressionLevel, "Intermediate Expression")
      }
      else if (100 < sampleGenes[[test]][row, "FPKM"] & sampleGenes[[test]][row, "FPKM"]<= 200) {
        expressionLevel <- append(expressionLevel, "High Expression")
      }
      else if (sampleGenes[[test]][row, "FPKM"] > 200) {
        expressionLevel <- append(expressionLevel, "V.High Expression")
      }
    }
    sampleGenes[[test]] <- cbind(sampleGenes[[test]], data.frame(ExpressionLevel = expressionLevel))
  }
}
