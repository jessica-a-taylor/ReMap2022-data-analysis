# Filter and merge all PlantExp data into one dataframe.
tissue <- list.files(path = "Data\\PlantExp data")

PlantExp <- data.frame()

for (t in tissue) {
  filenames <- list.files(path = paste("Data\\PlantExp data\\", t, sep = ""),pattern="*.tsv")
  
  for (file in filenames) {
    data <- as.data.frame(read_tsv(paste("Data\\PlantExp data\\", t, "\\", file, sep = ""), show_col_types = FALSE))
    data <- data[c(which(data$geneId %in% allGenes)),]
    data <- data[,-c(2,3,5)]
    
    dataDF <- data.frame()
    
    for (gene in unique(data$geneId)) {
      df <- data[data$geneId==gene,]
      
        
      dataDF <- rbind(dataDF, data.frame(geneId = df$geneId[1],
                                         FPKM = mean(df$FPKM),
                                         tissue = t))
    }
  }
  PlantExp <- rbind(PlantExp, dataDF)
}

# For each tissue type, calculate the mean expression of each gene.
expressionHash <- hash()

for (t in tissue) {
  data <- PlantExp[PlantExp$tissue == t,]
  
  dataDF <- data.frame()
  
  for (gene in unique(data$geneId)) {
    df <- data[data$geneId==gene,]
  
    if (nrow(df) > 1) {
      df <- data.frame(geneId = df$geneId[1],
                       FPKM = mean(df$FPKM),
                       tissue = t)
    }
    dataDF <- rbind(dataDF, df)
  }
  expressionHash[[t]] <- dataDF
}

