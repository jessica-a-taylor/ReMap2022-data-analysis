# Get RNA-seq expression data for plants exposed to Hpa.
Hpa_data <- as.data.frame(read_xlsx("Data\\RNA-seq data\\Hpa experiments.xlsx"))
  
NLR_data <- as.data.frame(read_csv("Data\\RNA-seq data\\NLRs_genes.csv", skip = 1))
Hpa_data <- NLR_data[which(NLR_data$Sample %in% Hpa_data$Sample),]
Hpa_data <- Hpa_data[which(Hpa_data$Genotype=="Wild type"),]

Hpa_data <- Hpa_data[,-c(1:2, ncol(Hpa_data):(ncol(Hpa_data)-3))]
 

expressionData_Hpa <- data.frame()

for (t in unique(Hpa_data$Treatment)) {
  df <- Hpa_data[Hpa_data$Treatment == t,]
  
  for (col in 1:ncol(df[,c(which(!is.na(str_match(colnames(df), "^AT[1-9]+G[0-9]+$"))))])) {
    expressionData_Hpa <- rbind(expressionData_Hpa, data.frame(Gene = colnames(df)[col],
                                                               FPKM = mean(df[,col]),
                                                               Treatment = t))
  }
}
  
expressionLevel <- c()

for (row in 1:nrow(expressionData_Hpa)) {
  if (0 <= expressionData_Hpa[row, "FPKM"] & expressionData_Hpa[row, "FPKM"] <= 1) {
    expressionLevel <- append(expressionLevel, "No Expression")
  }
  else if (1 < expressionData_Hpa[row, "FPKM"] & expressionData_Hpa[row, "FPKM"] <= 50) {
    expressionLevel <- append(expressionLevel, "Low Expression")
  }
  else if (50 < expressionData_Hpa[row, "FPKM"] & expressionData_Hpa[row, "FPKM"] <= 100) {
    expressionLevel <- append(expressionLevel, "Intermediate Expression")
  }
  #else if (100 < PlantExpData[row, "FPKM"] & PlantExpData[row, "FPKM"]<= 200) {
  #  expressionLevel <- append(expressionLevel, "High Expression")
  #}
  else if (expressionData_Hpa[row, "FPKM"] > 100) {
    expressionLevel <- append(expressionLevel, "High Expression")
  }
}

expressionData_Hpa <- cbind(expressionData_Hpa, data.frame(ExpressionLevel = expressionLevel))
  

  
# Get PlantExp expression data for plants exposed to Hpa.
# Get the expression data for each experimental condition.
treatments <- list.files(path = "Data\\PlantExp data\\Hpa infected")

expressionData_Hpa <- data.frame()
for (t in treatments) {
  filenames <- list.files(path = paste("Data\\PlantExp data\\Hpa infected\\", t, sep = ""),pattern="*.tsv")
  expressionData <- data.frame()
  
  # Merge the data from all files for each experimental condition.
  for (file in filenames) {
    expressionData <- rbind(expressionData, as.data.frame(read_tsv(paste("Data\\PlantExp data\\Hpa infected\\",t, "\\" , file, sep = ""), show_col_types = FALSE)))
  }
  expressionData <- expressionData[c(which(expressionData$geneId %in% sampleGenes[["NLRs"]]$Gene)),]
  expressionData <- expressionData[,-c(2,3,5)]
  
  for (gene in unique(expressionData$geneId)) {
    df <- expressionData[expressionData$geneId==gene,]
    expressionData_Hpa <- rbind(expressionData_Hpa, data.frame(Gene = gene,
                                                             FPKM = mean(df$FPKM),
                                                             Condition = t))
  }
}


# Add a column to 'expressionData_Hpa' with the expression level of each gene.
expressionLevel <- c()

for (row in 1:nrow(expressionData_Hpa)) {
  if (0 <= expressionData_Hpa[row, "FPKM"] & expressionData_Hpa[row, "FPKM"] <= 1) {
    expressionLevel <- append(expressionLevel, "No Expression")
  }
  else if (1 < expressionData_Hpa[row, "FPKM"] & expressionData_Hpa[row, "FPKM"] <= 50) {
    expressionLevel <- append(expressionLevel, "Low Expression")
  }
  else if (50 < expressionData_Hpa[row, "FPKM"] & expressionData_Hpa[row, "FPKM"] <= 100) {
    expressionLevel <- append(expressionLevel, "Intermediate Expression")
  }
  #else if (100 < PlantExpData[row, "FPKM"] & PlantExpData[row, "FPKM"]<= 200) {
  #  expressionLevel <- append(expressionLevel, "High Expression")
  #}
  else if (expressionData_Hpa[row, "FPKM"] > 100) {
    expressionLevel <- append(expressionLevel, "High Expression")
  }
}

expressionData_Hpa <- cbind(expressionData_Hpa, data.frame(ExpressionLevel = expressionLevel))

