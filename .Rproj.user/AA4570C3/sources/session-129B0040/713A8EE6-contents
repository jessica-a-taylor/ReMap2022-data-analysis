# Function to get filtered expression data for each set of sample genes in each tissue. 

expressionFiltered <- function(bigExpressionData, sampleGenes) {
  
  # Filter each expression data for Col-0, leaf/root tissue, and wild-type conditions.
  treatments <- c("Mock", "untreated", "control", "Not treated", "mock", "unstressed control", "none",
                  "normal condition","healthy", "NA")
  
  tissue <- c("leaves", "root", "rosette leaf", "aerial seedling")
  
  for (test in names(bigExpressionData)) {
    # Remove first two columns.
    bigExpressionData[[test]] <- bigExpressionData[[test]][,-c(1:2)]
    
    bigExpressionData[[test]] <- bigExpressionData[[test]][bigExpressionData[[test]]$Ecotype=="Col-0",]
    bigExpressionData[[test]] <- bigExpressionData[[test]][bigExpressionData[[test]]$Genotype=="wild type",]
    
    bigExpressionData[[test]] <- bigExpressionData[[test]][c(which(bigExpressionData[[test]]$Treatment %in% treatments)),]
    bigExpressionData[[test]] <- bigExpressionData[[test]][c(which(bigExpressionData[[test]]$Tissue %in% tissue)),]
    
    # Get the mean expression in each tissue type.
    meanExpression <- data.frame(matrix(ncol = length(bigExpressionData[[test]])-8, nrow = 0))
    
    for(t in tissue) {
      df <- bigExpressionData[[test]][bigExpressionData[[test]]$Tissue==t,]
      
      tissueMean <- c()
      for (col in colnames(df)[1:(length(bigExpressionData[[test]])-8)]) {
        tissueMean <- append(tissueMean, mean(df[,col]))
      }
      meanExpression <- rbind(meanExpression, tissueMean) 
    }
    colnames(meanExpression) <- c(colnames(bigExpressionData[[test]][1:(length(bigExpressionData[[test]])-8)]))
    
    bigExpressionData[[test]] <- meanExpression
    bigExpressionData[[test]]$Tissue <- tissue
  }
  
  
  # Rewrite sampleGenes to match genes with downloaded expression data.
  for (test in names(sampleGenes)) {
    sampleGenes[[test]] <- withoutTEs[c(which(withoutTEs$Gene %in% colnames(bigExpressionData[[test]]))),]
  }
  
  # Create separate dataframes for the mean expression level in leaves and roots.
  leafTissue <- c("leaves", "rosette leaf", "aerial seedling")
  
  for (test in names(bigExpressionData)) {
    leafMeans <- bigExpressionData[[test]][bigExpressionData[[test]]$Tissue %in% leafTissue,][,-length(bigExpressionData[[test]])]
    
    leafExpression <- data.frame(Gene = character(),
                                 Expression = numeric())
    
    for (col in 1:ncol(leafMeans)) {
      leafExpression <- rbind(leafExpression, data.frame(Gene = colnames(leafMeans)[col],
                                                         Expression = mean(leafMeans[,col])))
    }
    
    rootMeans <- bigExpressionData[[test]][bigExpressionData[[test]]$Tissue =="root",][,-length(bigExpressionData[[test]])]
    
    rootExpression <- data.frame(Gene = character(),
                                 Expression = numeric())
    
    for (col in 1:ncol(rootMeans)) {
      rootExpression <- rbind(rootExpression, data.frame(Gene = colnames(rootMeans)[col],
                                                         Expression = mean(rootMeans[,col])))
    }
    
    bigExpressionData[[paste("leaf_", test, sep = "")]] <- leafExpression
    bigExpressionData[[paste("root_", test, sep = "")]] <- rootExpression
  }
  
  # Create a hash for leaf-expressed and root-expressed genes.
  leafExpression <- hash()
  rootExpression <- hash()
  
  # Sort leaf- and root-expressed genes into groups based on expression level.
  tissues <- c("leaf", "root")
  
  for (t in tissues) {
    for (test in names(sampleGenes)) {
      df <- sampleGenes[[test]]
      
      expressionList <- c()
      for (gene in sampleGenes[[test]]$Gene) {
        geneExpression <- bigExpressionData[[paste(t, "_", test, sep = "")]][bigExpressionData[[paste(t, "_", test, sep = "")]]$Gene==gene,]
        
        expressionList <- append(expressionList, geneExpression$Expression)
      }
      df <- cbind(df, data.frame(Expression = expressionList))
      
      tissueExpression <- hash(NoExpression = data.frame(),
                               LowExpression = data.frame(),
                               MedExpression = data.frame(),
                               HighExpression = data.frame(),
                               V.HighExpression = data.frame())
      
      
      for (row in 1:nrow(df)) {
        if (0 <= df[row, "Expression"] & df[row, "Expression"] <= 5) {
          tissueExpression[["NoExpression"]] <- rbind(tissueExpression[["NoExpression"]], df[row,])
        }
        else if (5 < df[row, "Expression"] & df[row, "Expression"] <= 30) {
          tissueExpression[["LowExpression"]] <- rbind(tissueExpression[["LowExpression"]], df[row,])
        }
        else if (30 < df[row, "Expression"] & df[row, "Expression"] <= 60) {
          tissueExpression[["MedExpression"]] <- rbind(tissueExpression[["MedExpression"]], df[row,])
        }
        else if (60 < df[row, "Expression"] & df[row, "Expression"]<= 90) {
          tissueExpression[["HighExpression"]] <- rbind(tissueExpression[["HighExpression"]], df[row,])
        }
        else if (df[row, "Expression"] > 90) {
          tissueExpression[["V.HighExpression"]] <- rbind(tissueExpression[["V.HighExpression"]], df[row,])
        }
      }
      if (t == "leaf") {
        leafExpression[[test]] <- tissueExpression
      }
      else rootExpression[[test]] <- tissueExpression
    }
  }
  sampleGenes[["leafExpression"]] <- leafExpression
  sampleGenes[["rootExpression"]] <- rootExpression
  
  return(sampleGenes)
}

