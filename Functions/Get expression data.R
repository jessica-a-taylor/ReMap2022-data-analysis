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
    meanExpression <- data.frame(matrix(ncol = 155, nrow = 0))
    
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
    
    rootMeans <- meanExpression[bigExpressionData[[test]]$Tissue =="root",][,-length(bigExpressionData[[test]])]
    
    rootExpression <- data.frame(Gene = character(),
                                 Expression = numeric())
    
    for (col in 1:ncol(rootMeans)) {
      rootExpression <- rbind(rootExpression, data.frame(Gene = colnames(rootMeans)[col],
                                                         Expression = mean(rootMeans[,col])))
    }
    
    bigExpressionData[[paste("leaf_", test, sep = "")]] <- leafExpression
    bigExpressionData[[paste("root_", test, sep = "")]] <- rootExpression
  }

    # Make a list for gene expression in leaves and roots. 
    # Considered not expressed if the expression level in below 0.
    for (test in names(bigExpressionData)[c(11:21)]) {
      LeafActive <- bigExpressionData[[test]][-c(which(bigExpressionData[[test]]$Expression==0)),]
      LeafSilent <- bigExpressionData[[test]][c(which(bigExpressionData[[test]]$Expression==0)),]
      
      sampleGenes[[paste("LeafActive_", test, sep = "")]] <- withoutTEs[c(which(withoutTEs$Gene %in% LeafActive$Gene)),]
      sampleGenes[[paste("LeafSilent_", test, sep = "")]] <- withoutTEs[c(which(withoutTEs$Gene %in% LeafSilent$Gene)),]
    }
    
    for (test in names(bigExpressionData)[c(23:33)]) {
      RootActive <- bigExpressionData[[test]][-c(which(bigExpressionData[[test]]$Expression==0)),]
      RootSilent <- bigExpressionData[[test]][c(which(bigExpressionData[[test]]$Expression==0)),]
      
      sampleGenes[[paste("RootActive_", test, sep = "")]] <- withoutTEs[c(which(withoutTEs$Gene %in% RootActive$Gene)),]
      sampleGenes[[paste("RootSilent_", test, sep = "")]] <- withoutTEs[c(which(withoutTEs$Gene %in% RootSilent$Gene)),]
    }
  return(sampleGenes)
}






