

# Function to get filtered expression data for each set of sample genes in each tissue. 
expressionFiltered <- function(dataToUse) {
  
  dataToUse[[paste(test, "_", t, sep = "")]] <- hash()
  
  
  otherExpressionData <- as.data.frame(read_xlsx("Data\\result_NLRs.xlsx"))
  
  treatments <- c("Mock", "untreated", "control", "Not treated", "mock", "unstressed control", "none",
                  "normal condition","healthy", "NA")
  
  tissue <- c("leaves", "rosette leaf", "aerial seedling")
  
  # Remove first two columns.
  otherExpressionData <- otherExpressionData[,-c(1:2, 162:165)]
  
  otherExpressionData <- otherExpressionData[otherExpressionData$Ecotype=="Col-0",]
  otherExpressionData <- otherExpressionData[otherExpressionData$Genotype=="wild type",]
  
  otherExpressionData <- otherExpressionData[c(which(otherExpressionData$Treatment %in% treatments)),]
  otherExpressionData <- otherExpressionData[c(which(otherExpressionData$Tissue %in% tissue)),]
  
  otherExpressionData <- otherExpressionData[,-c(156:159)]
  
  otherNLRexpression <- data.frame()
  
  for (col in 1:ncol(otherExpressionData)) {
    otherNLRexpression <- rbind(otherNLRexpression, data.frame(Gene = colnames(otherExpressionData)[col],
                                                               Expression = mean(otherExpressionData[,col])))
  }
  
  NLRlevel <- c()
  
  for (gene in sampleGenes[["NLRs"]]$Gene) {
    if (0 <= otherNLRexpression[otherNLRexpression$Gene==gene, "Expression"] & otherNLRexpression[otherNLRexpression$Gene==gene, "Expression"] <= 10) {
      NLRlevel <- append(NLRlevel, "No Expression")
    }
    else if (10 < otherNLRexpression[otherNLRexpression$Gene==gene, "Expression"] & otherNLRexpression[otherNLRexpression$Gene==gene, "Expression"] <= 50) {
      NLRlevel <- append(NLRlevel, "Low Expression")
    }
    else if (50 < otherNLRexpression[otherNLRexpression$Gene==gene, "Expression"] & otherNLRexpression[otherNLRexpression$Gene==gene, "Expression"] <= 100) {
      NLRlevel <- append(NLRlevel, "Intermediate Expression")
    }
    else if (100 < otherNLRexpression[otherNLRexpression$Gene==gene, "Expression"] & otherNLRexpression[otherNLRexpression$Gene==gene, "Expression"]<= 200) {
      NLRlevel <- append(NLRlevel, "High Expression")
    }
    else if (otherNLRexpression[otherNLRexpression$Gene==gene, "Expression"] > 200) {
      NLRlevel <- append(NLRlevel, "V.High Expression")
    }
  }
  
  sampleGenes[["NLRs"]] <- cbind(sampleGenes[["NLRs"]], data.frame(Level = NLRlevel))
  
  # Sort genes to hashes based on expression level.
  for (level in exLevel) {
    df <- PlantExpData[PlantExpData$ExpressionLevel==level,]
    
    dataToUse[[paste(test, "_", t, sep = "")]][[level]] <- dataToUse[[test]][c(which(dataToUse[[test]]$Gene %in% df$geneId)),]
    
    dataToUse[[paste(test, "_", t, sep = "")]][[level]] <- cbind(dataToUse[[paste(test, "_", t, sep = "")]][[level]], df[,c(2:4)])
  } 
}
