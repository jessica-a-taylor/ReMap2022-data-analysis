otherExpressionAnalysis <- function(sampleGenes) {
  otherExpressionData <- as.data.frame(read_xlsx("Data\\result_NLRs.xlsx"))
  
  sampleGenes[["otherNLRs"]] <- hash()
  
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
  FPKM <- c()
  
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
    FPKM <- append(FPKM, otherNLRexpression[otherNLRexpression$Gene==gene, "Expression"])
  }
  
  sampleGenes[["NLRs"]] <- cbind(sampleGenes[["NLRs"]], data.frame(RNA_seqLevel = NLRlevel))
  sampleGenes[["NLRs"]] <- cbind(sampleGenes[["NLRs"]], data.frame(RNA_seqExpression = FPKM))
  
  # Sort genes to hashes based on expression level.
  for (level in unique(sampleGenes[["NLRs"]]$Level)) {
    sampleGenes[["otherNLRs"]][[level]] <- sampleGenes[["NLRs"]][sampleGenes[["NLRs"]]$Level==level,]
  } 
  return(sampleGenes)
}
