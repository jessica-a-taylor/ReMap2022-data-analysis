expressionList <- c()
for (gene in sampleGenes[["NLRs"]]$Gene) {
  df <- bigExpressionData[["leaf_NLRs"]][bigExpressionData[["leaf_NLRs"]]$Gene==gene,]
  
  expressionList <- append(expressionList, df$Expression)
}

sampleGenes[["NLRs"]] <- cbind(sampleGenes[["NLRs"]], data.frame(Expression = expressionList))


expressionLevels <- hash(NoExpression = data.frame(),
                    LowExpression = data.frame(),
                    MedExpression = data.frame(),
                    HighExpression = data.frame(),
                    V.HighExpression = data.frame())

for (row in 1:nrow(sampleGenes[["NLRs"]])) {
  if (sampleGenes[["NLRs"]][row, "Expression"] == 0) {
    expressionLevels[["NoExpression"]] <- rbind(expressionLevels[["NoExpression"]], sampleGenes[["NLRs"]][row,])
  }
  else if (0 < bigExpressionData[["leaf_NLRs"]][row, "Expression"] & bigExpressionData[["leaf_NLRs"]][row, "Expression"] <= 10) {
    expressionLevels[["LowExpression"]] <- rbind(expressionLevels[["LowExpression"]], sampleGenes[["NLRs"]][row,])
  }
  else if (10 < bigExpressionData[["leaf_NLRs"]][row, "Expression"] & bigExpressionData[["leaf_NLRs"]][row, "Expression"] <= 50) {
    expressionLevels[["MedExpression"]] <- rbind(expressionLevels[["MedExpression"]], sampleGenes[["NLRs"]][row,])
  }
  else if (50 < bigExpressionData[["leaf_NLRs"]][row, "Expression"] & bigExpressionData[["leaf_NLRs"]][row, "Expression"]<= 100) {
    expressionLevels[["HighExpression"]] <- rbind(expressionLevels[["HighExpression"]], sampleGenes[["NLRs"]][row,])
  }
  else if (bigExpressionData[["leaf_NLRs"]][row, "Expression"] > 100) {
    expressionLevels[["V.HighExpression"]] <- rbind(expressionLevels[["V.HighExpression"]], sampleGenes[["NLRs"]][row,])
  }
}


