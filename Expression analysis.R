expressionList <- c()
for (gene in sampleGenes[["NLRs"]]$Gene) {
  df <- bigExpressionData[["leaf_NLRs"]][bigExpressionData[["leaf_NLRs"]]$Gene==gene,]
  
  expressionList <- append(expressionList, df$Expression)
}

Leaf_NLRs <- sampleGenes[["NLRs"]]
Leaf_NLRs <- cbind(Leaf_NLRs, data.frame(Expression = expressionList))


LeafExpression <- hash(NoExpression = data.frame(),
                    LowExpression = data.frame(),
                    MedExpression = data.frame(),
                    HighExpression = data.frame(),
                    V.HighExpression = data.frame())


for (row in 1:nrow(Leaf_NLRs)) {
  if (Leaf_NLRs[row, "Expression"] == 0) {
    LeafExpression[["NoExpression"]] <- rbind(LeafExpression[["NoExpression"]], Leaf_NLRs[row,])
  }
  else if (0 < Leaf_NLRs[row, "Expression"] & Leaf_NLRs[row, "Expression"] <= 10) {
    LeafExpression[["LowExpression"]] <- rbind(LeafExpression[["LowExpression"]], Leaf_NLRs[row,])
  }
  else if (10 < Leaf_NLRs[row, "Expression"] & Leaf_NLRs[row, "Expression"] <= 50) {
    LeafExpression[["MedExpression"]] <- rbind(LeafExpression[["MedExpression"]], Leaf_NLRs[row,])
  }
  else if (50 < Leaf_NLRs[row, "Expression"] & Leaf_NLRs[row, "Expression"]<= 100) {
    LeafExpression[["HighExpression"]] <- rbind(LeafExpression[["HighExpression"]], Leaf_NLRs[row,])
  }
  else if (Leaf_NLRs[row, "Expression"] > 100) {
    LeafExpression[["V.HighExpression"]] <- rbind(LeafExpression[["V.HighExpression"]], Leaf_NLRs[row,])
  }
}


rm(Leaf_NLRs, expressionList, df)