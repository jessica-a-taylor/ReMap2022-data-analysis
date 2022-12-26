Ding_ExpressionData <- as.data.frame(read_xlsx("Data\\ACRs Ding et al., 2021.xlsx", sheet = 1))

Ding_ExpressionData <- Ding_ExpressionData[which(Ding_ExpressionData$Gene %in% sampleGenes[["NLRs"]] $Gene),]

controlExpression <- c()
ETIexpression <- c()

for (row in 1:nrow(Ding_ExpressionData)) {
  controlExpression <- append(controlExpression, mean(as.numeric(Ding_ExpressionData[row, 2:4])))
  ETIexpression <- append(ETIexpression, mean(as.numeric(Ding_ExpressionData[row, 5:7])))
}

Ding_ExpressionData <- data.frame(Gene = Ding_ExpressionData$Gene,
                                  Control = controlExpression,
                                  ETI = ETIexpression)

# Perform ReMap analysis to determine the chromatin modification enrichment in each R-gene.
# Determine which genes appear to be the most responsive to infection - are there any trend in their chromatin modifications.


# Which R-genes overlap with ACRs?

# Which TFs are the R-genes associated with - is there a correlation between the enrichment of particlar chromatin 
# modifications and particular TFs?

# Are there similarities in chromatin modification and TF encrichment between co-expressed genes?

# Which R-genes are asssociated with the nuclear envelope?

