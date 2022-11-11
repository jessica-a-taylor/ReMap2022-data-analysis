# Get 10 sets of random genes and store in a hash.

geneSets <- function(dataToUse) {
  testData <- hash()
  
  for (n in c(1:10)) {
    df <- dataToUse[c(sample(nrow(dataToUse), 200)),]
    df <- df[order(df$Gene),]
    testData[[paste("control", n, sep = "")]] <- df
  }
  return(testData)
}

