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

existingSets <- function(dataToUse) {
  testData <- hash()
  
  files <- list.files(path = paste("Data\\RNA-seq data\\"), pattern= "*.csv")
  
  for (file in files[which(grepl("control", files))]) {
    
    findGenes <- as.data.frame(read_csv(paste("Data\\RNA-seq data\\", file, sep = ""), skip = 1))
    
    testData[[str_match(file, "^([a-zA-Z]+[0-9]+)_[a-zA-Z]+.[a-zA-Z]+$")[,2]]] <- dataToUse[c(which(dataToUse$Gene %in% colnames(findGenes))),]
  }
  return(testData)
}