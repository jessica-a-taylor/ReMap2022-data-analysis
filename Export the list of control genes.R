# Export the list of control genes.
for (test in names(sampleGenes)[grepl("control", names(sampleGenes))]) {
  write(paste(sampleGenes[[test]]$Gene, collapse= ",", sep = ""), file = paste("Data\\RNA-seq data\\",test, "_genes.txt", sep = ""))
}
