RNA_seqAnalysis <- function(dataToUse, exLevel) {
  
  treatments <- c("Mock", "untreated", "control", "Not treated", "mock", "unstressed control", "none",
                  "normal condition","healthy", "NA")
  
  tissue <- c("leaves", "root", "seedlings")
  
  expressionData <- hash()
  
  for (test in names(dataToUse)) {
    print(test)
    
    expressionData[[test]] <- as.data.frame(read_csv(paste("Data\\RNA-seq data\\", test, "_genes.csv", sep = ""), skip = 1))
    
    expressionData[[test]] <- expressionData[[test]][,-c(1:2, ncol(expressionData[[test]]):(ncol(expressionData[[test]])-3))]
    
    expressionData[[test]] <- expressionData[[test]][expressionData[[test]]$Ecotype=="Col-0",]
    expressionData[[test]] <- expressionData[[test]][expressionData[[test]]$Genotype=="wild type",]
    
    expressionData[[test]] <- expressionData[[test]][c(which(expressionData[[test]]$Treatment %in% treatments)),]
    expressionData[[test]] <- expressionData[[test]][c(which(expressionData[[test]]$Tissue %in% tissue)),]
    
    
    for (t in tissue) {
      df <- expressionData[[test]][expressionData[[test]]$Tissue == t,]
      
      expressionData_tissue <- data.frame()
      
      for (col in 1:ncol(df[,c(which(!is.na(str_match(colnames(df), "^AT[1-9]+G[0-9]+$"))))])) {
      expressionData_tissue <- rbind(expressionData_tissue, data.frame(Gene = colnames(df)[col],
                                                                         Expression = mean(df[,col]),
                                                                         Tissue = t))
      }
      
      expressionData_tissue <- expressionData_tissue[order(expressionData_tissue$Gene),]
      
      expressionLevel <- c()
      
      for (gene in expressionData_tissue$Gene) {
        if (0 <= expressionData_tissue[expressionData_tissue$Gene==gene, "Expression"] & expressionData_tissue[expressionData_tissue$Gene==gene, "Expression"] <= 1) {
          expressionLevel <- append(expressionLevel, "No Expression")
        }
        else if (1 < expressionData_tissue[expressionData_tissue$Gene==gene, "Expression"] & expressionData_tissue[expressionData_tissue$Gene==gene, "Expression"] <= 50) {
          expressionLevel <- append(expressionLevel, "Low Expression")
        }
        else if (50 < expressionData_tissue[expressionData_tissue$Gene==gene, "Expression"] & expressionData_tissue[expressionData_tissue$Gene==gene, "Expression"] <= 100) {
          expressionLevel <- append(expressionLevel, "Intermediate Expression")
        }
        #else if (50 < expressionData_tissue[expressionData_tissue$Gene==gene, "Expression"] & expressionData_tissue[expressionData_tissue$Gene==gene, "Expression"]<= 100) {
         # expressionLevel <- append(expressionLevel, "High Expression")
        #}
        else if (expressionData_tissue[expressionData_tissue$Gene==gene, "Expression"] > 100) {
          expressionLevel <- append(expressionLevel, "High Expression")
        }
      }
      
      expressionData_tissue <- cbind(expressionData_tissue, data.frame(ExpressionLevel = expressionLevel))

      # Sort genes to hashes based on tissue and expression level.
     for (level in exLevel) {
        df1 <- expressionData_tissue[expressionData_tissue$ExpressionLevel==level,]
        df1 <- df1[c(which(df1$Gene %in% dataToUse[[test]]$Gene)),]
        
        dataToUse[[paste(test, "_", t, sep = "")]][[level]] <- dataToUse[[test]][c(which(dataToUse[[test]]$Gene %in% df1$Gene)),]
        
        dataToUse[[paste(test, "_", t, sep = "")]][[level]] <- cbind(dataToUse[[paste(test, "_", t, sep = "")]][[level]], df1[,c(2:4)])
      }   
      write.csv(expressionData_tissue, file = paste("Data\\RNA-seq data\\", t,"\\Expression data.csv", sep = ""))
    }
  }
  return(dataToUse)
}
