# Calculate the difference in the proportion of each R-gene region covered by each modification between R-genes of the same cluster.

modificationDifference_Clusters <- data.frame()

for (mod in modList[[modMiniList]]) {
  df <- allResultsProportions_ClusterAnalysis[allResultsProportions_ClusterAnalysis$Modification == mod,]
  
  for (r in unique(df$Region)) {
    df1 <- df[df$Region == r,]
    
    for (cluster in unique(df1[grepl("cluster", df1$Expression),"Expression"])) {
      df2 <- df1[df1$Expression == cluster,]
      
      # For all pairwise comparisons,
      for (i in 1:nrow(df2)) {
        for (j in 1:nrow(df2)) {
          if (i != j & j > i) {
            # calculate the % difference in the proportion of each gene region covered by each modification,
            # using the abs() function to make all values positive.
           modificationDifference_Clusters <- rbind(modificationDifference_Clusters, 
                                                            data.frame(Region = r,
                                                                       Modification = mod,
                                                                       axisGroup = df2$axisGroup[1],
                                                                       ModificationDifference = abs(df2[i, "Proportion"]-df2[j, "Proportion"]),
                                                                       Cluster = cluster))
          }
          else next
        }
      }
    }
  }
  print(mod)
}

write.csv(modificationDifference_Clusters, file = paste(modMiniList, " - Cluster modification comparisons.csv", sep = ""))

# Calculate the difference in the proportion of each gene region covered by each modification between all R-genes.
modificationDifference_all <- data.frame()


for (mod in modList[[modMiniList]]) {
  df <- allResultsProportions_ClusterAnalysis[allResultsProportions_ClusterAnalysis$Modification == mod,]
  
  for (r in unique(df$Region)) {
    df1 <- df[df$Region == r,]
    
    for (i in 1:nrow(df1)) {
      for (j in 1:nrow(df1)) {
        if (i != j & j > i) {
          modificationDifference_all <- rbind(modificationDifference_all, 
                                              data.frame(Region = r,
                                                         Modification = mod,
                                                         axisGroup = df1$axisGroup[1],
                                                         ModificationDifference = abs(df1[i, "Proportion"]-df1[j, "Proportion"])))
        }
      }
    }
    print(r)
  }
  print(mod)
}

write.csv(modificationDifference_all, file = paste("Data\\",modMiniList, " - All R-gene modification comparisons.csv", sep = ""))