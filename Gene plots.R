library(ggplot2)
library(gggenes)
library(readr)
library(readxl)

genomicData <- as.data.frame(read_csv("Data\\Protein coding genes.csv"))
genomicData <- genomicData[,-1]

clusterGenes <- as.data.frame(read_xlsx("Data\\Arabidopsis NLRs.xlsx"))
clusterGenes <- clusterGenes[grepl("cluster", clusterGenes$Clustering),-c(5:6, 8:9)]

clusterGenes$Clustering <- unlist(strsplit(clusterGenes$Clustering, " cluster"))

allClusterGenes <- clusterGenes

for (cluster in unique(clusterGenes$Clustering)) {
  df <- clusterGenes[clusterGenes$Clustering==cluster,]
  
  withinCluster <- genomicData[which(genomicData$seqnames==df$Chromosome[1] & 
                                       genomicData$start > min(df$start) & genomicData$end < max(df$end)),]
  withinCluster <- withinCluster[-which(withinCluster$Gene %in% df$Gene),]
  
  allClusterGenes <- rbind(allClusterGenes, data.frame(Gene = withinCluster$Gene,
                                                 Chromosome = withinCluster$seqnames,
                                                 start = withinCluster$start,
                                                 end = withinCluster$end,
                                                 Clustering = rep(cluster, times = length(withinCluster$Gene))))
}

allClusterGenes <- allClusterGenes[order(allClusterGenes$Gene),]

chrom <- allClusterGenes[allClusterGenes$Chromosome==1,]
labs1 <- ifelse(chrom$Gene %in% 
                  c(chrom[grepl("AT1G27170", chrom$Clustering), "Gene"], chrom[grepl("AT1G63350", chrom$Clustering), "Gene"],
           chrom[grepl("RLM1", chrom$Clustering), "Gene"][c(1,3,5,8,11,14,17,18,21,25,29,31,33)], 
           chrom[grepl("RPP39", chrom$Clustering), "Gene"][c(1:8,10,11,13:15)],
           chrom[grepl("RPP7", chrom$Clustering), "Gene"][c(1,3,5,8,10:12,15,17,18,21,23)], 
           chrom[grepl("RPS5", chrom$Clustering), "Gene"][c(1,3:5,7:11)],
           chrom[grepl("SOC3", chrom$Clustering), "Gene"], 
           chrom[grepl("TN", chrom$Clustering), "Gene"][c(1,3:10,12,13)],
           chrom[grepl("WRR4", chrom$Clustering), "Gene"]), chrom$Gene, "")

labs2 <- ifelse(chrom$Gene %in%
                  c(chrom[grepl("RLM1", chrom$Clustering), "Gene"][c(2,4,9,12,15,19,23,26,30,32)], 
                    chrom[grepl("RPP7", chrom$Clustering), "Gene"][c(2,4,6,9,13,16,19,22,24)],
                    chrom[grepl("RPS5", chrom$Clustering), "Gene"][c(2,6)],
                    chrom[grepl("RPP39", chrom$Clustering), "Gene"][c(9,12)],
                    chrom[grepl("TN", chrom$Clustering), "Gene"][c(2,11)]), chrom$Gene, "")

labs3 <- ifelse(chrom$Gene %in%
                  c(chrom[grepl("RLM1", chrom$Clustering), "Gene"][c(7,10,13,16,22,28)], 
                    chrom[grepl("RPP7", chrom$Clustering), "Gene"][c(7,14,20)]), chrom$Gene, "")

labs4 <- ifelse(chrom$Gene %in%
                  c(chrom[grepl("RLM1", chrom$Clustering), "Gene"][c(6,20,24,27)]), chrom$Gene, "")

plot <- ggplot(chrom, aes(xmin = start, xmax = end, y = Clustering, label = Gene, fill = ifelse(Gene %in% clusterGenes$Gene, Gene, NA))) +
  geom_gene_arrow() + geom_text(mapping = aes(x = end - ((end-start)/2), y = 1.3, label = labs1), size = 4) +
  geom_text(mapping = aes(x = end - ((end-start)/2), label = labs2), size = 4, vjust=2.2) +
  geom_text(mapping = aes(x = end - ((end-start)/2), label = labs3), size = 4, vjust=-3) +
  geom_text(mapping = aes(x = end - ((end-start)/2), label = labs4), size = 4, vjust=3.8) +
  facet_wrap(~ factor(Clustering, unique(chrom$Clustering)), ncol = 1, scales = "free") + labs(x = "Position", y = "Cluster") +
  theme_genes() + theme(legend.position="none", axis.text = element_text(color = "black", size = 12),
                        axis.title = element_text(color = "black", size = 16), plot.margin = unit(c(.5,1,.5,.5), "cm")) +
  scale_fill_discrete(na.value = "transparent")

ggsave("Graphs\\Gene plots\\chrom1 clusters.png", plot = plot, width = 20, height = 12)


chrom <- allClusterGenes[allClusterGenes$Chromosome==2,]

plot <- ggplot(chrom, aes(xmin = start, xmax = end, y = Clustering, label = Gene, fill = Gene)) +
  geom_gene_arrow() + geom_text(mapping = aes(x = end - ((end-start)/2), y = 1.2, label = Gene), size = 4) +
  labs(x = "Position", y = "Cluster") + scale_fill_manual(values = c("coral1","white", "coral")) +
  theme_genes() + theme(legend.position="none", axis.text = element_text(color = "black", size = 12),
                        axis.title = element_text(color = "black", size = 16), plot.margin = unit(c(.5,1,.5,.5), "cm")) 

ggsave("Graphs\\Gene plots\\chrom2 clusters.png", plot = plot, width = 20, height = 2)


chrom <- allClusterGenes[allClusterGenes$Chromosome==3,]

plot <- ggplot(chrom, aes(xmin = start, xmax = end, y = Clustering, label = Gene, fill = ifelse(Gene %in% clusterGenes$Gene, Gene, NA))) +
  geom_gene_arrow() + geom_text(mapping = aes(x = end - ((end-start)/2), y = 1.3, label = Gene), size = 4) +
  facet_wrap(~ factor(Clustering, unique(chrom$Clustering)), ncol = 1, scales = "free") + labs(x = "Position", y = "Cluster") +
  theme_genes() + theme(legend.position="none", axis.text = element_text(color = "black", size = 12),
                        axis.title = element_text(color = "black", size = 16), plot.margin = unit(c(.5,1,.5,.5), "cm")) +
  scale_fill_discrete(na.value = "transparent")

ggsave("Graphs\\Gene plots\\chrom3 clusters.png", plot = plot, width = 20, height = 6.6)


chrom <- allClusterGenes[allClusterGenes$Chromosome==4,]

labs1 <- ifelse(chrom$Gene %in% 
                  c(chrom[grepl("AT4G09360", chrom$Clustering), "Gene"], 
                    chrom[grepl("DSC1", chrom$Clustering), "Gene"],
                    chrom[grepl("RPP4", chrom$Clustering), "Gene"][c(1,3,5,7,9,10,12,13)], 
                    chrom[grepl("RPP2", chrom$Clustering), "Gene"],
                    chrom[grepl("AT4G27190", chrom$Clustering), "Gene"],
                    chrom[grepl("AT4G36140", chrom$Clustering), "Gene"]), chrom$Gene, "")

labs2 <- ifelse(chrom$Gene %in%
                  c(chrom[grepl("RPP4", chrom$Clustering), "Gene"][c(2,4,6,8,11)]), chrom$Gene, "")

plot <- ggplot(chrom, aes(xmin = start, xmax = end, y = Clustering, label = Gene, fill = ifelse(Gene %in% clusterGenes$Gene, Gene, NA))) +
  geom_gene_arrow() + geom_text(mapping = aes(x = end - ((end-start)/2), y = 1.3, label = labs1), size = 4) +
  geom_text(mapping = aes(x = end - ((end-start)/2), label = labs2), size = 4, vjust=2.2) +
  facet_wrap(~ factor(Clustering, unique(chrom$Clustering)), ncol = 1, scales = "free") + labs(x = "Position", y = "Cluster") +
  theme_genes() + theme(legend.position="none", axis.text = element_text(color = "black", size = 12),
                        axis.title = element_text(color = "black", size = 16), plot.margin = unit(c(.5,1,.5,.5), "cm")) +
  scale_fill_discrete(na.value = "transparent")


ggsave("Graphs\\Gene plots\\chrom4 clusters.png", plot = plot, width = 20, height = 7.8)


chrom <- allClusterGenes[allClusterGenes$Chromosome==5,]

labs1 <- ifelse(chrom$Gene %in% 
                  c(chrom[grepl("CHS3", chrom$Clustering), "Gene"][c(1:6,8,9)], 
                    chrom[grepl("AT5G18350", chrom$Clustering), "Gene"],
                    chrom[grepl("AT5G40090", chrom$Clustering), "Gene"], 
                    chrom[grepl("RSG2", chrom$Clustering), "Gene"][c(1,3)],
                    chrom[grepl("TTR1", chrom$Clustering), "Gene"][c(1,2,4,7,9,12,14,18,20)],
                    chrom[grepl("RPS4", chrom$Clustering), "Gene"],
                    chrom[grepl("AT5G45440", chrom$Clustering), "Gene"],
                    chrom[grepl("AT5G47250", chrom$Clustering), "Gene"],
                    chrom[grepl("AT5G48770", chrom$Clustering), "Gene"][c(1,3)],
                    chrom[grepl("NRG1", chrom$Clustering), "Gene"],
                    chrom[grepl("SSI4", chrom$Clustering), "Gene"][c(seq(from = 2, to = 86, by = 6))],
                    chrom[grepl("RPS6", chrom$Clustering), "Gene"][c(1,3,5,9,11,13,15,17,21,25,27)]), chrom$Gene, "")

labs2 <- ifelse(chrom$Gene %in%
                  c(chrom[grepl("CHS3", chrom$Clustering), "Gene"][c(7)],
                    chrom[grepl("RSG2", chrom$Clustering), "Gene"][c(2)],
                    chrom[grepl("TTR1", chrom$Clustering), "Gene"][c(3,6,8,10,13,15,17,19)],
                    chrom[grepl("AT5G48770", chrom$Clustering), "Gene"][c(2)],
                    chrom[grepl("SSI4", chrom$Clustering), "Gene"][c(seq(from = 3, to = 87, by = 6))],
                    chrom[grepl("RPS6", chrom$Clustering), "Gene"][c(2,4,6,10,12,14,16,19,22,24,26,28)]), chrom$Gene, "")

labs3 <- ifelse(chrom$Gene %in%
                c(chrom[grepl("SSI4", chrom$Clustering), "Gene"][c(seq(from = 1, to = 85, by = 6))],
                chrom[grepl("RPS6", chrom$Clustering), "Gene"][c(7,20,23)],
                chrom[grepl("TTR1", chrom$Clustering), "Gene"][c(5,11,16)]), chrom$Gene, "")
                  
labs4 <- ifelse(chrom$Gene %in%
                c(chrom[grepl("SSI4", chrom$Clustering), "Gene"][c(seq(from = 4, to = 88, by = 6))],
                chrom[grepl("RPS6", chrom$Clustering), "Gene"][c(8,18)]), chrom$Gene, "")

labs5 <- ifelse(chrom$Gene %in%
                  c(chrom[grepl("SSI4", chrom$Clustering), "Gene"][c(seq(from = 5, to = 89, by = 6))]), chrom$Gene, "")

labs6 <- ifelse(chrom$Gene %in%
                  c(chrom[grepl("SSI4", chrom$Clustering), "Gene"][c(seq(from = 6, to = 84, by = 6))]), chrom$Gene, "")

plot <- ggplot(chrom, aes(xmin = start, xmax = end, y = Clustering, label = Gene, fill = ifelse(Gene %in% clusterGenes$Gene, Gene, NA))) +
  geom_gene_arrow() + geom_text(mapping = aes(x = end - ((end-start)/2), y = 1.15, label = labs1), size = 4) +
  geom_text(mapping = aes(x = end - ((end-start)/2), label = labs2), size = 4, vjust=2.2) +
  geom_text(mapping = aes(x = end - ((end-start)/2), label = labs3), size = 4, vjust=-2.6) +
  geom_text(mapping = aes(x = end - ((end-start)/2), label = labs4), size = 4, vjust=3.6) +
  geom_text(mapping = aes(x = end - ((end-start)/2), label = labs5), size = 4, vjust=-4) +
  geom_text(mapping = aes(x = end - ((end-start)/2), label = labs6), size = 4, vjust=5.3) +
  facet_wrap(~ factor(Clustering, unique(chrom$Clustering)), ncol = 1, scales = "free") + labs(x = "Position", y = "Cluster") +
  theme_genes() + theme(legend.position="none", axis.text = element_text(color = "black", size = 12),
                        axis.title = element_text(color = "black", size = 16), plot.margin = unit(c(.5,1,.5,.5), "cm")) +
  scale_fill_discrete(na.value = "transparent")

ggsave("Graphs\\Gene plots\\chrom5 clusters.png", plot = plot, width = 24, height = 20)