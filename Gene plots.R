library(ggplot2)
library(gggenes)

clusterGenes <- as.data.frame(read_xlsx("Data\\Arabidopsis NLRs.xlsx"))
clusterGenes <- clusterGenes[grepl("cluster", clusterGenes$Clustering),-c(5:10, 12)]

chrom1 <- clusterGenes[clusterGenes$Chromosome==1,]
labs1 <- ifelse(chrom1$Gene %in% 
                  c(chrom1[grepl("AT1G27170", chrom1$Clustering), "Gene"], chrom1[grepl("AT1G63350", chrom1$Clustering), "Gene"],
           chrom1[grepl("RLM1", chrom1$Clustering), "Gene"][c(1,3,4,6,7)], chrom1[grepl("RPP39", chrom1$Clustering), "Gene"],
           chrom1[grepl("RPP7", chrom1$Clustering), "Gene"][c(1,3,4,8)], chrom1[grepl("RPS5", chrom1$Clustering), "Gene"],
           chrom1[grepl("SOC3", chrom1$Clustering), "Gene"], chrom1[grepl("TN", chrom1$Clustering), "Gene"][c(1,2,4:9,11,12)],
           chrom1[grepl("WRR4", chrom1$Clustering), "Gene"]), chrom1$Gene, "")

labs2 <- ifelse(!chrom1$Gene %in%
                  c(chrom1[grepl("RLM1", chrom1$Clustering), "Gene"][c(2,5)], chrom1[grepl("RPP7", chrom1$Clustering), "Gene"][c(2)],
           chrom1[grepl("TN", chrom1$Clustering), "Gene"][c(3,10)]), "",chrom1$Gene)

plot <- ggplot(chrom1, aes(xmin = start, xmax = end, y = Clustering, label = Gene, fill = Gene)) +
  geom_gene_arrow() + geom_text(mapping = aes(x = end - ((end-start)/2), y = 1.3, label = labs1), size = 3) +
  geom_text(mapping = aes(x = end - ((end-start)/2), label = labs2), size = 3, vjust=2.2) +
  facet_wrap(~ factor(Clustering, unique(chrom1$Clustering)), ncol = 1, scales = "free") + labs(x = "Position", y = "Cluster") +
  theme_genes() + theme(legend.position="none", axis.text = element_text(color = "black", size = 10),
                        axis.title = element_text(color = "black", size = 14))