source('03_kallisto_to_DESeq2.R')


DESeq2Table <- DESeq(ddsG)

# Pulls gene count info from table
GeneCounts <- counts(DESeq2Table)
# Just keeps genes where there's at least one count
idx.nz <- apply(GeneCounts, 1, function(x) { all(x > 0)})
sum(idx.nz)

# Calculate size (sequence depth) correction
DESeq2Table <- estimateSizeFactors(DESeq2Table)
sizeFactors(DESeq2Table)

# Checks to see whether size factor correction is working.
# Expect to see all plots looking about the same
library(geneplotter)
multidensity( counts(DESeq2Table, normalized = T)[idx.nz ,],
              xlab="mean counts", xlim=c(0, 1000))

#
multiecdf( counts(DESeq2Table, normalized = T)[idx.nz ,],
           xlab="mean counts", xlim=c(0, 1000))

# heatmap. Taken from DESeq2 tutorial
rld <- rlogTransformation(DESeq2Table, blind=TRUE)
colnames(rld) <- data.frame(colnames(rld)) %>% left_join(.,SampleTable,by=c("colnames.rld."="unique.ID")) %>% mutate(ID=paste(Genotype,Diet,lane.ID,sample.ID,sep='_')) %>% .[['ID']]
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
library("RColorBrewer")
library("pheatmap")
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# looks crazy busy, because technical replicates are not collapsed together
# demonstrates that, at least, the tech replicates are related
# so it's safe to collpase them together
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)         
# collapse technical replicates and repeat heatmap
ddsColl <- collapseReplicates(dds, dds$sample.ID, dds$unique.ID)
DESeq2Table <- DESeq(ddsColl)
rld <- rlogTransformation(DESeq2Table, blind=TRUE)
colnames(rld) <- data.frame(colnames(rld)) %>% mutate(sample.ID=colnames.rld.) %>% right_join(.,SampleTable) %>% dplyr::select(sample.ID, Genotype, Diet) %>% distinct()  %>% mutate(ID=paste(Genotype,Diet,sample.ID,sep='_')) %>% .[['ID']]
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
# not seeing any structure that matches biological knowldge (Genotype and/or Diet not grouping together)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)         


# tsne clustering
library(ggrepel)
library(Rtsne)
library(tidyr)
tsne_out <- Rtsne(as.matrix(t(assay(rld))),perplexity = 5)
tsne_plot <- data.frame(tsne_out$Y)
tsne_plot$ID <- colnames(rld)
tsne_plot %>% mutate(ID=gsub("B12_Deficient","B12Deficient",colnames(rld))) %>% separate(ID,c("Genotype","Diet","B12","Sample"),sep='_') %>% 
  ggplot(.,aes(x=X1,y=X2,label=Sample,colour=Genotype,shape=Diet)) + 
  geom_text_repel() + geom_point(size=3) + theme_bw() + xlab("") + ylab("") + ggtitle("t-sne Clustering")


#PCA with 1000 most variable genes
# something is...different about sample 7, we may want to remove it. Not certain. 
# also, we see SOME difference in PC1 (26% of variance between genotypes)
# I've checked PC3 and PC4 and not much interesting is happening there
library(matrixStats)
ntop = 1000
Pvars <- rowVars(assay(rld))
select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop, 
                                                      length(Pvars)))]
PCA <- prcomp(t(assay(rld)[select, ]), scale = F)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
dataGG = data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2], 
                    PC3 = PCA$x[,3], PC4 = PCA$x[,4], 
                    Diet = colData(rld)$Diet,
                    Genotype = colData(rld)$Genotype,
                    Samples=colData(rld)$sample.ID)

ggplot(dataGG, aes(PC1, PC2, label=Samples, color=Genotype, shape=Diet)) +
  geom_point(size=3) + geom_text_repel() +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + theme_bw()