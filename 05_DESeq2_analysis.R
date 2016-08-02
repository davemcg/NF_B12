source('~/git/NF_B12/03_kallisto_to_DESeq2.R')

library(IHW)
library(fdrtool)

# collapse technical replicates (G is by genotype, d is by diet, MF is by genotype and diet (multi-factorial))
ddsCollG <- collapseReplicates(ddsG, ddsG$sample.ID, ddsG$unique.ID)
ddsCollD <- collapseReplicates(ddsD, ddsD$sample.ID, ddsD$unique.ID)
ddsCollMF <- collapseReplicates(ddsMF, ddsMF$sample.ID, ddsMF$unique.ID)

process <- function(dds,contrast, design){
  DESeq2Table <- DESeq(dds)
  DESeq2Table <- estimateDispersions(DESeq2Table)
  DESeq2Table <-  nbinomWaldTest(DESeq2Table)
  DESeq2Res <- results(DESeq2Table, contrast=contrast,design = ~ design,filterFun=ihw)

  print(table(DESeq2Res$padj < 0.1))
  # 101 initially

  # hist of pvalues
  hist(DESeq2Res$pvalue, col = "lavender", main = "WT vs KO, no correction", xlab = "p-values")

  # Filter by padj and padjust to get rid of NA
  DESeq2Res <- DESeq2Res[ !is.na(DESeq2Res$padj), ]
  DESeq2Res <- DESeq2Res[ !is.na(DESeq2Res$pvalue), ]

  # Remove adjusted values, to add new values later
  DESeq2Res <- DESeq2Res[, -which(names(DESeq2Res) == "padj")]
  # Use fdtool's zscore to restimate adjusted pvals
  FDR.DESeq2Res <- fdrtool(DESeq2Res$stat, statistic= "normal", plot = T)
  # add new pvalues back
  DESeq2Res[,"padj"]  <- p.adjust(FDR.DESeq2Res$pval, method = "BH")

  hist(FDR.DESeq2Res$pval, col = "royalblue4", 
     main = "WT vs KO, corrected null model", xlab = "CORRECTED p-values")

  table(DESeq2Res[,"padj"] < 0.1)

  # MA plot
  plotMA(DESeq2Res)

  # output data
  stats <- data.frame(DESeq2Res)
  gene_counts <- counts(DESeq2Table, normalized=TRUE)
  # merge together
  out <- merge(stats,gene_counts,by="row.names")
  colnames(out)[1] <- 'Gene'
  return(out)
}

cd320.ko_vs_ctrl<-process(ddsCollG,c('Genotype','KO','Control'),'~Genotype')
B12neg_vs_B12pos<-process(ddsCollD,c('Diet','B12_Deficient','Normal'),'~Diet')
cd320.koB12neg_vs_ctrlB12pos<-process(ddsCollMF,c('GenoDiet','KO.B12_Deficient','Control.Normal'),'~GenoDiet')
cd320.koB12neg_vs_ctrlB12neg<-process(ddsCollMF,c('GenoDiet','KO.B12_Deficient','Control.B12_Deficient'),'~GenoDiet')

#add gtf annotation
source('gtf_annotate.R')
cd320.ko_vs_ctrl <- gtf_merge(cd320.ko_vs_ctrl)
B12neg_vs_B12pos <- gtf_merge(B12neg_vs_B12pos)
cd320.koB12neg_vs_ctrlB12pos <- gtf_merge(cd320.koB12neg_vs_ctrlB12pos)
cd320.koB12neg_vs_ctrlB12neg <- gtf_merge(cd320.koB12neg_vs_ctrlB12neg)

# rename and reorder counts columns (B12_01 to B12_16)
renamer <- function(input) {
  colnames(input)[grep('B12',colnames(input))] <- 
  SampleTable %>% dplyr::select(sample.ID,Genotype,Diet) %>% distinct() %>% mutate(ID=paste(Genotype,Diet,sample.ID,sep='.')) %>% arrange(sample.ID) %>% .[['ID']]
  ordered_names <- SampleTable %>% dplyr::select(sample.ID,Genotype,Diet) %>% distinct() %>% mutate(ID=paste(Genotype,Diet,sample.ID,sep='.')) %>% arrange(ID) %>% .[['ID']]
  ordered_names <- c('Gene','Source','Chr','Start','End','baseMean','log2FoldChange','lfcSE','stat','pvalue','weight','padj',ordered_names)
  input %>% dplyr::select_(.dots=ordered_names)
}
cd320.ko_vs_ctrl <- renamer(cd320.ko_vs_ctrl)
B12neg_vs_B12pos <- renamer(B12neg_vs_B12pos)
cd320.koB12neg_vs_ctrlB12pos <- renamer(cd320.koB12neg_vs_ctrlB12pos)
cd320.koB12neg_vs_ctrlB12neg <- renamer(cd320.koB12neg_vs_ctrlB12neg)

# make volcano plots

