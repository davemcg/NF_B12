library(ggplot2)
library(ggrepel)


volcano_plotter <- function(df, title) {
  df %>% 
    mutate(Significance=ifelse(padj<0.05 & abs(log2FoldChange)>0.25,"FDR<0.05 & log2FC>0.25",ifelse(padj<0.05,"FDR<0.05","Not Significant"))) %>% 
    mutate(Significance = factor(Significance,levels=c('Not Significant','FDR<0.05','FDR<0.05 & log2FC>0.25'))) %>% 
    ggplot(aes(x=log2FoldChange,y=-log10(pvalue),label=Gene)) + geom_point(aes(colour=Significance)) +
     geom_text_repel(data=subset(df, abs(log2FoldChange) > 0.25 & padj<0.05)) +
     scale_colour_manual(values=c("gray","pink","darkred")) +
    ggtitle(paste("Volcano Plot,",title)) + theme_bw()
}


#df$Significant <- ifelse(df$padj < 0.05, "FDR < 0.05", "Not Significant")
#df$Significant[df$padj < 0.05 & abs(df$log2FoldChange) > 2] <- "FDR < 0.05 & log2FC > 2"
#df$Significant <- factor(df$Significant,levels=c("Not Significant","FDR < 0.05", "FDR < 0.05 & log2FC > 2"))
#ggplot(data=df,aes(x=log2FoldChange,y=-log10(pvalue))) + 
#  geom_point(aes(colour=Significant)) +
#  scale_colour_manual(values=c("gray","pink","darkred")) + 
#  geom_text_repel(data=subset(df,abs(log2FoldChange) > 2 & padj<0.05), 
#                  aes(label=Gene_Name)) +
#  geom_hline(aes(yintercept=-log10(0.05/nrow(df))),linetype="dotted") +
#  geom_vline(aes(xintercept=-2),linetype="dotted") +
#  geom_vline(aes(xintercept=2),linetype="dotted") +
#  scale_x_continuous(breaks=c(-6,-4,-2,0,2,4,6)) +
#  ggtitle("Volcano Plot, Nlz2 KO vs WT") +
#  theme_bw()

#pdf(file="~/Documents/PROJECTS/brooks/nlz2-rna-seq/nlz2_volcano_DESeq_WT_vs_Hom.pdf",width=10,height=12)
#plot
#dev.off()