# kallisto import via tximport to DESeq2
library(tximport)
#  tximport_1.0.3
library(DESeq2)
# DESeq2_1.12.3 
library(data.table)
library(readr)
library(dplyr)

# grab all kallisto quantification files
# currently technical replicates are separate
files <- paste(list.files(path='/Volumes/ThunderBay/PROJECTS/brody/NF_B12_mouse/kallisto/', pattern=".*.*_kallisto", full.names=TRUE, include.dirs=TRUE, recursive=TRUE), '/abundance.tsv', sep='')

# annotation info to jump from transcript to gene
# kallisto quantifies on the transcript level, but it is much 
# cleaner to work on the gene level, statistically and mentally
# the mouse fasta file has the file info in it, so I can parse it
# to grab the gene name
anno <- fread(files[1]) # any file will work
anno$Gene <- sapply(anno$target_id,function(x) strsplit(x,'\\|')[[1]][6])

# actually merge tx specific counts to gene level
txi <- tximport(files, type = "kallisto", tx2gene = anno, reader = read_tsv, countsFromAbundance = c("lengthScaledTPM"))
TPM <- data.frame(txi$counts)

# name again with folder names
names <- sapply(files, function(x) strsplit(x,"\\/|_kallisto")[[1]][9])
colnames(TPM) <- names

######
# quick side-step into clustering with t-sne
library(ggrepel)
library(Rtsne)
tsne_out <- Rtsne(as.matrix(log2(t(TPM)+1)),perplexity = 7)
tsne_plot <- data.frame(tsne_out$Y)
tsne_plot$Run <- colnames(TPM)

decoder <- fread('NF_12.decoder.matrix')
tsne_plot <- left_join(tsne_plot,decoder,by=c("Run"="unique.ID"))
ggplot(tsne_plot,aes(x=X1,y=X2,label=sample.ID,colour=sample.ID,shape=lane.ID)) + 
  geom_text_repel() + geom_point(size=3) + theme_bw() + xlab("") + ylab("") + ggtitle("t-sne Clustering")
#######
