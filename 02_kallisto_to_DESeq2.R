# kallisto import via tximport to DESeq2
library(tximport)
#  tximport_1.0.3
library(DESeq2)
# DESeq2_1.12.3 
library(data.table)
library(readr)
library(dplyr)
library(readxl)

# grab all kallisto quantification files
# currently technical replicates are separate
files <- paste(list.files(path='/Volumes/ThunderBay/PROJECTS/brody/NF_B12_mouse/kallisto', pattern=".*.*_kallisto", full.names=TRUE, include.dirs=TRUE, recursive=TRUE), '/abundance.tsv', sep='')

# Annotation info to jump from transcript to gene.
# Kallisto quantifies on the transcript level, but it is much 
# cleaner to work on the gene level, statistically and mentally.
# The mouse fasta file has the file info in it, so I can parse it
# to grab the gene name
anno <- fread(files[1]) # any file will work
anno$Gene <- sapply(anno$target_id,function(x) strsplit(x,'\\|')[[1]][6])
anno <- data.frame(anno)
anno <- anno %>% dplyr::select(target_id, Gene)

# actually merge tx specific counts to gene level
txi <- tximport(files, type = "kallisto", tx2gene = anno, reader = read_tsv, countsFromAbundance = c("lengthScaledTPM"))
TPM <- data.frame(txi$counts)

# name again with folder names
names <- sapply(files, function(x) strsplit(x,"\\/|_kallisto")[[1]][8])
colnames(TPM) <- names

######
# quick side-step into clustering with t-sne
# library(ggrepel)
#library(Rtsne)
#tsne_out <- Rtsne(as.matrix(log2(t(TPM)+1)),perplexity = 15)
#tsne_plot <- data.frame(tsne_out$Y)
#tsne_plot$Run <- colnames(TPM)

decoder <- fread('NF_12.decoder.matrix')
#tsne_plot <- left_join(tsne_plot,decoder,by=c("Run"="unique.ID"))
#ggplot(tsne_plot,aes(x=X1,y=X2,label=sample.ID,colour=sample.ID,shape=lane.ID)) + 
#  geom_text_repel() + geom_point(size=3) + theme_bw() + xlab("") + ylab("") + ggtitle("t-sne Clustering")
#######


# OK, business time. Time to send to prep for DESeq2
bernard_anno <- read_excel('~/git/NF_B12/key_for_randomization.xlsx')
bernard_anno$sample.ID <- sprintf("B12_%02d",as.integer(sapply(bernard_anno$`Sample #`, function(x) strsplit(x,'-')[[1]][2])))
bernard_anno <- bernard_anno %>% mutate(Diet = ifelse(grepl('B12',Diet),'Normal','B12_Deficient')) %>% dplyr::select(sample.ID, Genotype, Diet, Hematocrit, Age, Replicate)
SampleTable <- left_join(decoder,bernard_anno) %>% mutate(GenoDiet = paste(Genotype,Diet,sep="."))
row.names(SampleTable)<-SampleTable$unique.ID
txi.deseq2 <- tximport(files, type = "kallisto", tx2gene = anno, reader = read_tsv)
colnames(txi.deseq2$counts) <- names

ddsG <- DESeqDataSetFromTximport(txi, SampleTable, ~Genotype)
ddsD <- DESeqDataSetFromTximport(txi, SampleTable, ~Diet)
ddsMF <- DESeqDataSetFromTximport(txi, SampleTable, ~GenoDiet)




