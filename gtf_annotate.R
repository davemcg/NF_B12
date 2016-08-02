library(data.table)
library(dplyr)

gtf <- fread("~/GenomicData/gencode.vM8.primary_assembly.annotation.gtf")

gtf$Info <- sapply(gtf$V9, function(x) strsplit(x,";")[[1]])
gtf$Gene <- sapply(gtf$Info, function(x) {unlist(x)[grepl("gene_name",unlist(x))]})
gtf$Gene <- sapply(gtf$Gene, function(x) strsplit(x,"\"")[[1]][2])

gtf$Gene.ID <- sapply(gtf$Info, function(x) {unlist(x)[grepl("gene_id",unlist(x))]})
gtf$Gene.ID <- sapply(gtf$Gene.ID, function(x) strsplit(x,"\"")[[1]][2])

#gtf$Transcript.ID <- sapply(gtf$Info, function(x) {unlist(x)[grepl("transcript_id",unlist(x))]})



#gtf_core <- gtf %>% dplyr::filter(V3=="transcript") %>% select(V1,V4,V5, V2, Gene, Gene.ID, Transcript.ID)
gtf_core <- gtf %>% dplyr::filter(V3=="gene") %>% dplyr::select(V1,V4,V5, V2, Gene, Gene.ID)

#gtf_core$Transcript.ID <- sapply(gtf_core$Transcript.ID, function(x) strsplit(x,"\"")[[1]][2])
colnames(gtf_core) <- c("Chr","Start","End", "Source","Gene_Name","Gene_ID")

# edit to match whatever comparison you are doing
gtf_merge <- function(input){
  out <- merge(input,gtf_core,by.x="Gene",by.y="Gene_Name")
  out %>% dplyr::select(Gene, Source, Chr, Start, End, baseMean:B12_16)
}
