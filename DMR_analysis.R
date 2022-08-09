# # annotate DMRs
# library(genomation)
# library(methylKit)
# library(AnnotationHub)
# library(GenomicFeatures) 
# library(CompGO)
# # checking out annotation hub
# ah <- AnnotationHub()
# ah <- subset(ah, species == "Homo sapiens" & genome=="GRCh38")
# ah
# # human <- query(ah, c("Homo sapiens"))
# # human
# grch38_granges <- query(ah, "GRanges")
# head(grch38_granges, 65)
# 
# # AH55233 | Homo_sapiens.GRCh38.89.chr.gtf 
# # AH68818 | Homo_sapiens.GRCh38.95.abinitio.gtf            
# # AH68819 | Homo_sapiens.GRCh38.95.chr.gtf                 
# # AH68820 | Homo_sapiens.GRCh38.95.chr_patch_hapl_scaff.gtf
# # AH68821 | Homo_sapiens.GRCh38.95.gtf 
# 
# hg38_gencode <- ah[['AH68819']]
# head(hg38_gencode, 3)

# dna2 <- ah[['AH49722']]
# for makeTxDbFromGRanges
# 
# txdb2 <- makeTxDbFromGRanges(hg38_gencode)
# # DMR.test$chr <- DMR.test
# bed.sample <- DMR.test[,c("chr", "start", "end")]
# bed.sample$chr <- gsub("chr", "", bed.sample$chr)
# range = GRanges(seqnames=bed.sample$chr, 
#                 IRanges(start=as.numeric(bed.sample$start), end=as.numeric(bed.sample$end)))
# x2 = annotateBedFromDb(gRanges = range, db = txdb2)

# read DMRs
library(methyAnalysis)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

fns <- list.files("~/test/DMR_topup//", pattern=".bed")
for (fn in fns){
  if (file.size(paste0("~/test/DMR_topup/",fn)) == 0){
    next
  }
  print(fn)
DMR.test <- read.csv(paste0("~/test/DMR_topup/",fn), 
                     sep = "\t", header = F)
head(DMR.test)
if (dim(DMR.test)[1] <= 1){
  next
}
DMR.test$chr <- stringr::str_split_fixed(DMR.test$V1, ":",2)[,1]
tmp.vr <- stringr::str_split_fixed(DMR.test$V1, ":",2)[,2]
DMR.test$start <- stringr::str_split_fixed(tmp.vr, "-", 2)[,1]
DMR.test$end <- stringr::str_split_fixed(tmp.vr, "-", 2)[,2]

####### example 2 SUCCESSFUL #####

bed.sample <- DMR.test[,c("chr", "start", "end")]
# bed.sample$chr <- gsub("chr", "", bed.sample$chr)
range = GRanges(seqnames=bed.sample$chr, 
                IRanges(start=as.numeric(bed.sample$start), end=as.numeric(bed.sample$end)))

sigDMRInfo.ann2 <- annotateDMRInfo(range, 'TxDb.Hsapiens.UCSC.hg38.knownGene')

data.annot <- data.frame(sigDMRInfo.ann2$sigDMRInfo)

write.table(data.annot, 
            file=paste0("~/test/DMR_topup/", fn, "_annotated.txt"),
            sep = "\t",
            quote = F,
            row.names = F)
}

# install.packages("pathfindR")

# ah2 = AnnotationHub()
# epiFiles <- query(ah2, "EpigenomeRoadMap")
# epiFiles

# #####
# library(tximport)
# library(readr)
# library(DESeq2)
# library(biomaRt)
# 
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# library(CompGO)
# txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
# bed.sample <- DMR.test[,c("chr", "start", "end")]
# range = GRanges(seqnames=bed.sample$chr, IRanges(start=as.numeric(bed.sample$start), end=as.numeric(bed.sample$end)))
# x = annotateBedFromDb(gRanges = range, db = txdb)
# 
# # making our own txdb obtject
# download.file("https://www.encodeproject.org/files/gencode.v24.primary_assembly.annotation/@@download/gencode.v24.primary_assembly.annotation.gtf.gz", "gencode.v24.primary_assembly.annotation.gtf.gz")
# library(GenomicFeatures)
# txdb = makeTxDbFromGFF('gencode.v24.primary_assembly.annotation.gtf.gz')
# y = annotateBedFromDb(gRanges = range, db = txdb)
# y
# # annotates the ENGS ids for genes 
# # hugo names for the genes can be obtained from other source/e.g. ENSG to HUGO map
# # data table for ENSG to HUGO
# ensble.hugo <- read.csv("~/Projects/GSI/Projects/CMPP/Annotation/cBioWrap/files/ensemble_conversion.txt", sep = "\t", header =F)
# row.names(ensble.hugo) <- ensble.hugo$V1
# 
# # y to data frame
# df <- data.frame(seqnames=seqnames(y),
#                  start=start(y),
#                  end=end(y),
#                  strand=strand(y),
#                  gene_id=y$gene_id)
# df$HUGO <- ensble.hugo[df$gene_id.value,c("V2")]
# 
# # df data frame contains all hugo names
# 
# ### METHOD 2 ####
# library(Organism.dplyr)
# library(org.Hs.eg.db)
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# 
# src <- src_organism("TxDb.Hsapiens.UCSC.hg38.knownGene")
# 
# src <- src_ucsc("Homo sapiens")
# 
# 
# gr <- GRangesFilter(GenomicRanges::GRanges("chr1:44000000-55000000"))
# transcripts(src, filter=~(symbol %startsWith% "SNORD" & gr) | symbol == "ADA")
# 
# 
# 
# # read the gene BED file
# gene.obj=readTranscriptFeatures(system.file("extdata", "refseq.hg38.bed.txt", 
#                                             package = "methylKit"))
# promoters=regionCounts(myobj,gene.obj$promoters)
# 
# head(promoters[[1]])
# # pathways analysis on DMRs
# 


# PATHWAYS

library(pathfindR)
library(GenomicRanges)
library(methyAnalysis)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
fns <- list.files("~/Projects/GSI/Projects/CMPP/DMR_topup/", pattern=".txt")
# fns <- list.files("~/test/DMR_topup/", pattern=".txt")
fns <- fns[!grepl(".bed_annotated.txt", fns)]
fns <- fns[!grepl("_rel_pth", fns)]
fns <- fns[grepl("topup3", fns)]
for (fn in fns){
# df.sig <- read.csv(paste0("~/test/DMR_topup/", fn),
# #                    sep = "\t", as.is = T)
# df.sig <- read.csv(paste0("~/Projects/GSI/Projects/CMPP/DMR_topup/", fn),
#                      sep = "\t", as.is = T)
# df.sig <- df.sig[df.sig$padj <= 0.05,] # establish significance
# if (is.na(unique(df.sig$log2FoldChange))){
#   next
# }
# # cleanse all NAs
# df.sig <- df.sig[! is.na(df.sig$log2FoldChange),]
# # df.sig 
# regs <- data.frame(row.names(df.sig))
# if (dim(regs)[1] < 1){
#   next
# }
# regs$chr <- stringr::str_split_fixed(regs$row.names.df.sig., ":",2)[,1]
# tmp.vr <- stringr::str_split_fixed(regs$row.names.df.sig., ":",2)[,2]
# regs$start <- stringr::str_split_fixed(tmp.vr, "-", 2)[,1]
# regs$end <- stringr::str_split_fixed(tmp.vr, "-", 2)[,2]
# 
# regs <- regs[!grepl('NA',regs$chr),]
# bed.sample <- regs[,c("chr", "start", "end")]
# # head(bed.sample)
# 
# range = GRanges(seqnames=bed.sample$chr,
#                 IRanges(start=as.numeric(bed.sample$start),
#                         end=as.numeric(bed.sample$end)))

# re-create cpm data like file and save it 
# sigDMRInfo.ann3 <- annotateDMRInfo(range, 'TxDb.Hsapiens.UCSC.hg38.knownGene')
# read this from the annotated.txt files

# data.annot <- data.frame(sigDMRInfo.ann3$sigDMRInfo)
# head(data.annot)
# 
# # reannotated the CPM object
# 
# data.annot$defintion <- paste0(data.annot$seqnames, ":",
#                                data.annot$start, "-", data.annot$end)
# 
# row.names(data.annot) <- data.annot$defintion

# df.sig$Gene_Symbol <- data.annot[row.names(df.sig),]$GeneSymbol
# # write.table(df.sig, "~/test/")
# print (head(df.sig))
# # break
# 
# colnames(df.sig) <- c("logFC", "p",	"FDR_p", "Gene_symbol")
# df.sig <- df.sig[,c("Gene_symbol", "logFC", "FDR_p")]
# df.sig <- df.sig[!is.na(df.sig$FDR_p),] # remove NAs from FDRs
fnd <- gsub(".txt", "", fn)
# write.table(df.sig, paste0("~/test/DMR_topup/", fnd, "_rel_pth.txt"),
#              sep= "\t", row.names = F, quote = F)
df.sig <- read.csv(paste0("~/Projects/GSI/Projects/CMPP/DMR_topup/", fnd, "_rel_pth.txt"), sep = "\t", as.is = T)
head(df.sig)
pdf(paste0("~/Projects/GSI/Projects/CMPP/Annotation/res2/", fnd, "_rel_pth.pdf"))
output_df <- run_pathfindR(df.sig,
                           output_dir = paste0("~/Projects/GSI/Projects/CMPP/Annotation/res2/", fnd),
                           iterations = 5)
dev.off()
pdf(paste0("~/Projects/GSI/Projects/CMPP/Annotation/res2/", fnd, "_clustered.pdf"))
clustered_df <- cluster_enriched_terms(output_df)
dev.off()
}


