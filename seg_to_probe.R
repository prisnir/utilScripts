library(CNTools)
args = commandArgs(trailingOnly=TRUE)
#### FUNCTION #####
preProcCNA <- function(segfile, genebed, gain, amp, htz, hmz, genelist){
  
  ## small fix segmentation data
  segData <- read.delim(segfile, header=TRUE) # segmented data already
  # /Volumes/gsiprojects/dicklab/MATS/data/WG/sequenza/res
  #sample	chr	start	end	value
  names(segData) <- c("ID", "chrom", "loc.start", "loc.end", "seg.mean")
  segData <- segData[,c("chrom", "loc.start", "loc.end","seg.mean", "ID")]
  segData$chrom <- gsub("chr", "", segData$chrom)
  
  # thresholds
  print("setting thresholds")
  gain=as.numeric(gain)
  amp=as.numeric(amp)
  htz=as.numeric(htz)
  hmz=as.numeric(hmz)
  
  # get the gene info
  print("getting gene info")
  geneInfo <- read.delim(genebed, sep="\t", header=F)
  names(geneInfo) <- c("chrom", "start", "end", "gene_def")
  
  # make CN matrix gene level
  print("converting seg")
  cnseg <- CNSeg(segData)
  rdByGene <- getRS(cnseg, by="gene", imput=FALSE, XY=FALSE, geneMap=geneInfo, what="median")
  reducedseg <- rs(rdByGene)
  df_cna <- data.frame(reducedseg)
  # df_cna <- merge(df_cna, geneInfo, by=c("chrom", "start", "end")) 
  
  # # some reformatting and return log2cna data
  # df_cna <- subset(reducedseg[,c(5, 6:ncol(reducedseg))], !duplicated(reducedseg[,c(5, 6:ncol(reducedseg))][,1]))
  # colnames(df_cna) <- c("Hugo_Symbol", colnames(df_cna)[2:ncol(df_cna)])
  
  # set thresholds and return 5-state matrix
  print("thresholding cnas")
  df_cna_thresh <- df_cna
  df_cna_thresh[,c(5:ncol(df_cna))] <- sapply(df_cna_thresh[,c(5:ncol(df_cna))], as.numeric)
  
  # threshold data
  for (i in 5:ncol(df_cna_thresh))
  {
    df_cna_thresh[,i] <- ifelse(df_cna_thresh[,i] > amp, 2,
                                ifelse(df_cna_thresh[,i] < hmz, -2,
                                       ifelse(df_cna_thresh[,i] > gain & df_cna_thresh[,i] <= amp, 1,
                                              ifelse(df_cna_thresh[,i] < htz & df_cna_thresh[,i] >= hmz, -1, 0)
                                       )
                                )
    )
  }
  
  head(df_cna_thresh)
  
  # fix rownames of log2cna data
  # rownames(df_cna) <- df_cna$Hugo_Symbol
  # df_cna$Hugo_Symbol <- NULL
  # df_cna <- signif(df_cna, digits=4)
  
  # fix rownames of thresholded data
  # row.names(df_cna_thresh) <- df_cna_thresh[,1]
  # df_cna_thresh <- df_cna_thresh[,-1] # matrix where row names are genes, samples are columns
  # 
  # subset if gene list given
  # if (exists("genelist")) {
  #   keep_genes <- readLines(genelist)
  #   df_cna <- df_cna[row.names(df_cna) %in% keep_genes,]
  #   df_cna_thresh <- df_cna_thresh[row.names(df_cna_thresh) %in% keep_genes,]
  # }
  
  # return the list of dfs
  CNAs=list()
  CNAs[[1]] <- segData
  CNAs[[2]] <- df_cna
  CNAs[[3]] <- df_cna_thresh
  return(CNAs)
  
}


genebed="/Volumes/gsiprojects/dicklab/MATS/data/reference/intervals/variant_intervals_nochr.bed"
# format this 
genebed.data <- read.csv(genebed, sep = "\t", header =F)
names(genebed.data) <- c("chrom", "spos", "epos", "def")
head(genebed.data)

rep_vals <- data.frame(table(genebed.data$def))
rep_vals$new_def <- paste0(rep_vals$Var1, "_", rep_vals$Freq)
row.names(rep_vals) <- rep_vals$Var1

genebed.data$new_def <- rep_vals[genebed.data$def,]$new_def
genebed.data$chrom <- paste0("chr",genebed.data$chrom)
head(genebed.data)

# save this 
write.table(genebed.data[,c("chrom", "spos", "epos", "new_def")], 
            file = "/Volumes/gsiprojects/dicklab/MATS/data/reference/intervals/variant_intervals_chr_formatted.bed",
            sep = "\t", row.names = F, col.names=F, quote = F)

genebed <- "/Volumes/gsiprojects/dicklab/MATS/data/reference/intervals/variant_intervals_nochr_formatted.bed"
gain=0.3
amp=0.7
htz=-0.3
hmz=-0.7
genelist="/Volumes/gsiprojects/dicklab/MATS/data/reference/intervals/genes_of_interest.txt" # read targeted gene list

cbiodir <- "~/Projects/GSI/Projects/MATS/WG/CNV_summary"

# segfile=args[1] # read seg file
# ident=args[2] # identifier

segfile <- "/Volumes/gsiprojects/dicklab/MATS/data/WG/sequenza/res/sequenza_combination.seg" # read seg file
ident <- "MATS_WG" # identifier

CNAs <- preProcCNA(segfile, genebed, gain, ampl, htzd, hmzd, genelist)
# write cbio files
print("writing seg file")
write.table(CNAs[[1]], file=paste0(cbiodir, "/data_segments_", ident, ".txt"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(CNAs[[2]],
            file=paste0(cbiodir, "/data_log2CNA_", ident, ".txt"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(CNAs[[3]],
            file=paste0(cbiodir, "/data_CNA_", ident, ".txt"), sep="\t", row.names=FALSE, quote=FALSE)

####### SUBSET TO MATS08/MATS_0002 ############
rel_cols <- colnames(df_cna)[grepl("MATS_0002",colnames(df_cna))]
df_cna_MATS08 <- df_cna[,c("chrom", "start", "end", "def", rel_cols)]
df_cna_thresh_MATS08 <- df_cna_thresh[,c("chrom", "start", "end", "def", rel_cols)]

# write to file
write.table(df_cna_MATS08,
            file=paste0(cbiodir, "/data_log2CNA_", "MATS08", ".txt"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(df_cna_thresh_MATS08,
            file=paste0(cbiodir, "/data_CNA_", "MATS08", ".txt"), sep="\t", row.names=FALSE, quote=FALSE)

# 
####################################


