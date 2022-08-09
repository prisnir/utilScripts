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
    geneInfo <- read.delim(genebed, sep="\t", header=TRUE)
    
    # make CN matrix gene level
    print("converting seg")
    cnseg <- CNSeg(segData)
    rdByGene <- getRS(cnseg, by="gene", imput=FALSE, XY=FALSE, geneMap=geneInfo, what="median")
    reducedseg <- rs(rdByGene)
    
    # some reformatting and return log2cna data
    df_cna <- subset(reducedseg[,c(5, 6:ncol(reducedseg))], !duplicated(reducedseg[,c(5, 6:ncol(reducedseg))][,1]))
    colnames(df_cna) <- c("Hugo_Symbol", colnames(df_cna)[2:ncol(df_cna)])
    
    # set thresholds and return 5-state matrix
    print("thresholding cnas")
    df_cna_thresh <- df_cna
    df_cna_thresh[,c(2:ncol(df_cna))] <- sapply(df_cna_thresh[,c(2:ncol(df_cna))], as.numeric)
    
    # threshold data
    for (i in 2:ncol(df_cna_thresh))
    {
        df_cna_thresh[,i] <- ifelse(df_cna_thresh[,i] > amp, 2,
                                    ifelse(df_cna_thresh[,i] < hmz, -2,
                                           ifelse(df_cna_thresh[,i] > gain & df_cna_thresh[,i] <= amp, 1,
                                                  ifelse(df_cna_thresh[,i] < htz & df_cna_thresh[,i] >= hmz, -1, 0)
                                           )
                                    )
        )
    }
    
    # fix rownames of log2cna data
    rownames(df_cna) <- df_cna$Hugo_Symbol
    df_cna$Hugo_Symbol <- NULL
    df_cna <- signif(df_cna, digits=4)
    
    # fix rownames of thresholded data
    row.names(df_cna_thresh) <- df_cna_thresh[,1]
    df_cna_thresh <- df_cna_thresh[,-1] # matrix where row names are genes, samples are columns
    
    # subset if gene list given
    if (exists("genelist")) {
        keep_genes <- readLines(genelist)
        df_cna <- df_cna[row.names(df_cna) %in% keep_genes,]
        df_cna_thresh <- df_cna_thresh[row.names(df_cna_thresh) %in% keep_genes,]
    }
    
    # return the list of dfs
    CNAs=list()
    CNAs[[1]] <- segData
    CNAs[[2]] <- df_cna
    CNAs[[3]] <- df_cna_thresh
    return(CNAs)
    
}


genebed="/u/prath/gencode_v33_hg38_genes_nochrom.bed" # read appropriate gene to coordinated conversion file
gain=0.3
ampl=0.7
htzd=-0.3
hmzd=-0.7
genelist="/u/prath/targeted_genelist.txt" # read targeted gene list

cbiodir <- "/.mounts/labs/gsiprojects/dicklab/MATS/data/WG/sequenza/res"

segfile=args[1] # read seg file
ident=args[2] # identifier

CNAs <- preProcCNA(segfile, genebed, gain, ampl, htzd, hmzd, genelist)
# write cbio files
print("writing seg file")
write.table(CNAs[[1]], file=paste0(cbiodir, "/data_segments_", ident, ".txt"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(data.frame("Hugo_Symbol"=rownames(CNAs[[2]]), CNAs[[2]], check.names=FALSE),
            file=paste0(cbiodir, "/data_log2CNA_", ident, ".txt"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(data.frame("Hugo_Symbol"=rownames(CNAs[[3]]), CNAs[[3]], check.names=FALSE),
            file=paste0(cbiodir, "/data_CNA_", ident, ".txt"), sep="\t", row.names=FALSE, quote=FALSE)

# write the truncated data_CNA file (remove genes which are all zero) for oncoKB annotator
df_CNA <- CNAs[[3]][apply(CNAs[[3]], 1, function(row) !all(row == 0 )),]
write.table(data.frame("Hugo_Symbol"=rownames(df_CNA), df_CNA, check.names=FALSE),
            file=paste0(cbiodir, "/data_CNA_short_", ident, ".txt"), sep="\t", row.names=FALSE, quote=FALSE)

####################################


