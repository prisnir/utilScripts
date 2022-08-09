library(DESeq2)

groupName <- "group3"
outputPath <- "/.mounts/labs/gsiprojects/external/svhaAU/RICH/Analysis/diffExpr/"



outputPath <- paste0(outputPath, "/", groupName)
dir.create(outputPath, showWarnings = F, recursive = T)

ens_id_map <- read.csv("/u/prath/cBioWrap/files/ensemble_conversion.txt", sep = "\t", header = F)
row.names(ens_id_map) <- ens_id_map$V1
ens_id_map

countMat <- read.csv("/.mounts/labs/gsiprojects/external/svhaAU/RICH/rsem/RICH2.genes_all_samples_COUNT.txt", sep = "\t")
row.names(countMat) <- countMat$gene_id
countMat <- countMat[,-1]



# annotations
annotDF <- read.csv("/.mounts/labs/gsiprojects/external/svhaAU/RICH/rich_master_annot.txt", sep = "\t")
row.names(annotDF) <- gsub("-", ".", annotDF$internal.sample.name)
annotDF[,groupName] <- gsub(" ", "_", annotDF[,groupName])
annotDF[,groupName] <- gsub("[+]", "_plus_", annotDF[,groupName])
annotDF

analysesConditions <- annotDF[colnames(countMat),]
analysesConditions$condition <- analysesConditions[,groupName]
# analysesConditions <- analysesConditions[analysesConditions$group1 %in% condition1,]
# subset count matrix
countData <- data.matrix(countMat)
# countData <- countData[,colnames(countData) %in% row.names(analysesConditions)]
# # refine the annotations DF
# analysesConditions <- data.frame(condition = analysesConditions[colnames(countData),"group1"])
# rownames(analysesConditions) <- colnames(countData)

analysesConditions

dds <- DESeqDataSetFromMatrix(countData=round(countData), colData=analysesConditions, design= ~ condition)

dds <- dds[ rowSums(counts(dds)) >= 10, ]


dds2 <- DESeq(dds)

# for various groups of conditions
for (cond1 in unique(analysesConditions$condition)){
  for (cond2 in unique(analysesConditions$condition)){
    if (cond1 == cond2){ next }
    contrast <- c("condition", cond1, cond2)
    
    conditions_text <- paste0(cond1, "_VS_", cond2)
    fname=paste0(outputPath, "/resOrdered_SelectedStats_", conditions_text, ".txt")
    conditions_text2 <- paste0(cond2, "_VS_", cond1)
    fname2=paste0(outputPath, "/resOrdered_SelectedStats_", conditions_text2, ".txt")
    if (file.exists(fname2)){
      next
    }
    # calculate the res and resOrdered
    res <- results(dds2, contrast=contrast)
    resOrdered <- res[order(res$padj),]
    # add HUGO gene symbols
    resOrdered$geneSymbol <- ens_id_map[row.names(resOrdered),]$V2
    # resOrdered.sig <- resOrdered[resOrdered$padj < 0.05,]
    
    
    
    write.table(data.frame(resOrdered[,c("log2FoldChange", "pvalue", "padj")]) , 
                file=fname, 
                col.names = T, row.names = T, sep="\t", append=F, quote = F)
    
    fname=paste0(outputPath, "/resOrdered_SelectedStats_significant_", conditions_text, ".txt")
    rs <- data.frame(resOrdered)
    resOrdered.sig <- rs[rs$padj < 0.05,]
    write.table(data.frame(resOrdered.sig[,c("log2FoldChange", "pvalue", "padj")]) , 
                file=fname, 
                col.names = T, row.names = T, sep="\t", append=F, quote = F)
    
    # for the hugo genes
    fname=paste0(outputPath, "/resOrdered_SelectedStats_HUGO_", conditions_text, ".txt")
    resOrdered <- resOrdered[!is.na(resOrdered$geneSymbol),]
    write.table(data.frame(resOrdered[,c("geneSymbol", "log2FoldChange", "pvalue", "padj")]) , 
                file=fname, 
                col.names = T, row.names = T, sep="\t", append=F, quote = F)
    rs <- data.frame(resOrdered)
    rs <- rs[!is.na(rs),]
    resOrdered.sig <- rs[rs$padj < 0.05,]
    fname=paste0(outputPath, "/resOrdered_SelectedStats_significant_HUGO_", conditions_text, ".txt")
    write.table(data.frame(resOrdered.sig[,c("geneSymbol", "log2FoldChange", "pvalue", "padj")]) , 
                file=fname, 
                col.names = T, row.names = T, sep="\t", append=F, quote = F)
  }
  # break
}



