# preprocessing the Rcounts into a gene (hugo) x samples matrix
inRSEMCounts <- "~/Documents/TGL/TGL_support_project/TGL25/RNASEQ/Rcounts/RawCounts_Schwannomatosis_RNAseq.txt"

loadRCounts <- read.csv(inRSEMCounts, sep = "\t")
loadRCounts$ID <- as.character(loadRCounts$ID)
# loadRCounts$gene_id <- sapply(strsplit(loadRCounts$gene_id,"[.]"), `[`, 1)

#
hugoENSmap <- "~/Documents/TGL/TGL_support_project/TGL25/RNASEQ/ensemble_conversion.txt"
loadGeneMap <- read.csv(hugoENSmap, sep = "\t", header = F)
names(loadGeneMap) <- c("Ensembl.ID", "Gene.Symbol")
loadGeneMap$formatted.Ensembl.ID <- sapply(strsplit(as.character(loadGeneMap$Ensembl.ID),"[.]"), "[", 1)
head(loadGeneMap)
# loadGeneMap$Ensembl.ID <- as.character(loadGeneMap$Ensembl.ID)
# map gene_id to Ensembl.ID

write.table(loadGeneMap, "/Volumes/TGL/gsi/pipeline/data/TGL25/RNASEQ/ESTIMATE/ensGenes_2.txt", quote = F, col.names = F, row.names = F, sep = "\t")

rCountsGene <- merge(y = loadRCounts, x = loadGeneMap, by.y = "ID", by.x = "formatted.Ensembl.ID")
drops <- c("formatted.Ensembl.ID",  "Gene.Symbol")
rCountsGene <- rCountsGene[,!names(rCountsGene) %in% drops]

rCountsGene <- rCountsGene[!duplicated(rCountsGene),]
# rCountsGene$Gene.Symbol <- paste("\"", rCountsGene$Gene.Symbol, "\"", sep = "")
# names(rCountsGene) <- paste("\"", names(rCountsGene), "\"", sep = "")
head(rCountsGene)
names(rCountsGene)[1] <- "gene_id"
write.table(rCountsGene, file = "~/Documents/TGL/TGL_support_project/TGL25/RNASEQ/Rcounts/RCounts_formatted.txt", sep = "\t", 
            quote = F, row.names = F, col.names = T)
