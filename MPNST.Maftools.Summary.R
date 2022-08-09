
require(maftools)
library(ggplot2)
library(reshape2)



cbiofolder <- "/Users/prath/Documents/GSI/MPSNT/cbioportal/cBioWrap_20190909/"
muts.File <- paste0(cbiofolder, "data_mutations_extended.txt")
muts.Data <- read.csv(muts.File, header = TRUE, sep = "\t", as.is = T)
head(muts.Data)

# get VAF distribution for the project
ggplot(muts.Data) + 
  geom_density(aes(x=tumor_vaf), fill = "grey") + 
  geom_vline(xintercept = 0.05, color = "red", linetype = 2) + 
  geom_vline(xintercept = mean(muts.Data$tumor_vaf), linetype = 2, color = "blue") + 
  geom_vline(xintercept = median(muts.Data$tumor_vaf), linetype = 2, color = "green") + theme_bw() 

clin.File <- paste0(cbiofolder, "data_clinical_samples.txt")
clin.Data <- read.csv(clin.File, as.is = T, skip = 4, sep = "\t")
names(clin.Data)[2] <- "Tumor_Sample_Barcode"
clin.Data$MATCHED_CONTROL <- ifelse(is.na(clin.Data$PURITY), "No", "Yes")
head(clin.Data)


laml <- read.maf(muts.Data, clinicalData = clin.Data)
so.my.top.20.genes.are <- getGeneSummary(laml)$Hugo_Symbol[1:20]

# getting CNA for the same genes
cna.File <-  paste0(cbiofolder, "data_CNA.txt")
cna.Data <- read.csv(cna.File, header = TRUE, sep = "\t", as.is = T)
cna.Data <- cna.Data[cna.Data$Hugo_Symbol %in% so.my.top.20.genes.are,]
cna.Data <- melt(cna.Data)
cna.Data <- cna.Data[cna.Data$value %in% c(-2,2),]
cna.Data$CNV <- ifelse(cna.Data$value == -2, "Del", "Amp")
cna.Data <- cna.Data[,c("Hugo_Symbol", "variable", "CNV")]
names(cna.Data)[2] <- "Tumor_Sample_Barcode"
head(cna.Data)


laml <- read.maf(muts.Data, clinicalData = clin.Data, cnTable = cna.Data)

dev.off()
oncoplot(laml,
         showTumorSampleBarcodes = TRUE, 
         removeNonMutated = FALSE,
         top = 20, 
         clinicalFeatures = c("MATCHED_CONTROL","NF1_STATUS", "LOCATION"),
         sortByAnnotation = TRUE)


########################################
# separate analysis for matched cases

matched.cases <- clin.Data[clin.Data$MATCHED_CONTROL == "Yes",]$Tumor_Sample_Barcode

matched.muts.Data <- muts.Data[muts.Data$Tumor_Sample_Barcode %in% matched.cases,]
head(matched.muts.Data)

# get VAF distribution for the project
ggplot(matched.muts.Data) + 
  geom_density(aes(x=tumor_vaf), fill = "grey") + 
  geom_vline(xintercept = 0.05, color = "red", linetype = 2) + 
  geom_vline(xintercept = mean(muts.Data$tumor_vaf), linetype = 2, color = "blue") + 
  geom_vline(xintercept = median(muts.Data$tumor_vaf), linetype = 2, color = "green") + theme_bw() 

matched.clin.Data <- clin.Data[clin.Data$Tumor_Sample_Barcode %in% matched.cases,]
head(matched.clin.Data)


laml <- read.maf(matched.muts.Data, clinicalData = matched.clin.Data)
so.my.top.20.genes.are <- getGeneSummary(laml)$Hugo_Symbol[1:20]

# getting CNA for the same genes
cna.File <-  paste0(cbiofolder, "data_CNA.txt")
cna.Data <- read.csv(cna.File, header = TRUE, sep = "\t", as.is = T)
cna.Data <- cna.Data[cna.Data$Hugo_Symbol %in% so.my.top.20.genes.are, c("Hugo_Symbol", matched.cases)]
cna.Data <- melt(cna.Data)
cna.Data <- cna.Data[cna.Data$value %in% c(-2,2),]
cna.Data$CNV <- ifelse(cna.Data$value == -2, "Del", "Amp")
cna.Data <- cna.Data[,c("Hugo_Symbol", "variable", "CNV")]
names(cna.Data)[2] <- "Tumor_Sample_Barcode"
head(cna.Data)


laml <- read.maf(matched.muts.Data, clinicalData = matched.clin.Data, cnTable = cna.Data)

dev.off()
oncoplot(laml,
         showTumorSampleBarcodes = TRUE, 
         removeNonMutated = FALSE,
         top = 20, 
         clinicalFeatures = c("MATCHED_CONTROL","NF1_STATUS", "LOCATION"),
         sortByAnnotation = TRUE)


########################################
# separate analysis for unmatched cases
unmatched.cases <- clin.Data[clin.Data$MATCHED_CONTROL == "No",]$Tumor_Sample_Barcode

unmatched.muts.Data <- muts.Data[muts.Data$Tumor_Sample_Barcode %in% unmatched.cases,]
head(unmatched.muts.Data)

# get VAF distribution for the project
ggplot(unmatched.muts.Data) + 
  geom_density(aes(x=tumor_vaf), fill = "grey") + 
  geom_vline(xintercept = 0.05, color = "red", linetype = 2) + 
  geom_vline(xintercept = mean(muts.Data$tumor_vaf), linetype = 2, color = "blue") + 
  geom_vline(xintercept = median(muts.Data$tumor_vaf), linetype = 2, color = "green") + theme_bw() 

unmatched.clin.Data <- clin.Data[clin.Data$Tumor_Sample_Barcode %in% unmatched.cases,]
head(unmatched.clin.Data)


laml <- read.maf(unmatched.muts.Data, clinicalData = unmatched.clin.Data)
so.my.top.20.genes.are <- getGeneSummary(laml)$Hugo_Symbol[1:20]

# getting CNA for the same genes
cna.File <-  paste0(cbiofolder, "data_CNA.txt")
cna.Data <- read.csv(cna.File, header = TRUE, sep = "\t", as.is = T)
cna.Data <- cna.Data[cna.Data$Hugo_Symbol %in% so.my.top.20.genes.are, c("Hugo_Symbol", unmatched.cases)]
cna.Data <- melt(cna.Data)
cna.Data <- cna.Data[cna.Data$value %in% c(-2,2),]
cna.Data$CNV <- ifelse(cna.Data$value == -2, "Del", "Amp")
cna.Data <- cna.Data[,c("Hugo_Symbol", "variable", "CNV")]
names(cna.Data)[2] <- "Tumor_Sample_Barcode"
head(cna.Data)


laml <- read.maf(unmatched.muts.Data, clinicalData = unmatched.clin.Data, cnTable = cna.Data)

dev.off()
oncoplot(laml,
         showTumorSampleBarcodes = TRUE, 
         removeNonMutated = FALSE,
         top = 20, 
         clinicalFeatures = c("MATCHED_CONTROL","NF1_STATUS", "LOCATION"),
         sortByAnnotation = TRUE)


muts.Data[muts.Data$Hugo_Symbol == "NF1",c("Tumor_Sample_Barcode", "Hugo_Symbol", "HGVSp_Short" )]




