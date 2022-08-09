library(ggplot2)
library(reshape2)

sigs.file <- "/Volumes/gsiprojects/external/zadehglab/MPNST/data/MPNST/cbioportal/cBioWrap_20190909/output/supplementary_data/sigs/weights.txt"

sigs.data <- read.csv(sigs.file, sep = "\t", as.is = T)

head(sigs.data)


mastersheet <-  read.csv("/Volumes/gsiprojects/external/zadehglab/MPNST/mastersheet/MPNST.mastersheet.csv", as.is = T)
head(mastersheet)


sigs.data <- merge(x = sigs.data, y = mastersheet, by.x = "X", by.y = "external_name")
head(sigs.data)

rel.colnames <- c("X",
                  "Signature.1A",
                  "Signature.1B",
                  "Signature.2",
                  "Signature.3",
                  "Signature.4",
                  "Signature.5",
                  "Signature.6",
                  "Signature.7",
                  "Signature.8",
                  "Signature.9",
                  "Signature.10",
                  "Signature.11",
                  "Signature.12",
                  "Signature.13",
                  "Signature.14",
                  "Signature.15",
                  "Signature.16",
                  "Signature.17",
                  "Signature.18",
                  "Signature.19",
                  "Signature.20",
                  "Signature.21",
                  "Signature.R1",
                  "Signature.R2",
                  "Signature.R3",
                  "Signature.U1",
                  "Signature.U2" ,
                  "proj_id",
                  "miso_id",
                  "patient_name",
                  "matched_miso",
                  "matched_external")

melt.id <- c("X",
             "proj_id",
             "miso_id",
             "patient_name",
             "matched_miso",
             "matched_external")

signature_Def = c(
  # "Signature.1" = "spontaneous deamination of 5-methylcytosine",
  "Signature.1A" = "spontaneous deamination of 5-methylcytosine",
  "Signature.1B" = "spontaneous deamination of 5-methylcytosine",
  "Signature.2" = "activity of the AID/APOBEC family of cytidine deaminases",
  "Signature.3" = "failure of DNA double-strand break-repair by homologous recombination",
  "Signature.4" = "smoking and exposure to tobacco carcinogens",
  "Signature.5" = "unknown",
  "Signature.6" = "defective DNA mismatch repair and microsatellite instability",
  "Signature.7" = "ultraviolet light exposure",
  "Signature.8" = "unknown",
  "Signature.9" = "activity of AID during somatic hypermutation",
  "Signature.10" = "recurrent POLE somatic mutations",
  "Signature.11" = "alkylating agents",
  "Signature.12" = "unknown",
  "Signature.13" = "activity of the AID/APOBEC family of cytidine deaminases converting cytosine to uracil",
  "Signature.14" = "unknown",
  "Signature.15" = "defective DNA mismatch repair",
  "Signature.16" = "unknown",
  "Signature.17" = "unknown",
  "Signature.18" = "unknown",
  "Signature.19" = "unknown",
  "Signature.20" = "defective DNA mismatch repair",
  "Signature.21" = "unknown",
  "Signature.22" = "exposures to aristolochic acid",
  "Signature.23" = "unknown",
  "Signature.24" = "exposures to aflatoxin",
  "Signature.25" = "unknown",
  "Signature.26" = "defective DNA mismatch repair",
  "Signature.27" = "unknown",
  "Signature.28" = "unknown",
  "Signature.29" = "tobacco chewing",
  "Signature.30" = "unknown",
  "Signature.R1" = "unknown",
  "Signature.R2" = "unknown",
  "Signature.R3" = "unknown",
  "Signature.U1" = "unknown",
  "Signature.U2" = "unknown"
)

melt.sigs.data <- melt(sigs.data, id = melt.id)
melt.sigs.data$matched_external <- ifelse(is.na(melt.sigs.data$matched_external), "unmatched", melt.sigs.data$matched_external)
melt.sigs.data$matched_status <- ifelse(melt.sigs.data$matched_external == "unmatched", "unmatched", "matched")
melt.sigs.data$value <- as.numeric(melt.sigs.data$value)
melt.sigs.data <- melt.sigs.data[melt.sigs.data$variable %in% names(signature_Def), ]
melt.sigs.data$Description <- paste(melt.sigs.data$variable, signature_Def[melt.sigs.data$variable], sep = ":")
head(melt.sigs.data)



library(RColorBrewer)

colourCount = length(unique(melt.sigs.data$Description))
getPalette = colorRampPalette(brewer.pal(9, "Set3"))
clrs <- getPalette(colourCount)
# creat color maps
cmap <- c()
for (i in 1:colourCount){
  sgn <- unique(melt.sigs.data$Description)[[i]]
  clr <- clrs[[i]]
  # print (c(as.character(sgn), clr))
  cmap[as.character(sgn)] <- clr
}

melt.sigs.data$Description <- as.character(melt.sigs.data$Description)

# also add clinical info
melt.sigs.data <- merge(melt.sigs.data, clinical.data, by.x = "X", by.y= "Sample_ID")
melt.sigs.data$annot <-ifelse(melt.sigs.data$Tumor.Type == "Neurofibroma", "NF", ifelse(melt.sigs.data$Tumor.Type == "Premalignant_NF", "pre-NF", "MPNST"))
  
colScale = scale_colour_manual(name = "Mutational Signatures",values = cmap, aesthetics = "fill")


p <- ggplot(melt.sigs.data[!melt.sigs.data$variable %in% c("Signature.1A", "Signature.1B"),]) + geom_bar(aes(x = X, y = value, fill = Description), stat = "identity", position = "stack") + 
  facet_grid(.~matched_status+annot, space = "free", scales = "free") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  xlab("Sample IDs") + ylab("Relative strength of COSMIC mutational signatures") + colScale
p

ggsave("~/Documents/GSI/MPSNT/Mutational_sigs/MutSigs2.png", p, width = 20, height = 5)


p <- ggplot(melt.sigs.data) + geom_bar(aes(x = X, y = value, fill = Description), stat = "identity", position = "stack") + 
  facet_grid(.~matched_status+annot, space = "free", scales = "free") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  xlab("Sample IDs") + ylab("Relative strength of COSMIC mutational signatures") + colScale
p

ggsave("~/Documents/GSI/MPSNT/Mutational_sigs/MutSigs.png", p, width = 20, height = 5)


# MISC
matched.samples <- unique(melt.sigs.data[melt.sigs.data$matched_status == "matched", ]$X)

unmatched.samples <- unique(melt.sigs.data[melt.sigs.data$matched_status == "unmatched", ]$X)

data.segs <- read.csv("/Volumes/gsiprojects/external/zadehglab/MPNST/data/MPNST/cbioportal/cBioWrap_20190909/output/cbioportal_import_data/data_segments.txt",
                      sep = "\t", as.is = T)

data.segs.unmatched <- data.segs[data.segs$ID %in% unmatched.samples,]

data.segs.matched <- data.segs[data.segs$ID %in% matched.samples,]

write.table(data.segs.matched, "~/Documents/GSI/MPSNT/copy_numbers/MPNST.matched.seg", sep = "\t", quote = F, row.names = F)

write.table(data.segs.unmatched, "~/Documents/GSI/MPSNT/copy_numbers/MPNST.unmatched.seg", sep = "\t", quote = F, row.names = F)


