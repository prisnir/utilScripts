# to obtain TMB for this project
library(dplyr)

allowed.cancer.types <- c("angs", 
                          "cscc",
                          "desm",
                          "es",
                          "gbm",
                          "glioma",
                          "lgg",
                          "mbl",
                          "mpnst",
                          "nbl",
                          "pcpg",
                          "rms",
                          "sarc",
                          "skcm",
                          "ucs",
                          "um",
                          "uvm")

# function
# get sample ids 
computePGA <- function(mpnst.seg.data, thresh = 0.3){
  mpnst.seg.data <- mpnst.seg.data[,c("tcga_ID", "additional_desc", "ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean")]
  PGA <- c("tcga_ID", "additional_desc", "Sample_ID", "num_segs", "genome_length")
  samples <- unique(mpnst.seg.data$ID)
  for (sample in samples){
    sample.seg <- mpnst.seg.data[mpnst.seg.data$ID == sample,]
    length.altered_any <- 0
    length.altered_high <- 0
    length.altered_LOH <- 0
    genome.length <- 0
    for (chrm in unique(sample.seg$chrom)){
      chrom_seg <- sample.seg[sample.seg$chrom==chrm, ]
      chrom.length <- max(chrom_seg$loc.end) - min(chrom_seg$loc.start)
      altered_any <- chrom_seg[chrom_seg$seg.mean > thresh | chrom_seg$seg.mean < (-thresh),]
      dist.altered_any <- sum(altered_any$loc.end - altered_any$loc.start)
      length.altered_any <- length.altered_any + dist.altered_any
      genome.length <- genome.length + chrom.length
    }
    
    # cat(paste("\nLength of genome covered:", genome.length/1000000,"Mb\n"))
    pga <- (length.altered_any/genome.length)*100
    num_segs = nrow(sample.seg)
    # genome.length.mb <- paste0((genome.length/1000000),"Mb")
    tcga.ID <- unique(as.character(sample.seg$tcga_ID))
    additional.desc <- unique(as.character(sample.seg$additional_desc))
    pga.list <- c(tcga.ID, additional.desc, sample, num_segs, genome.length)
    PGA <- rbind(PGA, pga.list)
  }
  colnames(PGA) <- PGA[1,]
  PGA <- data.frame(PGA[-1,])
  PGA$fraction_genome_altered <- as.numeric(as.character(PGA$num_segs))/as.numeric(as.character(PGA$genome_length))
  return (PGA)
}

# compute the percent genome altered in MPNST cases 

mpnst.seg <- "/Volumes/gsiprojects/external/zadehglab/MPNST/data/MPNST/cbioportal/cBioWrap_20190909/output/cbioportal_import_data/data_segments.txt"
mpnst.seg.data <- read.csv(mpnst.seg, sep = "\t", as.is = T)
# add additional annotations 
clinical.file <- "~/Documents/GSI/MPSNT/clinical/clinical_info.txt"
clinical.data <- read.csv(clinical.file, sep ="\t", as.is = T)
clinical.data$Sample_ID <- paste0("Tumor_", clinical.data$Tumor.Bank.No.)
head(clinical.data)

mpnst.seg.data <- merge(x = mpnst.seg.data, y = clinical.data, by.x = "ID", by.y = "Sample_ID")
mpnst.seg.data$tcga_ID <- paste0("GSI ",ifelse(mpnst.seg.data$Tumor.Type == "Neurofibroma", "NF", ifelse(mpnst.seg.data$Tumor.Type == "Premalignant_NF", "pre-NF", "MPNST")))
mpnst.seg.data$additional_desc <- mpnst.seg.data$Tumor.Type

head(mpnst.seg.data)


gsi.mpnst.pga <- computePGA(mpnst.seg.data)
head(gsi.mpnst.pga)

# compute percent genome altered in TCGA cases 
tcga.seg.csv <- "/Volumes/gsiprojects/external/zadehglab/MPNST/TCGA/tcga.pan_can.cbiodir.map.seg.csv"
tcga.seg.map <- read.csv(tcga.seg.csv, header = F)
head(tcga.seg.map)

TCGA.pga <- c("tcga_ID", "additional_desc", "Sample_ID", "num_segs", "genome_length")
for (v3 in tcga.seg.map$V3){
  seg.fl <- paste0("/Volumes/gsiprojects/external/zadehglab/MPNST/TCGA/", v3)
  tcga.ID <- unique(tcga.seg.map[tcga.seg.map$V3 == v3,]$V1)
  additional.desc <- unique(tcga.seg.map[tcga.seg.map$V3 == v3,]$V2)
  if (file.exists(seg.fl)){
    tcga.seg <- read.csv(seg.fl, as.is = T, sep = "\t")
    # tcga.seg$tcga_ID <- tcga.ID
    # tcga.seg$additional_desc <- additional.desc
    tcga.seg$tcga_ID <- rep(tcga.ID, dim(tcga.seg)[1])
    tcga.seg$additional_desc <- rep(additional.desc, dim(tcga.seg)[1])
    TCGA.pga <- rbind(TCGA.pga, computePGA(tcga.seg))
  }
}

# 
# TCGA.combined.seg <- data.frame(TCGA.combined.seg[-1,])
# head(TCGA.combined.seg)
# 
# TCGA.pga <- computePGA(TCGA.combined.seg)
TCGA.pga <- data.frame(TCGA.pga[-1,])
TCGA.pga <- TCGA.pga[TCGA.pga$tcga_ID %in% allowed.cancer.types, ]
head(TCGA.pga)

combined.TCGA.mpnst.pga <- data.frame(rbind(gsi.mpnst.pga, TCGA.pga))
head(combined.TCGA.mpnst.pga)

# rank them as per tcga_ID
# compute percentile for each sample per cancer type
newTCGA.percentile <- c("tcga_ID",
                        "additional_desc",
                        "Sample_ID",
                        "num_segs",
                        "genome_length",          
                        "fraction_genome_altered",
                        "percentile_rank",
                        "logPct",
                        "logPGA",
                        "median.PGA")
for (cancertype in unique(combined.TCGA.mpnst.pga$tcga_ID)){
  TCGA.pga.ctype <- combined.TCGA.mpnst.pga[combined.TCGA.mpnst.pga$tcga_ID == cancertype,]
  head(TCGA.pga.ctype)
  TCGA.pga.ctype$fraction_genome_altered <- as.numeric(as.character(TCGA.pga.ctype$fraction_genome_altered))
  TCGA.pga.ctype <- mutate(TCGA.pga.ctype, percentile_rank = ntile(TCGA.pga.ctype$fraction_genome_altered,100))
  TCGA.pga.ctype$logPct <- log10(as.numeric(TCGA.pga.ctype$percentile_rank))
  TCGA.pga.ctype$logPGA <- log10(as.numeric(as.character(TCGA.pga.ctype$fraction_genome_altered)))
  TCGA.pga.ctype <- TCGA.pga.ctype[!is.na(TCGA.pga.ctype$Sample_ID),]
  TCGA.pga.ctype <- TCGA.pga.ctype[! is.na(TCGA.pga.ctype$fraction_genome_altered),]
  # TCGA.tmb.ctype$Mutation_burden <- ifelse( is.na(TCGA.tmb.ctype$Mutation_burden), 0, TCGA.tmb.ctype$Mutation_burden)
  pga.median <- (median(TCGA.pga.ctype$fraction_genome_altered))
  TCGA.pga.ctype$median.PGA <- rep(pga.median, dim(TCGA.pga.ctype)[1])
  newTCGA.percentile <- rbind(newTCGA.percentile, TCGA.pga.ctype)
}
newTCGA.percentile <- data.frame(newTCGA.percentile[-1,])
newTCGA.percentile$percentile_rank <- as.numeric(as.character(newTCGA.percentile$percentile_rank))
newTCGA.percentile$fraction_genome_altered <- as.numeric(as.character(newTCGA.percentile$fraction_genome_altered))
newTCGA.percentile$logPGA <- as.numeric(as.character(newTCGA.percentile$logPGA))
newTCGA.percentile$median.PGA <- as.numeric(as.character(newTCGA.percentile$median.PGA))
newTCGA.percentile$tcga_ID <- as.character(newTCGA.percentile$tcga_ID)
# levels(newTCGA.percentile$median.TMB) <-  sort(unique(newTCGA.percentile$median.TMB))
newTCGA.percentile$tcga_ID <-
  factor(newTCGA.percentile$tcga_ID, levels = unique(newTCGA.percentile$tcga_ID[order(newTCGA.percentile$median.PGA)]))
newTCGA.percentile$log10median <- log10(as.numeric(as.character(newTCGA.percentile$median.PGA)))



TCGA.pgaplot <- ggplot(data.frame(newTCGA.percentile), aes (x = percentile_rank , y = fraction_genome_altered)) + geom_point(size = 0.3) +
  facet_grid(.~tcga_ID, scales = "free_x", space = "free", labeller = label_wrap_gen(width = 5, multi_line = TRUE)) + 
  xlab("") + 
  ylab("Fraction genome altered (TCGA)") +
  theme_classic() + theme(strip.text.y = element_text(angle = 90), axis.text.x = element_blank()) + scale_y_log10()
TCGA.pgaplot
TCGA.pgaplot <- TCGA.pgaplot + geom_point(data=newTCGA.percentile, aes(x = percentile_rank , y = median.PGA), 
                                          size = 0.5, shape = "-", color = "orange")

TCGA.pgaplot
TCGA.pgaplot <- TCGA.pgaplot + geom_point(data = newTCGA.percentile[grepl("GSI", newTCGA.percentile$tcga_ID),], 
                                          aes (x = percentile_rank , y = fraction_genome_altered, color = tcga_ID), size = 0.6)

TCGA.pgaplot

TCGA.pgaplot <- TCGA.pgaplot + geom_point(data=newTCGA.percentile[grepl("GSI",newTCGA.percentile$tcga_ID),], 
                                          aes(x = percentile_rank , y = median.PGA), size = 0.7, shape = "-", color = "black")
TCGA.pgaplot


ggsave("~/Documents/GSI/MPSNT/TMB/allPGATCGA.png", TCGA.pgaplot, width = 15, height = 5)



#######