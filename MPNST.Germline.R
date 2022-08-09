library(ggplot2)

BASEDIR <- "~/Documents/GSI/MPSNT/Germline/"

# script to analyse germline calls 
variants.file <- "/Volumes/gsiprojects/external/zadehglab/MPNST/data/MPNST/EXOME/Seqware_GATKHaplotypeCaller/JointCalling/Suganth_genelist_MPNST_jointGT.annotated.tsv"
variants.data <- read.csv(variants.file, sep = "\t", as.is = T)
head(variants.data)
# read TSV file 
dim(variants.data)

variants.data$HGVS_p <- paste(variants.data$SYMBOL, stringr::str_split_fixed(variants.data$HGVSp, ":",2)[,2], sep = ":")
variants.data$HGVS_c <- paste(variants.data$SYMBOL, stringr::str_split_fixed(variants.data$HGVSc, ":",2)[,2], sep = ":")
head(variants.data)

names(variants.data)
sample.names <- names(variants.data)[52:length(names(variants.data))]

# clinvar

# pathogenic hits

clinvar.entry <- variants.data[variants.data$CLIN_SIG != "",] 

# get statistics about #variants with different levels of clinvar annotations

clinvar.annot <- unique(clinvar.entry$CLIN_SIG)

clinvar.stat <- c("CLNSIG", "#variants")
for (annot in clinvar.annot){
  df <- clinvar.entry[clinvar.entry$CLIN_SIG == annot,]
  clinvar.stat <- rbind(clinvar.stat, c(annot, dim(df)[1]))
}
colnames(clinvar.stat) <- clinvar.stat[1,]
clinvar.stat <- data.frame(clinvar.stat[-1,])
clinvar.stat
clinvar.stat$X.variants <- as.numeric(as.character(clinvar.stat$X.variants))
# plot this 
p <- ggplot(clinvar.stat) + geom_bar(aes(x = CLNSIG, y = X.variants), fill = "red", stat = "identity") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Clinvar ClinSig") +
  ylab("# variants")
p

ggsave(paste0(BASEDIR, "/clinvar_stats.png"), p)


# now prioritize all variants

# get all rare variants
cut.off = 0.02

rare.variants <- variants.data[variants.data$ExAC_AF < cut.off |
                                 variants.data$gnomAD_exomes_AF < cut.off |
                                 variants.data$gnomAD_genomes_AF < cut.off |
                                 variants.data$X1000Gp3_AF < cut.off |
                                 is.na(variants.data$X1000Gp3_AF) |
                                 is.na(variants.data$gnomAD_exomes_AF) |
                                 is.na(variants.data$gnomAD_genomes_AF) |
                                 is.na(variants.data$ExAC_AF), ]

no.of.rare.vars <- length(unique(rare.variants$HGVS_p))


rare.variants.non.syn <- rare.variants[rare.variants$Consequence != "synonymous_variant",]
no.of.rare.non.syn.vars <- length(unique(rare.variants.non.syn$HGVS_p))

# clinsig ones 

clin.variants <- rare.variants.non.syn[grepl("pathogenic", rare.variants.non.syn$CLIN_SIG, ignore.case = T),]

no.of.rare.clinsig.vars <- length(unique(clin.variants$HGVS_p))

vus.variants <- rare.variants.non.syn[grepl("uncertain_signficance", rare.variants.non.syn$CLIN_SIG, ignore.case = T),]

no.of.rare.vus.vars <- length(unique(vus.variants$HGVS_p))

# filter out all pathogenic vars 
rare.variants.for.prior <- rare.variants.non.syn[!grepl("benign", rare.variants.non.syn$CLIN_SIG, ignore.case = T),]

non.benign.vars <- length(unique(rare.variants.for.prior$HGVS_p))

# prioritizing criteria
# polyphen, sift, CADD

# for each sample get GT 

write.table(rare.variants.for.prior, paste0(BASEDIR, "/germline.for.prior.txt") , sep = "\t", quote = F, row.names = F)

# from the clinvar entry now get the sample names, clinvar annot, chrom, start, end, ref, alt, consequence,HGVSp, HGVSc  

df <- data.frame(t(c("Total_variants" = length(unique(variants.data$HGVS_p)),
         "Rare_variants" = no.of.rare.vars,
  "Non-synonymous_variants" = no.of.rare.non.syn.vars,
  "Non CLINVAR bning variants" = non.benign.vars,
  "Pathogenic_CLINVAR" = no.of.rare.clinsig.vars)))
df
write.table(df, paste0(BASEDIR, "/funnel.txt") , sep = "\t", quote = F, row.names = F)

# priotize all 1/1 variants

# get all samples with GT == 1/1
all.names <- names(rare.variants.for.prior)[!names(rare.variants.for.prior) %in% sample.names]
hom.samples <- c(names(rare.variants.for.prior)[!names(rare.variants.for.prior) %in% sample.names], "sample.name")
names(hom.samples) <- hom.samples
ch.samples <- hom.samples
check.for.happloinsuff <- hom.samples

for (gene in unique(rare.variants.for.prior$SYMBOL)){
  gene.df <- rare.variants.for.prior[rare.variants.for.prior$SYMBOL == gene,]
  # count total number of variants 
  gene.vars <- length(unique(gene.df$HGVS_p))
  # if count == 1; then count number of samples with GT 1/1
  for (sample in sample.names){
    samples.gt <- stringr::str_split_fixed(gene.df[, sample], ":", 4)[,1]
    if (length(unique(samples.gt)) == 1 &
        gene.vars == 1 &
        unique(samples.gt) == "1/1") {
      hif <- cbind(gene.df[,all.names], rep(sample, dim(gene.df)[1]))
      names(hif)[length(names(hif))] <- "sample"
      hom.samples <- rbind(hom.samples, gene.df[,c(all.names, hif)])
    }
    else if (length(samples.gt) > 1 &
             unique(samples.gt) == "0/1" &
             gene.vars > 1){
      hif <- cbind(gene.df[,all.names], rep(sample, dim(gene.df)[1]))
      names(hif)[length(names(hif))] <- "sample"
      ch.samples <- rbind(ch.samples, hif)
    } else if (length(unique(samples.gt)) == 1 &
               gene.vars == 1 &
               unique(samples.gt) == "0/1"){
      hif <- cbind(gene.df[,all.names], rep(sample, dim(gene.df)[1]))
      names(hif)[length(names(hif))] <- "sample"
      check.for.happloinsuff <- rbind(check.for.happloinsuff, hif)
    }
  }
}


ch.samples <- data.frame(ch.samples[-1,])
hom.samples <- data.frame(hom.samples[-1,])
check.for.happloinsuff <- data.frame(check.for.happloinsuff[-1,])

# write these to files

write.table(ch.samples,paste0(BASEDIR, "/compound_hets.txt") , sep = "\t", quote = F, row.names = F)
write.table(check.for.happloinsuff,paste0(BASEDIR, "/possibly_happloinsufficient.txt") , sep = "\t", quote = F, row.names = F)
write.table(hom.samples,paste0(BASEDIR, "/hom_recessive.txt") , sep = "\t", quote = F, row.names = F)





# 