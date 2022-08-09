# to obtain TMB for this project
library(dplyr)

mpnst.tmb <- read.csv("~/Documents/GSI/MPSNT/TMB/plots/tmb.txt", sep = "\t", as.is = T)
head(mpnst.tmb)
mpnst.tmb$Project <- "GSI-MPNST"
mpnst.tmb$Cancer.Type <- "mpnst"
mpnst.tmb$colorBy <- "GSI-MPNST"
mpnst.tmb <- mpnst.tmb[! grepl("Control", mpnst.tmb$Sample_ID),]


# annotate GSI-MPNST with project types (definitions from clinical annotations file)

clinical.file <- "~/Documents/GSI/MPSNT/clinical/clinical_info.txt"
clinical.data <- read.csv(clinical.file, sep ="\t", as.is = T)
clinical.data$Sample_ID <- paste0("Tumor_", clinical.data$Tumor.Bank.No.)
head(clinical.data)

mpnst.clin.tmb <- merge(x = mpnst.tmb, y = clinical.data, by = "Sample_ID")
mpnst.clin.tmb$colorBy <- paste0("GSI ",ifelse(mpnst.clin.tmb$Tumor.Type == "Neurofibroma", "NF", ifelse(mpnst.clin.tmb$Tumor.Type == "Premalignant_NF", "pre-NF", "MPNST")))


mpnst.tmb <- mpnst.clin.tmb

# obtain TMB for all TCGA projects 

tmb.for.all.tcga <- read.csv("/Volumes/gsiprojects/external/zadehglab/MPNST/TCGA/TCGA_exonic_tmb.tsv", sep = "\t", as.is = T)
head(tmb.for.all.tcga)

# plot TMB of MPNST vs all other TCGA projects 

combined.TMB <- tmb.for.all.tcga[,c("TCGA.Sample.Name", "TCGA.tmb", "CANCER.TYPE")]
names(combined.TMB) <- c("Sample_ID", "Mutation_burden", "Project")
combined.TMB$Cancer.Type <- stringr::str_split_fixed(combined.TMB$Project, "-", 2)[,1]
combined.TMB$colorBy <- combined.TMB$Cancer.Type
combined.TMB <- data.frame(rbind(combined.TMB, mpnst.tmb[,c("Sample_ID", "Mutation_burden", "Project", "Cancer.Type", "colorBy")]))

# restrict to cancer types 

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

combined.TMB <- combined.TMB[combined.TMB$Cancer.Type %in% allowed.cancer.types,]

newTCGA.percentile <- c("Sample_ID", "Mutation_burden", "Project", "Cancer.Type", "colorBy", "percentile_rank")

# compute percentile for each sample per cancer type
for (cancertype in unique(combined.TMB$colorBy)){
  TCGA.tmb.ctype <- combined.TMB[combined.TMB$colorBy == cancertype,]
  # TCGA.tmb.ctype$PERCENTILE <- ecdf(TCGA.tmb.ctype$TCGA.tmb)
  TCGA.tmb.ctype$Mutation_burden <- as.numeric(as.character(TCGA.tmb.ctype$Mutation_burden))
  TCGA.tmb.ctype <- mutate(TCGA.tmb.ctype, percentile_rank = ntile(TCGA.tmb.ctype$Mutation_burden,100))
  TCGA.tmb.ctype$logPct <- log10(as.numeric(TCGA.tmb.ctype$percentile_rank))
  TCGA.tmb.ctype$logTMB <- log10(as.numeric(as.character(TCGA.tmb.ctype$Mutation_burden)))
  TCGA.tmb.ctype <- TCGA.tmb.ctype[!is.na(TCGA.tmb.ctype$Sample_ID),]
  TCGA.tmb.ctype <- TCGA.tmb.ctype[! is.na(TCGA.tmb.ctype$Mutation_burden),]
  # TCGA.tmb.ctype$Mutation_burden <- ifelse( is.na(TCGA.tmb.ctype$Mutation_burden), 0, TCGA.tmb.ctype$Mutation_burden)
  tmb.median <- (median(TCGA.tmb.ctype$Mutation_burden))
  TCGA.tmb.ctype$median.TMB <- rep(tmb.median, dim(TCGA.tmb.ctype)[1])
  newTCGA.percentile <- rbind(newTCGA.percentile, TCGA.tmb.ctype)
}
newTCGA.percentile <- data.frame(newTCGA.percentile[-1,])
newTCGA.percentile$percentile_rank <- as.numeric(as.character(newTCGA.percentile$percentile_rank))
newTCGA.percentile$Mutation_burden <- as.numeric(as.character(newTCGA.percentile$Mutation_burden))
newTCGA.percentile$logTMB <- as.numeric(as.character(newTCGA.percentile$logTMB))
newTCGA.percentile$median.TMB <- as.numeric(as.character(newTCGA.percentile$median.TMB))
newTCGA.percentile$Cancer.Type <- as.character(newTCGA.percentile$Cancer.Type)
# levels(newTCGA.percentile$median.TMB) <-  sort(unique(newTCGA.percentile$median.TMB))
newTCGA.percentile$colorBy <-
  factor(newTCGA.percentile$colorBy, levels = unique(newTCGA.percentile$colorBy[order(newTCGA.percentile$median.TMB)]))
newTCGA.percentile$log10median <- log10(as.numeric(as.character(newTCGA.percentile$median.TMB)))


###
strwrap_strip_text = function(p, pad=0.05) { 
  # get facet font attributes
  th = theme_get()
  if (length(p$theme) > 0L)
    th = th + p$theme
  
  require("grid")
  grobs <- ggplotGrob(p)
  
  # wrap strip x text
  if ((class(p$facet)[1] == "grid" && !is.null(names(p$facet$cols))) ||
      class(p$facet)[1] == "wrap")
  {
    ps = calc_element("strip.text.x", th)[["size"]]
    family = calc_element("strip.text.x", th)[["family"]]
    face = calc_element("strip.text.x", th)[["face"]]
    
    if (class(p$facet)[1] == "wrap") {
      nm = names(p$facet$facets)
    } else {
      nm = names(p$facet$cols)
    }
    
    # get number of facet columns
    levs = levels(factor(p$data[[nm]]))
    npanels = length(levs)
    if (class(p$facet)[1] == "wrap") {
      cols = n2mfrow(npanels)[1]
    } else {
      cols = npanels
    }
    
    # get plot width
    sum = sum(sapply(grobs$width, function(x) convertWidth(x, "in")))
    panels_width = par("din")[1] - sum  # inches
    # determine strwrap width
    panel_width = panels_width / cols
    mx_ind = which.max(nchar(levs))
    char_width = strwidth(levs[mx_ind], units="inches", cex=ps / par("ps"), 
                          family=family, font=gpar(fontface=face)$font) / 
      nchar(levs[mx_ind])
    width = floor((panel_width - pad)/ char_width)  # characters
    
    # wrap facet text
    p$data[[nm]] = unlist(lapply(strwrap(p$data[[nm]], width=width, 
                                         simplify=FALSE), paste, collapse="\n"))
  }
  
  if (class(p$facet)[1] == "grid" && !is.null(names(p$facet$rows))) {  
    ps = calc_element("strip.text.y", th)[["size"]]
    family = calc_element("strip.text.y", th)[["family"]]
    face = calc_element("strip.text.y", th)[["face"]]
    
    nm = names(p$facet$rows)
    
    # get number of facet columns
    levs = levels(factor(p$data[[nm]]))
    rows = length(levs)
    
    # get plot height
    sum = sum(sapply(grobs$height, function(x) convertWidth(x, "in")))
    panels_height = par("din")[2] - sum  # inches
    # determine strwrap width
    panels_height = panels_height / rows
    mx_ind = which.max(nchar(levs))
    char_height = strwidth(levs[mx_ind], units="inches", cex=ps / par("ps"), 
                           family=family, font=gpar(fontface=face)$font) / 
      nchar(levs[mx_ind])
    width = floor((panels_height - pad)/ char_height)  # characters
    
    # wrap facet text
    p$data[[nm]] = unlist(lapply(strwrap(p$data[[nm]], width=width, 
                                         simplify=FALSE), paste, collapse="\n"))
  }
  
  invisible(p)
}
###

# plot ranks for each cancer type
TCGA.tmbplot <- ggplot(data.frame(newTCGA.percentile), aes (x = percentile_rank , y = Mutation_burden)) + geom_point(size = 0.3) +
  facet_grid(.~colorBy, scales = "free_x", space = "free", labeller = label_wrap_gen(width = 5, multi_line = TRUE)) + 
               xlab("") + 
               ylab("Mutation Burden (TCGA)") +
  theme_classic() + theme(strip.text.y = element_text(angle = 90), axis.text.x = element_blank()) + scale_y_log10()
TCGA.tmbplot
TCGA.tmbplot <- TCGA.tmbplot + geom_point(data=newTCGA.percentile, aes(x = percentile_rank , y = median.TMB), 
                                          size = 0.5, shape = "-", color = "orange")
TCGA.tmbplot <- TCGA.tmbplot + geom_point(data = newTCGA.percentile[newTCGA.percentile$Project == "GSI-MPNST",], 
                                      aes (x = percentile_rank , y = Mutation_burden, color = colorBy), size = 0.6)

TCGA.tmbplot <- TCGA.tmbplot + geom_point(data=newTCGA.percentile[newTCGA.percentile$Project == "GSI-MPNST",], 
                                          aes(x = percentile_rank , y = median.TMB), size = 0.6, shape = "-", color = "black")

TCGA.tmbplot


ggsave("~/Documents/GSI/MPSNT/TMB/allTMBTCGA.png", TCGA.tmbplot, width = 15, height = 5)






# for each Project ; rank TMB ; plot ; plot median ; 
library(ggplot2)

p <- ggplot(combined.TMB) + 
  geom_boxplot(aes(x = reorder(Project, Cancer.Type), y = Mutation_burden, fill = Cancer.Type)) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) + 
  facet_grid(.~colorBy, scales = "free", space = "free") +
  scale_y_log10()

p


p <- ggplot(combined.TMB) + 
  geom_boxplot(aes(x = colorBy, y = Mutation_burden, fill = Cancer.Type)) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) + 
  facet_grid(.~colorBy, scales = "free", space = "free") +
  scale_y_log10()

p


# save percentile rank of each MPNST sample with respect to all cancer types 

# 