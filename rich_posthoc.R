
postHocAnova <- function(pth = "Hallmark", outdir, groupName){
  dir.create(outdir, showWarnings = T, recursive = T)
  gsva_es <- read.csv(paste0("/Volumes/gsiprojects/external/svhaAU/RICH/ssGSEA/results/RICH2.genes_all_samples_FPKM.txt.",pth, ".ssGSEA.txt"),
                      sep = "\t")
  gsva_es$Geneset <- row.names(gsva_es)
  head(gsva_es)
  # annotations
  annotDF <- read.csv("/Volumes/gsiprojects/external/svhaAU/RICH/rich_master_annot.txt", sep = "\t")
  row.names(annotDF) <- gsub("-", ".", annotDF$internal.sample.name)
  annotDF$SAMPLE_ID <- row.names(annotDF)
  annotDF[,"condition"] <- gsub(" ", "_", annotDF[,groupName])
  annotDF[,"condition"] <- gsub("[+]", "_plus_", annotDF[,"condition"])
  annotDF
  # pth <- "Immune-Dahner"
  data.pth <- melt(gsva_es, id = "Geneset")
  colnames(data.pth) <- c("Geneset", "SAMPLE_ID", "Enrichment_Scores")
  # annotate conditions
  data.pth <- merge(data.pth, annotDF, by = "SAMPLE_ID")
  data.pth$conditions <- data.pth[, groupName]
  kw.sig <- c()
  pval <- c()
  aov.sig <- c()
  tuk <- c()
  for (gs in unique(data.pth$Geneset)){
    # compute P-value - Kruskal-Wallis Test
    kw <- kruskal.test(Enrichment_Scores ~ condition, data = data.pth[data.pth$Geneset == gs,])
    # res.aov <- aov(Enrichment_Scores ~ condition, data = data.pth[data.pth$Geneset == gs,])
    # P < 0.02 pathways (report this)
    pval[[gs]] <- kw$p.value
    # if (pval[[gs]] < 0.05){
    kw.sig[[gs]] <- kw
    aov.sig[[gs]] <- aov(Enrichment_Scores ~ condition, data = data.pth[data.pth$Geneset == gs,])
    tuk[[gs]] <- TukeyHSD(aov.sig[[gs]])
    # }
  }
  
  Total.genesets <- length(unique(data.pth$Geneset))
  print (paste0("Total number of genesets in ", pth, " = ", Total.genesets))
  
  # compute adjusted p-values 
  data.pth$pvalues <- pval[data.pth$Geneset]
  Significant.genesets <- length(unique(data.pth[data.pth$pvalues < 0.05,]$Geneset))
  print (paste0("Number of significant genesets (p < 0.05) according to non-parametric Kruskal-Walis test for ", pth, " = ", Significant.genesets))
  data.pth$padj <- p.adjust(data.pth$pvalues, method = "fdr")
  write.table(apply(data.pth,2,as.character), file =  paste0(outdir,"/", pth, "_Kruskal_Walis_pathways.txt"), 
              sep = "\t",
              quote = F,
              row.names = F)
  # subset geneset by significant 
  data.pth.sig <- data.pth[data.pth$padj < 0.05,]
  Significant.genesets.fdr <- length(unique(data.pth.sig$Geneset))
  print (paste0("Number of significant genesets post FDR correction (p.adjust < 0.05) according to non-parametric Kruskal-Walis test for ", pth, " = ", Significant.genesets.fdr))
  
  if (dim(data.pth.sig)[1] > 0){
    # write to file
    write.table(data.pth.sig, file =  paste0(outdir,"/", pth, "_Kruskal_Walis_Significant_pathways.txt"), 
                sep = "\t",
                quote = F,
                row.names = F)
    
    # plot this 
    p <- ggplot(data = data.pth.sig, aes(x = condition, y = Enrichment_Scores)) + 
      geom_boxplot(aes(color = condition)) +
      geom_jitter(size = 0.6) +
      theme_minimal() + 
      facet_wrap(~Geneset, scales = "free") + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 0.6)) + 
      ggtitle("Kruskal-Walis test significant (p.adjust < 0.05) genesets enriched per condtional group")
    p
    ggsave(file = paste0(outdir,"/", pth, "_Kruskal_Walis_Significant_pathways.png"), p,  width = 25, height = 20)
    
  }
  
  # F value cut-off
  all.genesets = unique(data.pth$Geneset)
  # matrix of F value and P-value from Anova
  anova.matrix <- matrix(NA, nrow=length(all.genesets), ncol=2)
  row.names(anova.matrix) <- all.genesets
  colnames(anova.matrix) <- c("F", "p.value")
  for (geneset in all.genesets){
    aov_cont <- aov.sig[[geneset]]
    anova.matrix[geneset, "F"] <- summary(aov_cont)[[1]]$F[1]
    anova.matrix[geneset, "p.value"] <- summary(aov_cont)[[1]]$P[1]
  }
  
  anova.matrix <- data.frame(anova.matrix)
  write.table(anova.matrix, file =  paste0(outdir,"/", pth, "_ANOVA_F_statistic.txt"), 
              sep = "\t",
              quote = F,
              row.names = T)
  
  # The F statistic is obtained from Anova for F([condition,Df], [Residuals,Df]) = fstatistic; p < Pr(>F)
  
  # matrix of post-hoc test
  all.genesets = unique(data.pth$Geneset)
  # matrix of F value and P-value from Anova
  posthoc.matrix <- matrix(NA, nrow=length(all.genesets), ncol=4)
  row.names(posthoc.matrix) <- all.genesets
  colnames(posthoc.matrix) <- unique(annotDF$condition)
  for (geneset in all.genesets){
    for (condition.compare in unique(annotDF$condition)){
      posthoc.matrix[geneset, condition.compare] <- tuk[[geneset]]$condition[,"p adj"][condition.compare]
    }
  }
  posthoc.matrix <- data.frame(posthoc.matrix)
  write.table(posthoc.matrix, file =  paste0(outdir,"/", pth, "_ANOVA_post_hoc.txt"), 
              sep = "\t",
              quote = F,
              row.names = T)
  
  # from the F-statistic anova matrix extract significant pathways/genesets
  anova.matrix.sig <- anova.matrix[anova.matrix$p.value < 0.05,]
  write.table(anova.matrix.sig, file =  paste0(outdir,"/", pth, "_significant_ANOVA_F_statistic.txt"), 
              sep = "\t",
              quote = F,
              row.names = T)
  
  Anova.Sig.Pathways <- dim(anova.matrix.sig)[1]
  print (paste0("Number of significant genesets (p.adjust < 0.05) according to one-way ANOVA test for ", pth, " = ", Anova.Sig.Pathways))
  
  
  # select the most relevant conditional groups from posthoc matrix
  posthoc.matrix.sig <- posthoc.matrix[row.names(anova.matrix.sig),]
  posthoc.matrix.sig$Genesets <- row.names(posthoc.matrix.sig)
  posthoc.matrix.sig <- melt(posthoc.matrix.sig)
  colnames(posthoc.matrix.sig) <- c("Genesets", "conditional_groups", "p.value")
  posthoc.matrix.sig$rank <- NA
  rank.order.anova.list <- c()
  rank.order.anova.top <- c()
  # rank condition groups based on increasing value of p-value for each pathway
  for (geneset in unique(posthoc.matrix.sig$Genesets)){
    ranked.matrix <- posthoc.matrix.sig[posthoc.matrix.sig$Genesets == geneset,]
    ranked.matrix$rank <- rank(ranked.matrix$p.value)
    posthoc.matrix.sig[posthoc.matrix.sig$Genesets == geneset,]$rank <- rank(ranked.matrix$p.value)
    # geneset.df <- posthoc.matrix.sig.rel[posthoc.matrix.sig.rel$Genesets == geneset,]
    rank.order.anova <- ranked.matrix$conditional_groups[order(ranked.matrix$rank)] 
    rank.order.anova.list[[geneset]] <- as.character(rank.order.anova)
    rank.order.anova.top[[geneset]] <- as.character(rank.order.anova[[1]])
  }
  
  rank.order.anova.list <- data.frame(t(data.frame(rank.order.anova.list)))
  write.table(rank.order.anova.list, file =  paste0(outdir,"/", pth, "_ordered_ANOVA_post_hoc_groups.txt"), 
              sep = "\t",
              quote = F,
              row.names = T)
  
  # plot the enrichment scores
  # plot this 
  posthoc.matrix.sig.anova <-data.pth[data.pth$Geneset %in% unique(posthoc.matrix.sig$Genesets),]
  anova.sig.enrichment <- ggplot(data = posthoc.matrix.sig.anova, aes(x = condition, y = Enrichment_Scores)) + 
    geom_boxplot(aes(color = condition)) +
    geom_jitter(size = 0.6) +
    theme_minimal() + 
    facet_wrap(~Geneset, scales = "free") + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 0.6))
  anova.sig.enrichment
  ggsave(file = paste0(outdir,"/", pth, "_ANOVA_Significant_pathways.png"), anova.sig.enrichment,  width = 25, height = 20)
  
  # plot p-values for each gene set
  anova.plot <- ggplot(data = posthoc.matrix.sig, aes(x = reorder(conditional_groups,rank), y = p.value)) + 
    geom_point(aes(color = conditional_groups)) + 
    geom_hline(yintercept = 0.05, color = "red", linetype = "dashed") +
    geom_line(group = 1) + 
    theme_bw() +
    scale_y_log10() +
    facet_wrap(~Genesets, scales = "free") +
    theme(axis.text.x = element_blank()) + 
    xlab("conditional ANOVA groups") +
    ggtitle("p-value - Tukey HSD (Honestly Significant Difference)")
  anova.plot
  ggsave(file = paste0(outdir,"/", pth, "_ANOVA_groups_pvalue_rank.png"), anova.plot,  width = 25, height = 20)
  
  
  # plot ranks for each gene set
  anova.rank.plot <- ggplot(data = posthoc.matrix.sig, aes(x = reorder(conditional_groups,rank), y = rank)) + 
    geom_point(aes(color = conditional_groups)) + 
    geom_line(group = 1) + 
    theme_bw() +
    facet_wrap(~Genesets, scales = "free") + 
    theme(axis.text.x = element_blank()) + 
    xlab("conditional ANOVA groups") + 
    ggtitle("Ranking conditional groups - Tukey HSD (Honestly Significant Difference) ")
  anova.rank.plot
  ggsave(file = paste0(outdir,"/", pth, "_ANOVA_groups_rank.png"), anova.rank.plot,  width = 25, height = 20)
  
  
  # retain significant condition ANOVA groups with p.value < 0.05
  # posthoc.matrix.sig.rel <- posthoc.matrix.sig[posthoc.matrix.sig$p.value < 0.05,]
  # write.table(posthoc.matrix.sig.rel, file =  paste0(outdir,"/", pth, "_significant_genesets_ANOVA_post_hoc_groups.txt"), 
  #             sep = "\t",
  #             quote = F,
  #             row.names = F)
  # # plot p-values for each gene set
  # anova.plot.sig <- ggplot(data = posthoc.matrix.sig.rel, aes(x = reorder(conditional_groups,rank), y = p.value)) + 
  #   geom_point(aes(color = conditional_groups)) + 
  #   geom_hline(yintercept = 0.05, color = "red", linetype = "dashed") +
  #   geom_line(group = 1) + 
  #   theme_bw() +
  #   scale_y_log10() +
  #   facet_wrap(~Genesets, scales = "free") +
  #   theme(axis.text.x = element_blank()) + 
  #   xlab("conditional ANOVA groups") +
  #   ggtitle("p-value - Tukey HSD (Honestly Significant Difference) for significantly variant conditional groups")
  # anova.plot.sig
  # ggsave(file = paste0(outdir,"/", pth, "_ANOVA_significantly_variant_groups_rank.png"), anova.plot.sig,  width = 25, height = 20)
  
}



#####


GMTMAP <- c("C7" = "/.mounts/labs/gsiprojects/external/svhaAU/RICH/ssGSEA/msigDB/c7.all.v7.4.symbols.gmt",
            "C6" = "/.mounts/labs/gsiprojects/external/svhaAU/RICH/ssGSEA/msigDB/c6.all.v7.4.symbols.gmt",
            "Hallmark" = "/.mounts/labs/gsiprojects/external/svhaAU/RICH/ssGSEA/msigDB/h.all.v7.4.symbols.gmt",
            "C2.CP" = "/.mounts/labs/gsiprojects/external/svhaAU/RICH/ssGSEA/msigDB/c2.cp.v7.4.symbols.gmt",
            "C2.CGP" = "/.mounts/labs/gsiprojects/external/svhaAU/RICH/ssGSEA/msigDB/c2.cgp.v7.4.symbols.gmt",
            "C5.GO" = "/.mounts/labs/gsiprojects/external/svhaAU/RICH/ssGSEA/msigDB/c5.go.v7.4.symbols.gmt",
            "C2.CP.Kegg" = "/.mounts/labs/gsiprojects/external/svhaAU/RICH/ssGSEA/msigDB/c2.cp.kegg.v7.4.symbols.gmt")
# test
pth = "C2.CGP"
outdir <- "~/Projects/GSI/Projects/RICH/ssGSEA/paired_group_analysis/res/"
dir.create(outdir,showWarnings = F)
groupName = "group1"
# end of test
postHocWithGoodPlots <- function(pth = "Hallmark", outdir, groupName = "group1"){
  dir.create(outdir, showWarnings = T, recursive = T)
  gsva_es <- read.csv(paste0("/Volumes/gsiprojects/external/svhaAU/RICH/ssGSEA/new_redo//RICH.genes_all_samples_FPKM.txt.",pth, ".ssGSEA.txt"),
                      sep = "\t")
  head(gsva_es)
  # annotations
  annotDF <- read.csv("/Volumes/gsiprojects/external/svhaAU/RICH/ssGSEA/new_redo/master_annot.txt", sep = "\t")
  row.names(annotDF) <- annotDF$effective.sample.name
  annotDF$SAMPLE_ID <- row.names(annotDF)
  annotDF[,"condition"] <- gsub(" ", "_", annotDF[,groupName])
  annotDF[,"condition"] <- gsub("[+]", "_plus_", annotDF[,"condition"])
  annotDF
  
  # plot better 
  library(ggpubr)
  library(reshape2)
  
  conditions <- unique(annotDF$condition)
  
  # write a global data frame
  cols.are <- c("Source","Geneset", 
                "kruskal.walis.p", "kruskal.walis.p.adj", 
                "anova.p", "anova.p.adj",
                paste0("parametric.pairwise.condition_",conditions[[1]], ".condition_", conditions[[2]],".p"), 
                paste0("parametric.pairwise.condition_",conditions[[1]], ".condition_", conditions[[3]],".p"), 
                paste0("parametric.pairwise.condition_",conditions[[1]], ".condition_", conditions[[4]],".p"), 
                paste0("parametric.pairwise.condition_",conditions[[2]], ".condition_", conditions[[3]],".p"), 
                paste0("parametric.pairwise.condition_",conditions[[2]], ".condition_", conditions[[4]],".p"), 
                paste0("parametric.pairwise.condition_",conditions[[3]], ".condition_", conditions[[4]],".p"), 
                paste0("parametric.pairwise.condition_",conditions[[1]], ".condition_", conditions[[2]],".padj"), 
                paste0("parametric.pairwise.condition_",conditions[[1]], ".condition_", conditions[[3]],".padj"), 
                paste0("parametric.pairwise.condition_",conditions[[1]], ".condition_", conditions[[4]],".padj"), 
                paste0("parametric.pairwise.condition_",conditions[[2]], ".condition_", conditions[[3]],".padj"), 
                paste0("parametric.pairwise.condition_",conditions[[2]], ".condition_", conditions[[4]],".padj"), 
                paste0("parametric.pairwise.condition_",conditions[[3]], ".condition_", conditions[[4]],".padj"),
                paste0("non-parametric.pairwise.condition_",conditions[[1]], ".condition_", conditions[[2]],".p"), 
                paste0("non-parametric.pairwise.condition_",conditions[[1]], ".condition_", conditions[[3]],".p"), 
                paste0("non-parametric.pairwise.condition_",conditions[[1]], ".condition_", conditions[[4]],".p"), 
                paste0("non-parametric.pairwise.condition_",conditions[[2]], ".condition_", conditions[[3]],".p"), 
                paste0("non-parametric.pairwise.condition_",conditions[[2]], ".condition_", conditions[[4]],".p"), 
                paste0("non-parametric.pairwise.condition_",conditions[[3]], ".condition_", conditions[[4]],".p"), 
                paste0("non-parametric.pairwise.condition_",conditions[[1]], ".condition_", conditions[[2]],".padj"), 
                paste0("non-parametric.pairwise.condition_",conditions[[1]], ".condition_", conditions[[3]],".padj"), 
                paste0("non-parametric.pairwise.condition_",conditions[[1]], ".condition_", conditions[[4]],".padj"), 
                paste0("non-parametric.pairwise.condition_",conditions[[2]], ".condition_", conditions[[3]],".padj"), 
                paste0("non-parametric.pairwise.condition_",conditions[[2]], ".condition_", conditions[[4]],".padj"), 
                paste0("non-parametric.pairwise.condition_",conditions[[3]], ".condition_", conditions[[4]],".padj"))
  
  all.path.deets <- matrix(NA, nrow = dim(gsva_es)[1], ncol =  length(cols.are))
  row.names(all.path.deets) <- row.names(gsva_es)
  colnames(all.path.deets) <- cols.are
  
  # define a pretty data frame
  gsva_es$Geneset <- row.names(gsva_es)
  good.df <- melt(gsva_es, id = "Geneset")
  colnames(good.df) <- c("Geneset", "SAMPLE_ID", "Enrichment_Score")
  good.df <- merge(good.df, annotDF[,c( "SAMPLE_ID", "condition")], by = "SAMPLE_ID")
  good.df$condition <- paste0("condition_",good.df$condition)
  
  # parametric test variables
  global.test.anova <- c()
  global.test.sig.parametric <- c()
  pair.wise.test.parametric <- c()
  post.hoc.plot.parametric <-c()
  
  # non-parametric test variables
  global.test.kruskal.walis <- c()
  global.test.sig.non.parametric <- c()
  pair.wise.test.non.parametric <- c()
  post.hoc.plot.non.parametric <-c()
  
  for (gs in unique(good.df$Geneset)){
    
    
    all.path.deets[gs, "Source"] <- pth
    all.path.deets[gs, "Geneset"] <- gs
    # ANOVA (parametric)
    global.test.anova[[gs]] <- compare_means(Enrichment_Score ~ condition,  data = good.df[good.df$Geneset == gs,] , method = "anova", p.adjust.method = "fdr")
    all.path.deets[gs, c("anova.p", "anova.p.adj")] <-  as.numeric(global.test.anova[[gs]][,2:3])
    # Kruskal-Walis (parametric)
    global.test.kruskal.walis[[gs]] <- compare_means(Enrichment_Score ~ condition,  data = good.df[good.df$Geneset == gs,] , method = "kruskal.test", p.adjust.method = "fdr")
    all.path.deets[gs, c("kruskal.walis.p", "kruskal.walis.p.adj")] <-  as.numeric(global.test.kruskal.walis[[gs]][,2:3])
    # non-parametric
    if (global.test.kruskal.walis[[gs]]$p.adj < 0.05){
      
      global.test.sig.non.parametric[[gs]] <- global.test.kruskal.walis[[gs]]
      # pair wise
      pair.wise.test.non.parametric[[gs]] <- compare_means(Enrichment_Score ~ condition,  data = good.df[good.df$Geneset == gs,], p.adjust.method = "fdr")
      
      # Pick up the most significant groups
      my_comparisons <- NA
      rel.groups <- data.frame(pair.wise.test.non.parametric[[gs]])[data.frame(pair.wise.test.non.parametric[[gs]])$p.adj < 0.05,]
      if (! is.null(rel.groups)){
        print (gs)
        groups.df <- rel.groups[,c("group1", "group2")]
        if (dim(groups.df)[1] > 0){
          my_comparisons <- list()
          for (g in 1:dim(groups.df)[1]){
            my_comparisons[[g]] <- as.character(groups.df[g,])
            npm <- paste0("non.parametric.pairwise.", as.character(groups.df[g,"group1"]), ".", as.character(groups.df[g,"group2"]))
            all.path.deets[gs, c(paste0(npm, ".p"), paste0(npm, ".padj"))] <-as.numeric(rel.groups[g,c("p", "p.adj")])
          }
        }
      }
      
      if (!is.na(my_comparisons)){
        # my_comparisons <- list(group1, group2, group3, group4, group5, group6)
        post.hoc.plot.non.parametric[[gs]] <- ggboxplot(good.df[good.df$Geneset == gs,], x = "condition", y = "Enrichment_Score",
                                                        color = "condition", palette = "jco") + 
          stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
          stat_compare_means() # Add global p-value 
        post.hoc.plot.non.parametric[[gs]] <- ggpar(post.hoc.plot.non.parametric[[gs]], main = paste(pth, gs, sep = " : "))
      }
    }
    
    # parametric
    if (global.test.anova[[gs]]$p.adj < 0.05){
      
      global.test.sig.parametric[[gs]] <- global.test.anova[[gs]]
      # pair wise
      pair.wise.test.parametric[[gs]] <- compare_means(Enrichment_Score ~ condition,  data = good.df[good.df$Geneset == gs,], method= "t.test", p.adjust.method = "fdr")
      
      # Pick up the most significant groups
      my_comparisons <- NA
      rel.groups <- data.frame(pair.wise.test.parametric[[gs]])[data.frame(pair.wise.test.parametric[[gs]])$p.adj < 0.05,]
      if (dim(rel.groups)[1] > 0){
        print (gs)
        groups.df <- rel.groups[,c("group1", "group2")]
        # groups.df$gs <- gs
        if (dim(groups.df)[1] > 0 | !is.na(groups.df$group1)){
          my_comparisons <- list()
          for (g in 1:dim(groups.df)[1]){
            my_comparisons[[g]] <- as.character(groups.df[g,])
            npm <- paste0("parametric.pairwise.", as.character(groups.df[g,"group1"]), ".", as.character(groups.df[g,"group2"]))
            if (!unique(c(paste0(npm, ".p"), paste0(npm, ".padj")) %in% colnames(all.path.deets))){
              npm <- paste0("parametric.pairwise.", as.character(groups.df[g,"group2"]), ".", as.character(groups.df[g,"group1"]))
            }
            all.path.deets[gs, c(paste0(npm, ".p"), paste0(npm, ".padj"))] <- as.numeric(rel.groups[g,c("p", "p.adj")])
          }
        } 
      }
      if (!is.na(my_comparisons[[1]])){
        # my_comparisons <- list(group1, group2, group3, group4, group5, group6)
        post.hoc.plot.parametric[[gs]] <- ggboxplot(good.df[good.df$Geneset == gs,], x = "condition", y = "Enrichment_Score",
                                                    color = "condition", palette = "jco") + 
          stat_compare_means(comparisons = my_comparisons, method = "t.test") + # Add pairwise comparisons p-value
          stat_compare_means(method = "anova") # Add global p-value # +
        # theme(axis.text.x = element_text(angle = 90, hjust = 1))
        post.hoc.plot.parametric[[gs]] <- ggpar(post.hoc.plot.parametric[[gs]], main = paste(pth, gs, sep = " : "))
      }
    }
  }
  # for all the plots; arrange them as a glob and tile together
  Total.number.of.genesets <- length(unique(good.df$Geneset))
  Total.number.of.significant.genesets.non.parametric <- length(names(global.test.sig.non.parametric))
  Total.number.of.significant.genesets.survived.post.hoc.non.parametric <- length(names(post.hoc.plot.non.parametric))
  
  Total.number.of.significant.genesets.parametric <- length(names(global.test.sig.parametric))
  Total.number.of.significant.genesets.survived.post.hoc.parametric <- length(names(post.hoc.plot.parametric))
  
  
  Total.number.of.significant.genesets.common.to.both.tests <- length(intersect(names(post.hoc.plot.non.parametric), names(post.hoc.plot.parametric)))
  
  print.this.line <- data.frame(c("For genesets ", pth), 
                                c("Total number of genesets = ", Total.number.of.genesets), 
                                c("Parametric test = ", "one-way Anova; t-test for pairwise"),
                                c("Total number of significant genesets = ", Total.number.of.significant.genesets.parametric),
                                c("Total number of significant genesets survived post-hoc =", Total.number.of.significant.genesets.survived.post.hoc.parametric),
                                c("Non-parametric test = ", "Kruskal Wallis; Wilcoxon for pairwise"),
                                c("Total number of significant genesets = ", Total.number.of.significant.genesets.non.parametric),
                                c("Total number of significant genesets survived post-hoc =", Total.number.of.significant.genesets.survived.post.hoc.non.parametric),
                                c("Total number of significant genesets common to both parametric and non-parametric tests =", Total.number.of.significant.genesets.common.to.both.tests))
  write.table(t(print.this.line), file =  paste0(outdir,"/", pth, "_summary_text.txt"), 
              sep = "\t",
              quote = F,
              row.names = F,
              col.names = F)
  
  # club together all plots
  pdf(paste0(outdir,"/", pth, "_non_parametric_summary_good_plots.pdf"))
  for (i in 1:Total.number.of.significant.genesets.survived.post.hoc.non.parametric) {
    print(post.hoc.plot.non.parametric[[i]])
  }
  dev.off()
  
  pdf(paste0(outdir,"/", pth, "_parametric_summary_good_plots.pdf"), width = 15, height = 10)
  for (i in 1:Total.number.of.significant.genesets.survived.post.hoc.parametric) {
    print(post.hoc.plot.parametric[[i]])
  }
  dev.off()
  
  
  # write the data frame to file
  all.path.deets <- data.frame(all.path.deets)
  write.table(all.path.deets, file =  paste0(outdir,"/", pth, "_statistical_summary.txt"), 
              sep = "\t",
              quote = F,
              row.names = F,
              col.names = T)
  
  
  int.path.deets <- all.path.deets[all.path.deets$Geneset %in% intersect(names(post.hoc.plot.non.parametric), names(post.hoc.plot.parametric)),]
  write.table(int.path.deets, file =  paste0(outdir,"/", pth, "_interesting_pathways_statistical_summary.txt"), 
              sep = "\t",
              quote = F,
              row.names = F,
              col.names = T)
}


#####



expr_data <- "/.mounts/labs/gsiprojects/external/svhaAU/RICH/rsem/RICH.genes_all_samples_FPKM.txt"
outdir <- "~/Projects/GSI/Projects/RICH/ssGSEA/paired_group_analysis/res//"
# outdir <- "/.mounts/labs/gsiprojects/external/svhaAU/RICH/ssGSEA/results"
dir.create(outdir, showWarnings = F, recursive = T, mode = "0777")
# ensFile <- "/u/prath/cBioWrap/files/ensemble_conversion.txt"
# 
# 
# 
# # run gsea for all pathways
# gsva_es <- run.ssGSEA(expr_data, ensFile, outdir)
outdir <- "~/Projects/GSI/Projects/RICH/ssGSEA/paired_group_analysis/res/"
# dir.create(outdir,showWarnings = F)
# groupName = "group1"

postHocWithGoodPlots(pth = "Hallmark", outdir, groupName = "group1")
postHocWithGoodPlots(pth = "C6", outdir, groupName = "group1")
postHocWithGoodPlots(pth = "C7", outdir, groupName = "group1")
postHocWithGoodPlots(pth = "C2.CP", outdir, groupName = "group1")
postHocWithGoodPlots(pth = "C2.CGP", outdir, groupName = "group1")
postHocWithGoodPlots(pth = "C5.GO", outdir, groupName = "group1")
postHocWithGoodPlots(pth = "C2.CP.Kegg", outdir, groupName = "group1")


# postHocAnova(pth = "C6", outdir = "~/Projects/GSI/Projects/RICH/ssGSEA/workDir2/group1/C6/" , groupName = "group1")


# postHocAnova(pth = "C7", outdir = "~/Projects/GSI/Projects/RICH/ssGSEA/workDir2/group1/C7/" , groupName = "group1")
# postHocAnova(pth = "C7", outdir = "~/Projects/GSI/Projects/RICH/ssGSEA/workDir2/group2/C7/" , groupName = "group2")
# postHocAnova(pth = "C7", outdir = "~/Projects/GSI/Projects/RICH/ssGSEA/workDir2/group3/C7/" , groupName = "group3")
# 
# # postHocWithGoodPlots(pth = "Hallmark", outdir = "~/Projects/GSI/Projects/RICH/ssGSEA/workDir/group1/HallMark/" , groupName = "group1")
# # postHocWithGoodPlots(pth = "C6", outdir = "~/Projects/GSI/Projects/RICH/ssGSEA/workDir/group1/C6/" , groupName = "group1")
# postHocWithGoodPlots(pth = "C7", outdir = "~/Projects/GSI/Projects/RICH/ssGSEA/workDir/group1/C7/" , groupName = "group1")
# 
# postHocWithGoodPlots(pth = "Hallmark", outdir = "~/Projects/GSI/Projects/RICH/ssGSEA/workDir/group2/HallMark/" , groupName = "group2")
# postHocWithGoodPlots(pth = "C6", outdir = "~/Projects/GSI/Projects/RICH/ssGSEA/workDir/group2/C6/" , groupName = "group2")
# postHocWithGoodPlots(pth = "C7", outdir = "~/Projects/GSI/Projects/RICH/ssGSEA/workDir/group2/C7/" , groupName = "group2")
# 
# 
# postHocWithGoodPlots(pth = "Hallmark", outdir = "~/Projects/GSI/Projects/RICH/ssGSEA/workDir/group3/HallMark/" , groupName = "group3")
# postHocWithGoodPlots(pth = "C6", outdir = "~/Projects/GSI/Projects/RICH/ssGSEA/workDir/group3/C6/" , groupName = "group3")
# postHocWithGoodPlots(pth = "C7", outdir = "~/Projects/GSI/Projects/RICH/ssGSEA/workDir/group3/C7/" , groupName = "group3")


# plot for specific pathways of interest
# my_comparisons <- list(group1, group2, group3, group4, group5, group6)
plot_selective <- function(pth, gene.sets.file, outdir){
  # pth <- "C5.GO"
  gene.sets <- readLines(gene.sets.file)
  # get good df
  dir.create(outdir, showWarnings = F, recursive = T)
  gsva_es <- read.csv(paste0("/Volumes/gsiprojects/external/svhaAU/RICH/ssGSEA/new_redo//RICH.genes_all_samples_FPKM.txt.",pth, ".ssGSEA.txt"),
                      sep = "\t")
  head(gsva_es)
  # annotations
  annotDF <- read.csv("/Volumes/gsiprojects/external/svhaAU/RICH/ssGSEA/new_redo/master_annot.txt", sep = "\t")
  row.names(annotDF) <- annotDF$effective.sample.name
  annotDF$SAMPLE_ID <- row.names(annotDF)
  annotDF[,"condition"] <- gsub(" ", "_", annotDF[,groupName])
  annotDF[,"condition"] <- gsub("[+]", "_plus_", annotDF[,"condition"])
  annotDF
  
  gsva_es$Geneset <- row.names(gsva_es)
  good.df <- melt(gsva_es, id = "Geneset")
  colnames(good.df) <- c("Geneset", "SAMPLE_ID", "Enrichment_Score")
  good.df <- merge(good.df, annotDF[,c( "SAMPLE_ID", "condition")], by = "SAMPLE_ID")
  good.df$condition <- paste0("condition_",good.df$condition)
  
  
  pdf(paste0(outdir,"/", pth, "_selective_parametric_summary_good_plots.pdf"), width = 15, height = 10)
  for (gs in gene.sets){
    good_plt<- ggboxplot(good.df[good.df$Geneset == gs,], x = "condition", y = "Enrichment_Score",
                         color = "condition", palette = "jco") + 
      stat_compare_means(comparisons = my_comparisons, method = "t.test") + # Add pairwise comparisons p-value
      stat_compare_means(method = "anova") # Add global p-value # +
    # theme(axis.text.x = element_text(angle = 90, hjust = 1))
    good_plt <- ggpar(good_plt, main = paste(pth, gs, sep = " : "))
    
    print (good_plt)
  }
  dev.off()
}

# pth <- "C5.GO"
# gene.sets <- readLines("~/Projects/GSI/Projects/RICH/ssGSEA/paired_group_analysis/GO_selection.txt")
plot_selective(pth = "C2.CGP", gene.sets.file = "~/Projects/GSI/Projects/RICH/ssGSEA/paired_group_analysis/Kw_CGP_selection.txt", outdir = "~/Projects/GSI/Projects/RICH/ssGSEA/paired_group_analysis/selected_KW/")

plot_selective(pth = "C2.CGP", gene.sets.file = "~/Projects/GSI/Projects/RICH/ssGSEA/paired_group_analysis/CGP_selection.txt", outdir = "~/Projects/GSI/Projects/RICH/ssGSEA/paired_group_analysis/selected_Anova/")
