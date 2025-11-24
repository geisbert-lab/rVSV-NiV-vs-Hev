#!/usr/bin/env Rscript

## setup -----------------------------------------------------------------------
rm(list=ls(all.names=TRUE))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(edgeR))
theme_set(ggpubr::theme_pubr() + theme(legend.position="right"))
options(stringsAsFactors=FALSE)

# helper variables
cols.dpi <- c("Baseline"="#ffffcc", 
              "4 DPI"="#a1dab4", 
              "7 DPI"="#41b6c4", 
              "10 DPI"="#225ea8",
              "Terminal"="black")
cols.reg <- c(Up="#d7191c",
              None="lightgrey",
              Down="#2c7bb6")
cols.treat <- c(Control="black", Vaccinated="#225ea8")
shapes.treat <- c(Control=21, Vaccinated=22)

# helper functions
normalize.counts <- function(counts, threshold=20, minimum.samples=3) {
  # remove probes with <minimum.samples above threshold 
  abovethld <- (counts > threshold)
  sampabove <- rowSums(abovethld)
  keepgenes <- (sampabove >= minimum.samples)
  counts <- counts[keepgenes, ]
  
  # format full matrix to get CPM and return
  counts %>% 
    DGEList() %>% 
    calcNormFactors() %>% 
    cpm()
}
run.pca <- function(logcpm.matrix, meta.matrix) {
  # run PCA
  pca <- logcpm.matrix %>%
         t() %>%
         prcomp()
  
  # extract % explained
  pcs <- summary(pca)$importance["Proportion of Variance", 1:2]
  pcs <- round(100*pcs)
  pcs <- paste0(names(pcs), " (", pcs, "%)")
  
  # format PCA matrix for downstream plotting
  pca <- pca$x %>%
         as.data.frame() %>%
         rownames_to_column("ID") %>%
         select(ID, PC1, PC2) %>%
         left_join(meta.matrix, by="ID")
  
  # return PCA object
  list(PCA=pca,
       PCs=pcs)
}
diffexpr <- function(counts, model, contrast, psig=0.05, fsig=1) {
  # calculate results
  de <- counts %>%
        voom(model) %>% # voom normalization
        lmFit(model) %>% # fit model
        contrasts.fit(contrast) %>% # incorporate contrast
        eBayes() %>% # perform test
        topTable(num=Inf) %>% # get all genes
        # format columns
        rownames_to_column("Gene") %>%
        dplyr::rename(lfc=logFC, padj=adj.P.Val) %>%
        select(Gene, lfc, padj)
  
  # add significance and annotate regulation
  de$Significant <- (de$padj < psig & abs(de$lfc) > fsig)
  de$Regulation <- "None"
  de[de$Significant & de$lfc > 0, "Regulation"] <- "Up"
  de[de$Significant & de$lfc < 0, "Regulation"] <- "Down"
  
  # return annotated DE matrix
  return(de)
}

## load data -------------------------------------------------------------------
# metadata: format DE variables as factors
meta <- read.csv("samplesheet.csv") %>%
        mutate(Treatment=factor(Treatment, levels=names(cols.treat)),
               Timepoint=factor(Timepoint, levels=names(cols.dpi)))
rownames(meta) <- meta$ID

# thresholded counts: remove + and - controls
cmat <- read.csv("counts-thresholded.csv",
                 row.names=1)
cmat <- cmat[!str_detect(rownames(cmat), "^POS_|NEG_"), ]

# align rows and columns
x <- intersect(rownames(meta), colnames(cmat))
meta <- meta[x, ]
cmat <- cmat[, x]
rm(x)

# save non-log CPM for CIBERSORT
cmat %>%
  normalize.counts() %>%
  cpm() %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  write.table("analysis/cibersortx-input.tsv", sep="\t", row.names=FALSE)

# sample plot
meta %>%
  ggplot(aes(DPI, NHP)) +
  geom_line(aes(group=NHP)) +
  geom_point(aes(shape=Treatment, fill=Timepoint), size=3) +
  scale_fill_manual(values=cols.dpi) +
  scale_shape_manual(values=shapes.treat) +
  scale_x_continuous(breaks=unique(meta$DPI)) +
  labs(x="Days postinfection",
       y=NULL,
       title="Sample schema") +
  guides(shape=guide_legend(override.aes=list(fill="black")),
         fill=guide_legend(override.aes=list(pch=21)))
ggsave("analysis/sampling.png",
       units="in", width=4, height=3)

## naive clustering ------------------------------------------------------------
# PCA of all samples
pca <- cmat %>%
       normalize.counts() %>%
       cpm(log=TRUE) %>%
       run.pca(meta)
pca$PCA %>%
  ggplot(aes(PC1, PC2)) +
  geom_hline(yintercept=0, linetype=3) +
  geom_vline(xintercept=0, linetype=3) +
  geom_point(aes(fill=Timepoint, shape=Treatment), size=3) +
  scale_shape_manual(values=shapes.treat) +
  scale_fill_manual(values=cols.dpi) +
  guides(shape=guide_legend(override.aes=list(fill="black")),
         fill=guide_legend(override.aes=list(pch=21))) +
  labs(x=pca$PCs[1],
       y=pca$PCs[2],
       title="PCA: all samples")
ggsave("analysis/pca-all.png",
       units="in", width=4, height=3)
rm(pca)

# UMAP
cmat %>%
  normalize.counts() %>%
  cpm(log=TRUE) %>%
  t() %>%
  umap::umap() %>%
  magrittr::extract2("layout") %>%
  as.data.frame() %>%
  rownames_to_column("ID") %>%
  rename(UMAP1=V1, UMAP2=V2) %>%
  left_join(meta, by="ID") %>%
  ggplot(aes(UMAP1, UMAP2)) +
  geom_point(aes(fill=Timepoint, shape=Treatment), size=3) +
  scale_shape_manual(values=shapes.treat) +
  scale_fill_manual(values=cols.dpi) +
  guides(shape=guide_legend(override.aes=list(fill="black")),
         fill=guide_legend(override.aes=list(pch=21))) +
  labs(title="UMAP: all samples")
ggsave("analysis/umap-all.png",
       units="in", width=4, height=3)

## controls over time ----------------------------------------------------------
# subset data and format times
m <- meta %>%
     filter(Treatment=="Control") %>%
     mutate(Timepoint=droplevels(Timepoint))
c <- cmat[, rownames(m)] %>%
     normalize.counts()

# PCA
pca <- c %>%
       cpm(log=TRUE) %>%
       run.pca(m)
pca$PCA %>%
  ggplot(aes(PC1, PC2)) +
  geom_hline(yintercept=0, linetype=3) +
  geom_vline(xintercept=0, linetype=3) +
  geom_point(aes(fill=Timepoint), shape=shapes.treat["Control"], size=3) +
  scale_fill_manual(values=cols.dpi) +
  labs(x=pca$PCs[1],
       y=pca$PCs[2],
       title="PCA: controls")
ggsave("analysis/pca-controls.png",
       units="in", width=4, height=3)

# make model matrix 
modmat <- model.matrix(~0+m$Timepoint)
# format column names manually
colnames(modmat) <- c("Baseline", "DPI4", "Terminal")
# run differential expression comparisons by time point
rmat <- list(DPI4=diffexpr(c, modmat, 
                          makeContrasts(DPI4 - Baseline, 
                                        levels=colnames(modmat))) %>%
                  mutate(Timepoint="4 DPI"),
             Terminal=diffexpr(c, modmat, 
                               makeContrasts(Terminal - Baseline, 
                                             levels=colnames(modmat))) %>%
                      mutate(Timepoint="Terminal")) %>%
        # combine matrices
        do.call(rbind, .)
write.csv(rmat, "analysis/results-controls.csv", row.names=FALSE)

# total DEGs
rmat %>%
  filter(Significant) %>%
  group_by(Timepoint, Regulation) %>%
  summarise(Genes=n(),
            .groups="drop") %>%
  ggplot(aes(Timepoint, Genes)) +
  geom_col(aes(fill=Regulation), col="black", position="dodge") +
  geom_text(aes(y=Genes+2, label=Genes, group=Regulation), 
            position=position_dodge(width=1),
            hjust=0.5) +
  scale_fill_manual(values=cols.reg) +
  labs(y="Total DE genes",
       title="Total DE genes: controls") +
  theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave("analysis/degenes-controls.png",
       units="in", width=3, height=3)

# volcano
topgenes <- rmat %>%
            filter(Significant) %>%
            group_by(Timepoint, Regulation) %>%
            top_n(n=8, wt=abs(lfc)) %>%
            ungroup()
rmat %>%
  ggplot(aes(lfc, -log10(padj))) +
  geom_point(aes(fill=Regulation, col=Regulation, size=Regulation), alpha=0.5) +
  scale_fill_manual(values=cols.reg) +
  scale_color_manual(values=cols.reg) +
  scale_size_manual(values=c(1, 0.5, 1)) +
  ggrepel::geom_text_repel(data=topgenes, aes(label=Gene), 
                           force=80, size=3) +
  facet_wrap(~Timepoint) +
  labs(x="Fold change (log2)",
       y="FDR-adj. p-value (-log10)",
       title="DE genes: controls")
ggsave("analysis/volcano-controls.png",
       units="in", width=6.5, height=3)

# clean up
rm(c, m, modmat, pca, rmat, topgenes)

## vaccinated over time --------------------------------------------------------
# subset data and format times
m <- meta %>%
     filter(Treatment=="Vaccinated") %>%
     mutate(Timepoint=droplevels(Timepoint))
c <- cmat[, rownames(m)] %>%
     normalize.counts()

# PCA
pca <- c %>%
       cpm(log=TRUE) %>%
       run.pca(m)
pca$PCA %>%
  ggplot(aes(PC1, PC2)) +
  geom_hline(yintercept=0, linetype=3) +
  geom_vline(xintercept=0, linetype=3) +
  geom_point(aes(fill=Timepoint), shape=shapes.treat["Vaccinated"], size=3) +
  scale_fill_manual(values=cols.dpi) +
  labs(x=pca$PCs[1],
       y=pca$PCs[2],
       title="PCA: vaccinated")
ggsave("analysis/pca-vaccinated.png",
       units="in", width=4, height=3)

# make model matrix 
modmat <- model.matrix(~0+m$Timepoint)
# format column names manually
colnames(modmat) <- c("Baseline", "DPI4", "DPI7", "DPI10")
# run differential expression comparisons by time point
rmat <- list(DPI4=diffexpr(c, modmat, 
                           makeContrasts(DPI4 - Baseline, 
                                         levels=colnames(modmat))) %>%
                  mutate(Timepoint="4 DPI"),
             DPI7=diffexpr(c, modmat, 
                           makeContrasts(DPI7 - Baseline, 
                                         levels=colnames(modmat))) %>%
                  mutate(Timepoint="7 DPI"),
             DPI10=diffexpr(c, modmat, 
                           makeContrasts(DPI10 - Baseline, 
                                         levels=colnames(modmat))) %>%
                   mutate(Timepoint="10 DPI")) %>%
        # combine matrices
        do.call(rbind, .) %>%
        # format timepoint
        mutate(Timepoint=factor(Timepoint, levels=levels(m$Timepoint)))
write.csv(rmat, "analysis/results-vaccinated.csv", row.names=FALSE)

# total DEGs
rmat %>%
  filter(Significant) %>%
  group_by(Timepoint, Regulation) %>%
  summarise(Genes=n(),
            .groups="drop") %>%
  # plot it!
  ggplot(aes(Timepoint, Genes)) +
  geom_col(aes(fill=Regulation), col="black", position="dodge") +
  geom_text(aes(y=Genes+4, label=Genes, group=Regulation), 
            position=position_dodge(width=1),
            hjust=0.5) +
  scale_fill_manual(values=cols.reg) +
  labs(y="Total DE genes",
       title="Total DE genes: vaccinated") +
  theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave("analysis/degenes-vaccinated.png",
       units="in", width=4, height=3)

# volcano
topgenes <- rmat %>%
            filter(Significant) %>%
            group_by(Timepoint, Regulation) %>%
            top_n(n=8, wt=abs(lfc)) %>%
            ungroup()
rmat %>%
  ggplot(aes(lfc, -log10(padj))) +
  geom_point(aes(fill=Regulation, col=Regulation, size=Regulation), alpha=0.5) +
  scale_fill_manual(values=cols.reg) +
  scale_color_manual(values=cols.reg) +
  scale_size_manual(values=c(1, 0.5, 1)) +
  ggrepel::geom_text_repel(data=topgenes, aes(label=Gene), 
                           force=80, size=3) +
  facet_wrap(~Timepoint) +
  labs(x="Fold change (log2)",
       y="FDR-adj. p-value (-log10)",
       title="DE genes: vaccinated")
ggsave("analysis/volcano-vaccinated.png",
       units="in", width=6.5, height=3)

# clean up
rm(c, m, modmat, pca, rmat, topgenes)

## treatment comparison by timepoint -------------------------------------------
# we will go timepoint-by-timepoint with norm, PCA, and DE

# Baseline
m <- filter(meta, Timepoint=="Baseline")
c <- cmat[, rownames(m)] %>%
     normalize.counts()
# run PCA and save for later 
pca <- list(baseline=c %>%
                     cpm(log=TRUE) %>%
                     run.pca(m))
# make the matrix and run DE
modmat <- model.matrix(~0+m$Treatment)
colnames(modmat) <- levels(m$Treatment)
rmat <- list(Baseline=diffexpr(c, modmat,
                               makeContrasts(Vaccinated - Control,
                                             levels=colnames(modmat))) %>%
                      mutate(Timepoint="Baseline"))

# 4 DPI
m <- filter(meta, Timepoint=="4 DPI")
c <- cmat[, rownames(m)] %>%
     normalize.counts()
# run PCA and save for later 
pca$`4 DPI` <- c %>%
               cpm(log=TRUE) %>%
               run.pca(m)
# make the matrix and run DE
modmat <- model.matrix(~0+m$Treatment)
colnames(modmat) <- levels(m$Treatment)
rmat$DPI4 <- diffexpr(c, modmat,
                      makeContrasts(Vaccinated - Control,
                                    levels=colnames(modmat))) %>%
             mutate(Timepoint="4 DPI")

# 7 DPI/terminal
m <- filter(meta, Timepoint %in% c("7 DPI", "Terminal"))
c <- cmat[, rownames(m)] %>%
     normalize.counts()
# run PCA and save for later 
pca$`7 DPI/Terminal` <- c %>%
                        cpm(log=TRUE) %>%
                        run.pca(m)
# make the matrix and run DE
modmat <- model.matrix(~0+m$Treatment)
colnames(modmat) <- levels(m$Treatment)
rmat$DPI7 <- diffexpr(c, modmat,
                      makeContrasts(Vaccinated - Control,
                                    levels=colnames(modmat))) %>%
             mutate(Timepoint="7 DPI/Terminal")

# plot all PCAs
pca <- names(pca) %>%
       lapply(function(i) {
         # extract the object
         pca.time <- magrittr::extract2(pca, i)
         # plot it
         pca.time$PCA %>%
           ggplot(aes(PC1, PC2)) +
           geom_hline(yintercept=0, linetype=3) +
           geom_vline(xintercept=0, linetype=3) +
           geom_point(aes(fill=Treatment, shape=Treatment), size=3) +
           scale_fill_manual(values=cols.treat) +
           scale_shape_manual(values=shapes.treat) +
           labs(x=pca.time$PCs[1],
                y=pca.time$PCs[2],
                title=paste("PCA:", i))
       })
ggsave("analysis/pca-dpi00.png", plot=pca[[1]],
       units="in", width=4, height=3)
ggsave("analysis/pca-dpi04.png", plot=pca[[2]],
       units="in", width=4, height=3)
ggsave("analysis/pca-dpi07.png", plot=pca[[3]],
       units="in", width=4, height=3)

# combine all results and save
rmat <- rmat %>%
        do.call(rbind, .) %>%
        # re-define up- and down-regulation
        # controls are the denominator = "down"
        mutate(Regulation=factor(Regulation, 
                                 levels=names(cols.reg),
                                 labels=c("Vaccinated", "Neither", "Control")),
               Timepoint=factor(Timepoint, 
                                levels=c("Baseline", "4 DPI", 
                                         "7 DPI/Terminal")))
write.csv(rmat, "analysis/results-comparison.csv", row.names=FALSE)

# total DEGs
rmat %>%
  filter(Significant) %>%
  group_by(Timepoint, Regulation) %>%
  summarise(Genes=n(),
            .groups="drop") %>%
  # add in missing zeros
  right_join(expand.grid(Timepoint=levels(rmat$Timepoint),
                         Regulation=names(cols.treat)),
             by=c("Timepoint", "Regulation")) %>%
  replace_na(list(Genes=0)) %>%
  # plot it!
  ggplot(aes(Timepoint, Genes)) +
  geom_col(aes(fill=Regulation), col="black", position="dodge") +
  geom_text(aes(y=Genes+4, label=Genes, group=Regulation), 
            position=position_dodge(width=1),
            hjust=0.5) +
  scale_fill_manual(values=cols.treat) +
  labs(y="Total DE genes",
       title="Total DE genes: comparison",
       fill="Higher in") +
  theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave("analysis/degenes-comparison.png",
       units="in", width=4, height=3)

# volcano
topgenes <- rmat %>%
            filter(Significant) %>%
            group_by(Timepoint, Regulation) %>%
            top_n(n=8, wt=abs(lfc)) %>%
            ungroup()
rmat %>%
  ggplot(aes(lfc, -log10(padj))) +
  geom_point(aes(fill=Regulation, col=Regulation, size=Regulation), alpha=0.5) +
  scale_fill_manual(values=c(cols.treat, Neither="lightgrey")) +
  scale_color_manual(values=c(cols.treat, Neither="lightgrey")) +
  scale_size_manual(values=c(1, 0.5, 1)) +
  ggrepel::geom_text_repel(data=topgenes, aes(label=Gene), 
                           force=80, size=3) +
  facet_wrap(~Timepoint) +
  labs(x="Fold change (log2)",
       y="FDR-adj. p-value (-log10)",
       title="DE genes: comparison")
ggsave("analysis/volcano-comparison.png",
       units="in", width=6.5, height=3)

# clean up
rm(c, m, modmat, pca, rmat, topgenes)

## DCQ with CIBERSORTx ---------------------------------------------------------
# read in the results file and format
cib.sort <- read.csv("analysis/cibersort-output.csv", check.names=FALSE) %>%
            filter(`P-value` < 0.05) %>%
            select(-`P-value`, -Correlation, -RMSE, 
                   -`Absolute score (sig.score)`) %>%
            rename(ID=Mixture) %>%
            # convert to "long" format
            reshape2::melt(id.vars="ID",
                           variable.name="Celltype",
                           value.name="Score") %>%
            # add metadata
            left_join(meta, by="ID")

# remove cell types that aren't detected (0s only)
shortlist <- cib.sort %>%
             group_by(Celltype) %>%
             summarise(Upper=max(Score),
                       .groups="drop") %>%
             filter(Upper > 0) %>%
             select(Celltype) %>%
             unlist()
cib.sort <- cib.sort %>%
            filter(Celltype %in% shortlist)

# controls (n=3)
c <- filter(cib.sort, Treatment=="Control")
# do cell types vary over time? Yes, monocytes
p.fd <- c %>%
        group_by(Celltype) %>%
        rstatix::friedman_test(Score ~ Timepoint | NHP) %>%
        rstatix::p_format(new.col=TRUE, add.p=TRUE) %>%
        mutate(p.format=paste("Friedman", p.format)) %>%
        filter(p < 0.05)
# use Dunn's test for post-hoc. Likely would have been significant at 4 DPI,
# but we have a very low n
p.dn <- c %>%
        filter(Celltype %in% p.fd$Celltype) %>%
        group_by(Celltype) %>%
        rstatix::dunn_test(Score ~ Timepoint, p.adjust.method="fdr") %>%
        rstatix::p_format(new.col=TRUE, digits=2, add.p=TRUE) %>%
        rstatix::add_y_position(scales="free_y")
# plot it!
c %>%
  filter(Celltype %in% p.fd$Celltype) %>%
  ggplot(aes(Timepoint, Score)) +
  geom_boxplot(aes(group=Timepoint, fill=Timepoint), alpha=0.5) +
  scale_fill_manual(values=cols.dpi) +
  ggpubr::stat_pvalue_manual(data=p.dn, label="p.adj.format", 
                             step.increase=0.12) +
  geom_text(data=p.fd, aes(y=1.4, x=1, label=p.format), hjust=0) +
  facet_wrap(~Celltype, scales="free_y") +
  scale_y_continuous(expand=c(0.1, 0.2)) +
  labs(x="Time point",
       y="Cell type score",
       title="CIBERSORTx: controls") +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        legend.position="none")
ggsave("analysis/cibersort-controls.png",
       units="in", width=3, height=3)

# vaccinated (n=5)
c <- filter(cib.sort, Treatment=="Vaccinated")
# do cell types vary over time? Yes, multiple!
p.fd <- c %>%
        group_by(Celltype) %>%
        rstatix::friedman_test(Score ~ Timepoint | NHP) %>%
        rstatix::p_format(new.col=TRUE, add.p=TRUE, digits=1) %>%
        mutate(p.format=paste("Friedman", p.format)) %>%
        filter(p <= 0.05)
# use Dunn's test for post-hoc; only keep sig. dif. comparisons
p.dn <- c %>%
        filter(Celltype %in% p.fd$Celltype) %>%
        group_by(Celltype) %>%
        rstatix::dunn_test(Score ~ Timepoint, p.adjust.method="fdr") %>%
        filter(p.adj <= 0.05) %>%
        rstatix::add_y_position(scales="free_y", step.increase=0) %>%
        mutate(y.position=1.15*y.position)
# plot it!
c %>%
  filter(Celltype %in% p.dn$Celltype) %>%
  ggplot(aes(Timepoint, Score)) +
  geom_boxplot(aes(group=Timepoint, fill=Timepoint), alpha=0.5) +
  scale_fill_manual(values=cols.dpi) +
  ggpubr::stat_pvalue_manual(data=p.dn, label="p.adj.signif", hide.ns=TRUE, 
                             step.group.by="Celltype", step.increase=0.12) +
  facet_wrap(~Celltype, scales="free_y") +
  scale_y_continuous(expand=expansion(mult=c(0.05, 0.1))) +
  labs(x="Time point",
       y="Cell type score",
       title="DCQ: vaccinated") +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        legend.position="none")
ggsave("analysis/cibersort-vaccinated.png",
       units="in", width=6.5, height=4)

# all
# update time point levels to compare D7 and terminal
cib.sort <- cib.sort %>%
            filter(Timepoint != "10 DPI") %>%
            mutate(Timepoint=factor(Timepoint, 
                                    labels=c("Baseline", "4 DPI", 
                                             "7 DPI/Terminal", 
                                             "7 DPI/Terminal")))
# re-filter cell types since we are using a subset of the data
# remove cell types that aren't detected (0s only)
shortlist <- cib.sort %>%
             group_by(Celltype) %>%
             summarise(Upper=max(Score),
                       .groups="drop") %>%
             filter(Upper > 0) %>%
             select(Celltype) %>%
             unlist()
cib.sort <- cib.sort %>%
            filter(Celltype %in% shortlist)
# only use groups with at least one datapoint > 0 for stats
p.mw <- cib.sort %>%
        group_by(Celltype, Timepoint) %>%
        summarise(Upper=max(Score),
                  .groups="drop") %>%
        filter(Upper > 0) %>%
        select(-Upper) %>%
        # add back in the data, keeping only valid groups
        left_join(cib.sort, by=c("Celltype", "Timepoint")) %>%
        group_by(Celltype, Timepoint) %>%
        # run the test and format for plotting
        rstatix::wilcox_test(Score ~ Treatment) %>%
        rstatix::add_significance() %>%
        rstatix::add_xy_position(x="Timepoint")
# plot it! All cell types
cib.sort %>%
  ggplot(aes(Timepoint, Score)) +
  geom_boxplot(aes(group=interaction(Timepoint, Treatment), fill=Treatment), 
               alpha=0.5) +
  scale_fill_manual(values=cols.treat) +
  ggpubr::stat_pvalue_manual(data=p.mw, label="p.signif", hide.ns=TRUE) +
  scale_y_continuous(expand=expansion(mult=c(0.05, 0.1))) +
  facet_wrap(~Celltype, scales="free_y", ncol=3) +
  labs(x="Time point",
       y="Cell type score",
       title="DCQ: comparison") +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        legend.position=c(0.8, -0.05))
ggsave("analysis/cibersort-comparison-all.png",
       units="in", width=6.5, height=10)
# plot it! Significantly different cell types only
p.mw <- filter(p.mw, p < 0.05)
cib.sort %>%
  filter(Celltype %in% p.mw$Celltype) %>%
  ggplot(aes(Timepoint, Score)) +
  geom_boxplot(aes(group=interaction(Timepoint, Treatment), fill=Treatment), 
               alpha=0.5) +
  scale_fill_manual(values=cols.treat) +
  ggpubr::stat_pvalue_manual(data=p.mw, label="p.signif", hide.ns=TRUE) +
  scale_y_continuous(expand=expansion(mult=c(0.05, 0.1))) +
  facet_wrap(~Celltype, scales="free_y", ncol=3) +
  labs(x="Time point",
       y="Cell type score",
       title="DCQ: comparison") +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        legend.position=c(0.8, -0.07))
ggsave("analysis/cibersort-comparison-subset.png",
       units="in", width=6.5, height=4)

# clean up
rm(cib.sort, c, p.dn, p.mw, p.fd, shortlist)

## fin -------------------------------------------------------------------------
sessionInfo()
