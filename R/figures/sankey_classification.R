#!/usr/bin/env Rscript
# sankey_classification.R
#
# Gene-classification Sankey (Figure S1 of Purzner et al. 2025, eLife
# 10.7554/eLife.105339.1). Shows how every transcript flows through the
# stepwise filter that defines the three developmental gene programs:
#
#   Step 1  expressed during the time course?            (expr)
#   Step 2  granule-neuron lineage vs cell-cycle /        (group_1, from Taylor
#           alternative-lineage contamination              scRNA-seq)
#   Step 3  dynamically expressed (>=2-fold)?              (avg_fc)
#   Step 4  time of peak expression: E15 / P0-P7 / P14-P56 (avg_max)
#   Step 5  per-dataset differential-expression p-value    (max_cat)
#   Step 6  drop late genes that rise mostly before P0      (t50_cat)
#   Step 7  separate by fold change                         (avg_fc)
#   Step 8  final groups: embryonic / proliferation / differentiation genes
#
# Ported verbatim (logic preserved) from
# GNP_MB_integrate/Ezh2_2022/3_gene_classification.Rmd. See analysis/ for the
# original Rmd.
#
# INPUT : data/gene_classification.tsv  (toc_t50 kinetics + scRNA-seq lineage
#         annotation; see data/README.md)
# OUTPUT: figures/sankey_gene_classification.pdf (+ .png)
#
# Requires ggsankey (GitHub only):
#   devtools::install_github("davidsjoberg/ggsankey")
#
# Usage:  Rscript R/figures/sankey_classification.R

suppressWarnings(suppressMessages({
  library(dplyr); library(forcats); library(stringr); library(ggplot2); library(ggsankey)
  args <- commandArgs(trailingOnly = FALSE)
  here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])))
  source(file.path(here, "..", "lib", "common.R"))
}))

cls <- read.table(file.path(MT_DATA, "gene_classification.tsv"),
                  sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                  quote = "", comment.char = "")

# max_cat_exclude: flags genes whose min falls in a window that makes the
# max-category call incompatible (self-contained derivation on max_cat/avg_min).
cls <- cls %>% dplyr::mutate(max_cat_exclude = case_when(
  max_cat == "E15"     &  avg_min >  6                   ~ "E15",
  max_cat == "E15"     &  avg_min <= 6                   ~ "E15_exclude",
  max_cat == "P0_7"    & (avg_min <  6 | avg_min > 20)   ~ "P0_P7",
  max_cat == "P0_7"    & (avg_min >= 6 | avg_min <= 20)  ~ "P0_P7_exclude",
  max_cat == "P14_P56" &  avg_min <  19                  ~ "P14_P56",
  max_cat == "P14_P56" &  avg_min >= 19                  ~ "P14_P56_exclude",
  TRUE                                                   ~ "other"))

# ---- build the step columns (step1..step7) --------------------------------
long_src <- cls %>%
  dplyr::filter(!is.na(expr)) %>%
  dplyr::mutate(max_cat2 = ifelse(max_cat == "unexpressed", "excluded", max_cat),
                max_cat2 = ifelse(max_cat2 == "P14_P56", "P14 to P56", max_cat2),
                max_cat2 = ifelse(max_cat2 == "P0_7", "P0 to P7", max_cat2)) %>%
  dplyr::mutate(step1 = recode(factor(expr), "TRUE" = "Expressed", "FALSE" = "Unexpressed")) %>%
  dplyr::mutate(step2 = ifelse(expr == FALSE, "excluded", group_1),
                step2 = ifelse(group_1 %in% c("Dynically expressed in granule neuron lineage",
                                               "Unchanged expression"), "Granule neuron lineage", step2),
                step2 = ifelse(max_cat == "unexpressed", "excluded", step2),
                step2 = recode(step2, "Granule neuron lineage" = "Granule neuron lineages",
                               "Contaminated genes" = "Alternative lineage",
                               "Cell cycle" = "Cell cycle", "excluded" = "excluded")) %>%
  dplyr::mutate(step3 = ifelse(avg_fc >= 1, "Dynamically expressed", "unchanged"),
                step3 = ifelse(expr == FALSE, "excluded", step3),
                step3 = ifelse(group_1 %in% c("Cell cycle", "Contaminated genes", "Unexpressed"),
                               "excluded", step3)) %>%
  dplyr::mutate(step4 = as.character(cut(avg_max, c(0,2,20,64), right = FALSE,
                                         labels = c("E15","P0 to P7","P14 to P56"))),
                step4 = ifelse(step3 == "Dynamically expressed", step4, step3),
                step4 = ifelse(expr == FALSE, "excluded", step4),
                step4 = ifelse(grepl("_exclude", max_cat_exclude), "incompatible max-min", step4),
                step4 = ifelse(step3 == "excluded", "excluded", step4)) %>%
  dplyr::mutate(step5 = ifelse(!(exclude_cat %in% c("other", "GN_expressed_not_changing")),
                               "excluded", max_cat2),
                step5 = ifelse(step4 %in% c("E15","P0 to P7","P14 to P56") & max_cat2 == "unchanged",
                               "p-val inadequate", step4),
                step5 = ifelse(grepl("_exclude", max_cat_exclude), "excluded", step5)) %>%
  dplyr::mutate(step5b = ifelse(t50_cat == "Early" & max_cat2 == "P14 to P56",
                                "increase before P0", step5),
                step5b = ifelse(step5b == "p-val inadequate", "unchanged", step5b),
                step5b = ifelse(expr == FALSE, "excluded", step5b),
                step5b = ifelse(group_1 %in% c("Cell cycle","Contaminated genes","Unexpressed"),
                                "excluded", step5b),
                step5b = ifelse(grepl("_exclude", max_cat_exclude), "excluded", step5b)) %>%
  dplyr::mutate(step6 = case_when(avg_fc < 1 ~ "0-1 log2FC",
                                  avg_fc >= 1 & avg_fc < 2 ~ "1-2 log2FC",
                                  avg_fc >= 2 ~ ">2 log2FC"),
                step6 = paste(step5, step6),
                step6 = ifelse(step5 %in% c("E15","P0 to P7","P14 to P56"), step6, "excluded"),
                step6 = ifelse(step5b == "increase before P0", "excluded", step6)) %>%
  dplyr::mutate(step7 = case_when(
    step6 == "P14 to P56 >2 log2FC" ~ "differentiation genes",
    step6 == "P0 to P7 >2 log2FC"   ~ "proliferation genes",
    step6 == "E15 >2 log2FC"        ~ "embryonic genes"))

long_sce <- ggsankey::make_long(long_src, step1, step2, step3, step4, step5, step5b, step6, step7)
dagg <- long_sce %>% dplyr::group_by(x, node) %>% dplyr::tally()
long_sce <- long_sce %>% dplyr::left_join(dagg, by = c("x" = "x", "node" = "node"))

# ---- node order + colours (stable mapping across all nodes) ---------------
re_order_sankey <- rev(c("Expressed","Granule neuron lineages","Dynamically expressed",
  "P14 to P56","P0 to P7","E15","Cell cycle","Alternative lineage","increase before P0",
  "incompatible max-min","p-val inadequate","P14 to P56 >2 log2FC","P14 to P56 1-2 log2FC",
  "P0 to P7 >2 log2FC","P0 to P7 1-2 log2FC","E15 >2 log2FC","E15 1-2 log2FC",
  "differentiation genes","proliferation genes","embryonic genes","unchanged","excluded",
  "Unexpressed", NA))
long_sce$node <- factor(long_sce$node, levels = re_order_sankey, ordered = TRUE)
long_sce <- long_sce %>% dplyr::group_by(node) %>%
  dplyr::arrange(desc(node), desc(next_node)) %>% dplyr::filter(!is.na(node))

sankey_colours <- c("white","gray40","gray60","lightblue4","springgreen4","red3","royalblue3",
  "springgreen1","springgreen4","red1","red3","royalblue1","royalblue3","lightblue2","lightblue2",
  "lightblue2","mediumpurple","mediumorchid","springgreen3","red2","royalblue2","dodgerblue3",
  "dodgerblue2","dodgerblue")
names(sankey_colours) <- re_order_sankey

sankey_x_labs_nl <- c(
  "Step 1:\nRemove genes \nunexpressed during \nthe time course",
  "Step 2:\nRemove genes \ncontaminated\nfrom other cell lines\nTaylor scRNAseq",
  "Step 3:\nRemove genes that\nchange less\nthen 2 fold",
  "Step 4:\nSplit changing genes\nbased on when they\narrive at max",
  "Step 5:\nRemove genes with \nadjusted p-val > 0.001 \nin differential expression \nfrom RNAseq or NeuroD1-TRAP",
  "Step 6:\nRemove differentiation\ngenes that\nincrease most\nprior to P0",
  "Step 7:\nSeperate\nby fold change",
  "Step 8:\nFinal groups")

fc_cat_labeller <- function(var) {
  var <- str_replace(var, "P14 to P56 ", ""); var <- str_replace(var, "P0 to P7 ", "")
  str_replace(var, "E15 ", "")
}

p <- long_sce %>%
  ggplot(aes(x = x, next_x = next_x, node = node, next_node = next_node,
             fill = node, label = paste0(fc_cat_labeller(node), "\nn=", n))) +
  geom_sankey(flow.alpha = .6, width = 0.25) +
  geom_sankey_text(size = 3, color = "black") +
  theme_sankey(base_size = 16) +
  scale_fill_manual(values = sankey_colours) +
  scale_x_discrete(labels = sankey_x_labs_nl) +
  xlab(NULL) +
  theme(legend.position = "none", axis.text.x = element_text(size = 9))

dir.create(MT_FIGURES, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(MT_FIGURES, "sankey_gene_classification.pdf"), p,
       width = 13, height = 7, bg = "white")
ggsave(file.path(MT_FIGURES, "sankey_gene_classification.png"), p,
       width = 13, height = 7, dpi = 120, bg = "white")
message("Wrote figures/sankey_gene_classification.{pdf,png}")
message("Final groups:")
print(table(long_src$step7, useNA = "no"))
