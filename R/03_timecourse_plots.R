#!/usr/bin/env Rscript
# 03_timecourse_plots.R
#
# Generate the standard timecourse figures from the shipped, precomputed
# artifacts (data/log2_norm_list_shift.RData + data/timecourse_query_index.tsv).
# No re-run of the merge is required.
#
# Outputs (into figures/):
#   kinetic_class_timelines.pdf   mean +/- SD trajectory of every changing
#                                 transcript grouped by kinetic class
#                                 (Max_* t50_*), both datasets overlaid, with
#                                 the median t50 marker per group.
#   example_genes_timeline.pdf    a small panel of illustrative marker genes.
#   t50_vs_log2fc.pdf             scatter of t50 vs signed fold change,
#                                 faceted by time-of-max category.
#
# Usage:  Rscript R/02_timecourse_plots.R

suppressWarnings(suppressMessages({
  library(ggplot2)
  library(dplyr)
  args <- commandArgs(trailingOnly = FALSE)
  here <- dirname(normalizePath(sub("^--file=", "",
            grep("^--file=", args, value = TRUE)[1])))
  source(file.path(here, "lib", "common.R"))
  mt_source_lib()
  source(file.path(MT_LIB, "mt_plot.R"))
}))

dir.create(MT_FIGURES, showWarnings = FALSE, recursive = TRUE)

expr <- mt_load_expr()
idx  <- mt_load_index()

# ---------------------------------------------------------------------------
# 1. Kinetic-class timelines: every changing transcript grouped by its
#    combined time-of-max / t50 class. This is the signature figure.
# ---------------------------------------------------------------------------
message("Figure 1/3: kinetic_class_timelines.pdf")

changing <- idx %>%
  filter(!kinetic_class %in% c("unexpressed", "unchanged"),
         !is.na(kinetic_class),
         transcript_id %in% rownames(expr[[1]]))

# Order classes by developmental time-of-max then t50 for a readable grid.
class_order <- c(
  "Max_E15 t50_Early",     "Max_E15 t50_Mid",     "Max_E15 t50_Late",
  "Max_P0_7 t50_Early",    "Max_P0_7 t50_Mid",    "Max_P0_7 t50_Late",
  "Max_P14_P56 t50_Early", "Max_P14_P56 t50_Mid", "Max_P14_P56 t50_Late"
)
class_order <- class_order[class_order %in% unique(changing$kinetic_class)]

g2c <- factor(changing$kinetic_class, levels = class_order)
names(g2c) <- changing$transcript_id

t50_tbl <- data.frame(gene = changing$transcript_id, avg_t50 = changing$t50)

p1 <- plot_group_timeline(
  expr_list = expr, meta_list = mt_meta_list(), gene2cluster = g2c,
  p_time = MT_PSEUDO_TIME, ylimits = c(-2.5, 2.5), include_n = TRUE,
  col_num = 3, text_size = 13, x_label = "remove",
  add_t50 = t50_tbl, time_ref = mt_time_ref()
)
ggsave(file.path(MT_FIGURES, "kinetic_class_timelines.pdf"), p1,
       width = 11, height = 10, limitsize = FALSE, bg = "white")
ggsave(file.path(MT_FIGURES, "kinetic_class_timelines.png"), p1,
       width = 11, height = 10, dpi = 110, limitsize = FALSE, bg = "white")

# ---------------------------------------------------------------------------
# 2. Example marker genes across the differentiation axis.
# ---------------------------------------------------------------------------
message("Figure 2/3: example_genes_timeline.pdf")

examples <- c("Ezh2", "Sox2", "Gli1", "Atoh1", "Neurod1", "Gabra6")
examples <- examples[examples %in% idx$gene]
p2 <- mt_plot_genes(examples, expr = expr, idx = idx,
                    ylimits = c(-3, 3), col_num = 3, text_size = 13)
ggsave(file.path(MT_FIGURES, "example_genes_timeline.pdf"), p2,
       width = 10, height = 7, limitsize = FALSE, bg = "white")
ggsave(file.path(MT_FIGURES, "example_genes_timeline.png"), p2,
       width = 10, height = 7, dpi = 110, limitsize = FALSE, bg = "white")

# ---------------------------------------------------------------------------
# 3. t50 vs signed fold change, faceted by time-of-max category.
# ---------------------------------------------------------------------------
message("Figure 3/3: t50_vs_log2fc.pdf")

sc <- idx %>%
  filter(expressed, !max_category %in% c("unexpressed", "unchanged"),
         log2fc > 1, !is.na(t50)) %>%
  mutate(max_category = factor(max_category, levels = c("E15", "P0_7", "P14_P56")),
         t50_disp = pmin(t50, 30))

p3 <- ggplot(sc, aes(x = t50_disp, y = signed_log2fc)) +
  geom_point(size = 0.4, alpha = 0.4) +
  facet_wrap(~max_category) +
  labs(x = "t50 (pseudo-time index, capped at 30)",
       y = "signed log2 fold change",
       title = "Transcript kinetics by time of peak expression") +
  theme_bw(base_size = 13)
ggsave(file.path(MT_FIGURES, "t50_vs_log2fc.pdf"), p3, width = 10, height = 4, bg = "white")

message("Done. Figures written to: ", MT_FIGURES)
