#Figure 5
setwd("")

library(pheatmap)
library(tidyverse)
library(grid)
library(gridExtra)
library(svglite)

CELL_WIDTH  <- 10
CELL_HEIGHT <- 3

SVG_WIDTH_LOLLIPOP  <- 16
SVG_HEIGHT_LOLLIPOP <- 7

pal <- colorRampPalette(c(
  "#FFFFFF", "#5EABA5", "#82AE8A",
  "#A9B06D", "#D1AF58", "#E4A04F", "#E68656"
))(100)

df <- read.csv("0-other-res-df-all-visits-quartile-group.csv", check.names = FALSE)

make_drug_heatmap_grob <- function(df, anchor_col, target_cols,
                                   col_labels = NULL,
                                   title) {
  
  anchor   <- paste0("X", anchor_col)
  all_cols <- paste0("X", c(anchor_col, target_cols))
  
  sub <- df[df[[anchor]] != 0, c("physresilience_group", all_cols), drop = FALSE]
  cat(sprintf("\n[%s] anchor '%s' — %d patients retained\n", title, anchor, nrow(sub)))
  if (nrow(sub) == 0) { message("No rows. Skipping."); return(nullGrob()) }
  
  get_order <- function(mat, idx) {
    if (length(idx) <= 1) return(idx)
    sub_mat <- mat[idx, , drop = FALSE]
    keep    <- apply(sub_mat, 1, function(x) !all(is.na(x)))
    if (sum(keep) <= 1) return(idx)
    hc <- hclust(dist(sub_mat[keep, ], method = "euclidean"), method = "complete")
    idx[keep][hc$order]
  }
  
  mat_df              <- sub[, all_cols, drop = FALSE]
  mat_df[mat_df == 0] <- NA
  mat                 <- as.matrix(log10(mat_df))
  if (!is.null(col_labels)) colnames(mat) <- col_labels
  
  rn            <- paste0("r", seq_len(nrow(mat)))
  rownames(mat) <- rn
  
  ann_row <- data.frame(
    physresilience_group = sub$physresilience_group,
    row.names            = rn,
    stringsAsFactors     = FALSE
  )
  ann_colors <- list(
    physresilience_group = c(High = "#1E6BB5", Low = "#D09200")
  )
  
  high_idx  <- which(sub$physresilience_group == "High")
  low_idx   <- which(sub$physresilience_group == "Low")
  row_order <- c(get_order(mat, high_idx), get_order(mat, low_idx))
  
  val_max    <- ceiling(max(mat, na.rm = TRUE))
  if (is.infinite(val_max) || is.na(val_max)) val_max <- 6
  breaks     <- seq(0, val_max, length.out = 101)
  leg_breaks <- seq(0, val_max, by = 1)
  
  tmp <- pheatmap(
    mat[row_order, ],
    main                 = title,
    color                = pal,
    breaks               = breaks,
    na_col               = "white",
    cluster_rows         = FALSE,
    cluster_cols         = FALSE,
    show_rownames        = FALSE,
    show_colnames        = TRUE,
    fontsize_col         = 9,
    angle_col            = 90,
    border_color         = NA,
    annotation_row       = ann_row[row_order, , drop = FALSE],
    annotation_colors    = ann_colors,
    annotation_names_row = FALSE,
    legend               = TRUE,
    legend_breaks        = leg_breaks,
    legend_labels        = as.character(leg_breaks),
    cellwidth            = CELL_WIDTH,
    cellheight           = CELL_HEIGHT,
    silent               = TRUE
  )
  
  return(tmp$gtable)
}

grob_verapamil <- make_drug_heatmap_grob(
  df,
  anchor_col  = "70755",
  target_cols = c("53885", "69728", "66093", "52678", "46673", "52067", "67307", "51479", "35698", "55248"),
  col_labels  = c("Verapamil",
                  "Dealkylation",
                  "Demethylation",
                  "2xDemethylation",
                  "Dealkylation + demethylation (rt: 4.04 min)",
                  "Dealkylation + demethylation (rt: 3.83 min)",
                  "Demethylation + glucuronidation",
                  "Hydroxylation",
                  "2xDemethylation + glucuronidation",
                  "Dealkylation + demethylation + glucuronidation",
                  "4xDemethylation"),
  title       = "Verapamil"
)

grob_metoprolol <- make_drug_heatmap_grob(
  df,
  anchor_col  = "38586",
  target_cols = c("18674", "22348", "27538",
                  "6472",  "16443", "25778", "21293",
                  "20824", "35588"),
  col_labels  = c("Metoprolol",
                  "Metoprolol benzoic acid",
                  "Demethylmetoprolol",
                  "Metoprolol acid",
                  "Demethylhydroxymetoprolol",
                  "Hydroxylmetoprolol acid",
                  "Hydroxymetoprolol (rt: 2.92 min)",
                  "Hydroxymetoprolol (rt: 2.62 min)",
                  "Hydroxylmetoprolol glucuronide",
                  "Metoprolol glucuronide"),
  title       = "Metoprolol"
)

n_cols_vera <- 11
n_cols_meto <- 10

w_vera  <- (n_cols_vera * CELL_WIDTH / 72) + 4
w_meto  <- (n_cols_meto * CELL_WIDTH / 72) + 4
total_w <- w_vera + w_meto
total_h <- 14   # adjust if rows are clipped

svglite("heatmap.svg", width = total_w, height = total_h)
grid.arrange(
  grob_verapamil,
  grob_metoprolol,
  ncol   = 2,
  widths = c(w_vera, w_meto)
)
dev.off()
cat(sprintf("heatmap.svg  (%.1f x %.1f in)\n", total_w, total_h))

get_log2fc <- function(df, anchor_col, target_cols, meta_labels, drug_name) {
  anchor  <- paste0("X", anchor_col)
  targets <- paste0("X", target_cols)
  
  sub <- df[df[[anchor]] != 0, c("physresilience_group", anchor, targets), drop = FALSE]
  cat(sprintf("\n[%s] anchor '%s' — %d patients retained\n", drug_name, anchor, nrow(sub)))
  if (nrow(sub) == 0) return(NULL)
  
  result <- lapply(seq_along(targets), function(i) {
    vals_high <- sub[[targets[i]]][sub$physresilience_group == "High"]
    vals_low  <- sub[[targets[i]]][sub$physresilience_group == "Low"]
    
    # Skip if either group has fewer than 3 non-zero values
    if (sum(vals_high != 0, na.rm = TRUE) < 3 |
        sum(vals_low  != 0, na.rm = TRUE) < 3) {
      cat(sprintf("  [SKIP] %s — High: %d, Low: %d non-zero values\n",
                  meta_labels[i],
                  sum(vals_high != 0, na.rm = TRUE),
                  sum(vals_low  != 0, na.rm = TRUE)))
      return(NULL)
    }
    
    ratio     <- sub[[targets[i]]] / sub[[anchor]]
    mean_high <- mean(ratio[sub$physresilience_group == "High"], na.rm = TRUE)
    mean_low  <- mean(ratio[sub$physresilience_group == "Low"],  na.rm = TRUE)
    log2fc    <- log2(mean_high / mean_low)
    
    cat(sprintf("  %-50s log2FC = %6.3f\n", meta_labels[i], log2fc))
    
    data.frame(
      drug       = drug_name,
      metabolite = meta_labels[i],
      log2fc     = log2fc
    )
  })
  
  bind_rows(result)
}

df_verapamil <- get_log2fc(
  df, "70755",
  c("53885","69728","66093","52678","46673","52067","67307","51479","35698","55248"),
  c("Dealkylation","Demethylation","2xDemethylation",
    "Dealkylation + demethylation (rt: 4.04 min)",
    "Dealkylation + demethylation (rt: 3.83 min)",
    "Demethylation + glucuronidation","Hydroxylation",
    "2xDemethylation + glucuronidation",
    "Dealkylation + demethylation + glucuronidation","4xDemethylation"),
  "Verapamil"
)

df_metoprolol <- get_log2fc(
  df, "38586",
  c("18674","22348","27538","6472","16443","25778","21293","20824","35588"),
  c("Metoprolol benzoic acid","Demethylmetoprolol","Metoprolol acid",
    "Demethylhydroxymetoprolol","Hydroxylmetoprolol acid",
    "Hydroxymetoprolol (rt: 2.92 min)","Hydroxymetoprolol (rt: 2.62 min)",
    "Hydroxylmetoprolol glucuronide","Metoprolol glucuronide"),
  "Metoprolol"
)

# ============================================================
# Combine results and prepare for plotting
# ============================================================
plot_df <- bind_rows(df_verapamil, df_metoprolol) %>%
  mutate(
    log2fc    = ifelse(is.finite(log2fc), log2fc, NA),
    direction = case_when(
      log2fc >= 0 ~ "Higher in high resilience",
      log2fc <  0 ~ "Higher in low resilience",
      TRUE        ~ NA_character_
    ),
    drug = factor(drug, levels = c("Verapamil", "Metoprolol"))
  ) %>%
  group_by(drug) %>%
  arrange(drug, desc(log2fc)) %>%
  mutate(metabolite = factor(metabolite, levels = unique(metabolite))) %>%
  ungroup()

print(plot_df)
cat(sprintf("\nRows in plot_df: %d\n", nrow(plot_df)))

if (nrow(plot_df) == 0) {
  message("plot_df is empty — no metabolites passed the n >= 3 filter. Check [SKIP] messages above.")
} else {
  
  y_range <- max(abs(plot_df$log2fc), na.rm = TRUE) * 1.2
  
  p_lollipop <- ggplot(plot_df, aes(x = metabolite, y = log2fc, color = direction)) +
    geom_segment(aes(x = metabolite, xend = metabolite, y = 0, yend = log2fc),
                 color = "grey30", linewidth = 0.6, na.rm = TRUE) +
    geom_point(size = 4, na.rm = TRUE) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    scale_color_manual(
      values = c("Higher in high resilience" = "#1E6BB5",
                 "Higher in low resilience"  = "#D09200"),
      name = "", na.translate = FALSE
    ) +
    facet_grid(. ~ drug, scales = "free_x", space = "free_x") +
    labs(
      x = "Metabolite",
      y = expression(log[2](High / Low))
    ) +
    ylim(-y_range, y_range) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid       = element_blank(),
      panel.border     = element_blank(),
      axis.line.x      = element_line(color = "black"),
      axis.line.y      = element_line(color = "black"),
      axis.text.x      = element_text(angle = 45, hjust = 1),
      legend.position  = "top",
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text       = element_text(face = "bold", size = 11),
      plot.margin      = margin(t = 10, r = 30, b = 10, l = 10, unit = "pt")
    )
  
  # -- Save as SVG --
  ggsave(
    filename = "lollipop_log2fc.svg",
    plot     = p_lollipop,
    width    = SVG_WIDTH_LOLLIPOP,
    height   = SVG_HEIGHT_LOLLIPOP,
    units    = "in"
  )
  cat("Saved: lollipop_log2fc.svg\n")
}

