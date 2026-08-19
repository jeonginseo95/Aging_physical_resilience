### Figure 5


library(ggplot2)
library(vegan)
library(dplyr)
library(lme4)
library(lmerTest)

BASE_DIR <- "C:/Users/sirju/OneDrive/Desktop/AgingPaper2/Tryptophan_metabolites"
OUT_DIR  <- file.path(BASE_DIR, "figures")
dir.create(OUT_DIR, showWarnings = FALSE)

ALL_VISITS <- "0-other-res-df-all-visits.csv"
GROUPS     <- "0-other-res-df-group.csv"
GNPS       <- "885c8595f0574824a2bdef9d6ee2db88-merged_results_with_gnps.tsv"

ID <- "record_ID"; AGE <- "host_age"; GROUP <- "physresilience_group"
MIN_MQSCORE   <- 0.70
MIN_DETECTION <- 0.70
FOCUS <- "Kynurenine"      # compound shown in the scatter and slope figures

PATHWAY <- paste0("tryptophan|tryptamin|kynuren|quinolinic|xanthurenic|picolinic|",
                  "anthranil|indol|indox|serotonin|melatonin|hydroxyindole|",
                  "nicotinamide|niacin|quinaldic|skatole")

COMPOUND_LABELS   <- c("L-Tryptophan", "Kynurenine", "Indole-3-lactate",
                       "Glucopyranosyl-L-tryptophan",
                       "5-Methylindole-3-carboxaldehyde")
COMPOUND_PATTERNS <- c("^l-tryptophan", "kynurenine",
                       "indole-?3-?lactic|indolelactic",
                       "glucopyranosyl.*tryptophan",
                       "methylindole-?3-?carboxaldehyde")

BLUE <- "#2567B5"; RED <- "#B87400"   # high resilience, low resilience
LEGEND <- c(High = "High physical resilience", Low = "Low physical resilience")

d <- read.csv(file.path(BASE_DIR, ALL_VISITS), check.names = FALSE)

meta_cols <- c(ID, "visit", "filename", AGE)
feat_cols <- setdiff(colnames(d), meta_cols)
feat_cols <- feat_cols[grepl("^[0-9]+$", feat_cols)]

ref <- read.csv(file.path(BASE_DIR, GROUPS), check.names = FALSE)
d[[GROUP]] <- ref[[GROUP]][match(d[[ID]], ref[[ID]])]
d <- d[!is.na(d[[GROUP]]), ]

stopifnot(all(tapply(d[[GROUP]], d[[ID]], function(x) length(unique(x))) == 1))

d$age_within <- d[[AGE]] - ave(d[[AGE]], d[[ID]], FUN = mean)
d$low        <- as.integer(d[[GROUP]] == "Low")

cat(nrow(d), "visits,", length(unique(d[[ID]])), "participants\n")
print(table(unique(d[, c(ID, GROUP)])[[GROUP]]))

feat <- as.matrix(d[, feat_cols])
colnames(feat) <- feat_cols
storage.mode(feat) <- "double"

cat("\nfeature matrix:", nrow(feat), "x", ncol(feat),
    " zeros:", sum(feat == 0), "\n")

rc <- as.matrix(vegan::decostand(feat, method = "rclr"))
colnames(rc) <- feat_cols
rc[feat == 0] <- NA

cat("non-missing rclr cells:", sum(!is.na(rc)), "of", length(rc), "\n")
cat("rclr range:", round(min(rc, na.rm = TRUE), 3), "to",
    round(max(rc, na.rm = TRUE), 3), "\n")

ex <- which(d[[ID]] == d[[ID]][1])
cat("\nexample participant", d[[ID]][1], "- kynurenine (scan 8541):\n")
print(data.frame(visit = d$visit[ex],
                 raw   = feat[ex, "8541"],
                 rclr  = round(rc[ex, "8541"], 4)), row.names = FALSE)


# ---------------------------------------------------------------- panel ------
lib <- read.delim(file.path(BASE_DIR, GNPS), check.names = FALSE)
scan_col <- grep("Scan", colnames(lib), value = TRUE)[1]
lib$scan <- as.character(lib[[scan_col]])

lib <- lib[order(-lib$MQScore), ]
lib <- lib[!duplicated(lib$scan), ]
lib <- lib[lib$MQScore >= MIN_MQSCORE, ]

hits <- lib[grepl(PATHWAY, tolower(lib$Compound_Name)) & lib$scan %in% feat_cols, ]
hits$detection <- colMeans(feat[, hits$scan, drop = FALSE] > 0)

name_of <- function(x) {
  x <- tolower(x)
  for (i in seq_along(COMPOUND_PATTERNS))
    if (grepl(COMPOUND_PATTERNS[i], x)) return(COMPOUND_LABELS[i])
  substr(x, 1, 40)
}
hits$compound <- vapply(hits$Compound_Name, name_of, character(1))

panel <- hits[hits$detection >= MIN_DETECTION, ]
panel <- panel[order(panel$compound, -panel$detection, -panel$MQScore), ]
panel <- panel[!duplicated(panel$compound), ]
panel <- panel[order(-panel$detection), ]

cat("\ncompounds tested:", nrow(panel), "\n")
print(data.frame(compound = panel$compound, scan = panel$scan,
                 detection = round(panel$detection, 3)), row.names = FALSE)

frames <- vector("list", nrow(panel))
names(frames) <- panel$compound

cat("\nfitting models\n")
for (i in seq_len(nrow(panel))) {
  df <- data.frame(y      = rc[, panel$scan[i]],
                   age    = d[[AGE]],
                   within = d$age_within,
                   group  = d[[GROUP]],
                   low    = d$low,
                   id     = d[[ID]],
                   stringsAsFactors = FALSE)
  df <- df[stats::complete.cases(df), ]
  df$id <- factor(df$id)
  
  n_obs <- nrow(df); n_id <- nlevels(df$id)
  cat(sprintf("  %-34s obs %4d  participants %4d  visits/participant %.2f\n",
              panel$compound[i], n_obs, n_id, n_obs / n_id))
  if (n_id >= n_obs)
    stop("only one visit per participant for ", panel$compound[i],
         " - check that rc rows line up with d rows")
  
  frames[[i]] <- list(model = lmer(y ~ within * low + (1 | id), data = df),
                      data  = df)
}

res <- do.call(rbind, lapply(seq_along(frames), function(i) {
  cf <- summary(frames[[i]]$model)$coefficients
  data.frame(
    compound     = panel$compound[i],
    scan         = panel$scan[i],
    detection_pc = round(100 * panel$detection[i]),
    n_visits     = nrow(frames[[i]]$data),
    slope_High   = cf["within", "Estimate"],
    slope_Low    = cf["within", "Estimate"] + cf["within:low", "Estimate"],
    slope_diff   = cf["within:low", "Estimate"],
    slope_SE     = cf["within:low", "Std. Error"],
    p            = cf["within:low", "Pr(>|t|)"],
    singular     = isSingular(frames[[i]]$model),
    stringsAsFactors = FALSE
  )
}))

res$FDR <- p.adjust(res$p, method = "BH")
res <- res[order(res$p), ]

cat("\nresults\n"); print(res, row.names = FALSE)
write.csv(res, file.path(OUT_DIR, "tryptophan_results.csv"), row.names = FALSE)


# forest plot
forest_df <- res
forest_df$compound <- factor(forest_df$compound,
                             levels = forest_df$compound[order(forest_df$slope_diff)])
forest_df$sig   <- forest_df$FDR < 0.05
forest_df$label <- ifelse(forest_df$sig, sprintf("FDR = %.3f", forest_df$FDR), "")

p_forest <- ggplot(forest_df, aes(x = slope_diff, y = compound)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbarh(aes(xmin = slope_diff - 1.96 * slope_SE,
                     xmax = slope_diff + 1.96 * slope_SE,
                     colour = sig), height = 0.18, linewidth = 0.7) +
  geom_point(aes(colour = sig, size = sig)) +
  geom_text(aes(x = slope_diff + 1.96 * slope_SE, label = label),
            hjust = -0.12, size = 3.2, colour = "#2B2B2B") +
  scale_colour_manual(values = c(`TRUE` = "#2B2B2B", `FALSE` = "#B4BAC2")) +
  scale_size_manual(values = c(`TRUE` = 2.6, `FALSE` = 2.0)) +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.30))) +
  labs(x = expression(Delta ~ "age slope (low - high)"), y = NULL,
       title = "Tryptophan-pathway panel") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 10.5),
        plot.title  = element_text(size = 11.5, hjust = 0))

ggsave(file.path(OUT_DIR, "forest_age_slope.svg"), p_forest,
       device = "svg", width = 6.4, height = 0.6 * nrow(res) + 2.0)


# scatter plot
fd   <- frames[[FOCUS]]$data
frow <- res[res$compound == FOCUS, ]

pt         <- unique(fd[, c("id", "group")])
n_by_group <- table(pt$group)
labs_group <- setNames(sprintf("%s (n = %d)", LEGEND[names(n_by_group)],
                               as.integer(n_by_group)), names(n_by_group))

p_scatter <- ggplot(fd, aes(x = age, y = y, colour = group, fill = group)) +
  geom_point(alpha = 0.28, size = 1.1, stroke = 0) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              linewidth = 1.2, alpha = 0.30) +
  scale_colour_manual(values = c(High = BLUE, Low = RED), labels = labs_group) +
  scale_fill_manual(values   = c(High = BLUE, Low = RED), labels = labs_group) +
  labs(x = "Chronological age (years)",
       y = paste0(FOCUS, " (rclr)"),
       title = sprintf("%s\nFDR = %.3f", FOCUS, frow$FDR)) +
  theme_classic(base_size = 13) +
  theme(legend.position = c(0.98, 0.02), legend.justification = c(1, 0),
        legend.title = element_blank(), legend.background = element_blank(),
        plot.title = element_text(size = 11.5, hjust = 0))

ggsave(file.path(OUT_DIR, paste0(gsub("[^A-Za-z0-9]+", "_", FOCUS), "_by_age.svg")),
       p_scatter, device = "svg", width = 5.6, height = 4.8)


# box plot
keep <- names(which(table(fd$id) >= 3))
sl <- lapply(keep, function(pid) {
  s <- fd[fd$id == pid, ]
  if (stats::sd(s$within) == 0) return(NULL)
  data.frame(id    = pid,
             group = s$group[1],
             slope = unname(coef(lm(y ~ within, data = s))[2]),
             stringsAsFactors = FALSE)
})
slopes <- do.call(rbind, sl)

wt <- wilcox.test(slope ~ group, data = slopes)
cat(sprintf("\nper-participant slope test: High median %.4f, Low median %.4f, Mann-Whitney P = %.4f\n",
            median(slopes$slope[slopes$group == "High"]),
            median(slopes$slope[slopes$group == "Low"]), wt$p.value))

slope_n   <- table(slopes$group)
slope_lab <- setNames(sprintf("%s\n(n=%d)", names(slope_n), as.integer(slope_n)),
                      names(slope_n))

p_slopes <- ggplot(slopes, aes(x = group, y = slope)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
  geom_jitter(aes(colour = group), width = 0.14, alpha = 0.55,
              size = 1.4, stroke = 0) +
  geom_boxplot(fill = NA, colour = "#2B2B2B", width = 0.55,
               outlier.shape = NA, linewidth = 0.5) +
  scale_colour_manual(values = c(High = BLUE, Low = RED)) +
  scale_x_discrete(labels = slope_lab) +
  labs(x = NULL, y = "Individual slope",
       title = sprintf("%s\nMann-Whitney P = %.4f", FOCUS, wt$p.value)) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
        plot.title = element_text(size = 11.5, hjust = 0))

ggsave(file.path(OUT_DIR, paste0(gsub("[^A-Za-z0-9]+", "_", FOCUS),
                                 "_participant_slopes.svg")),
       p_slopes, device = "svg", width = 4.4, height = 4.8)

cat("\nwrote figures to", OUT_DIR, "\n")



rc1 <- as.matrix(vegan::decostand(feat, method = "rclr", MARGIN = 1))
rc1[feat == 0] <- NA
implied1 <- log(feat[1, feat[1,] > 0]) - rc1[1, feat[1,] > 0]
summary(implied1)
