#Figure S2b
library(tidyverse)
library(ggplot2)

setwd("...")

df_all   <- read.csv("0-ranked_metabolites-physres.csv")
df_noAHT <- read.csv("0-ranked_metabolites-physres-no-antihypertensives.csv")


m_all   <- df_all$Metabolite
m_noAHT <- df_noAHT$Metabolite

overlap  <- intersect(m_all, m_noAHT)
only_all <- setdiff(m_all, m_noAHT)
only_no  <- setdiff(m_noAHT, m_all)

cat("=== Overall Overlap ===\n")
cat("File1 (all patients)        :", length(m_all), "\n")
cat("File2 (no antihypertensives):", length(m_noAHT), "\n")
cat("Overlap                     :", length(overlap),
    sprintf("(%.1f%%)\n", length(overlap)/length(m_all)*100))
cat("Only in all-patients        :", length(only_all), "\n")
cat("Only in no-antihypert.      :", length(only_no), "\n")



df_merged <- inner_join(
  df_all   %>% select(Metabolite, Weight, Direction, Stability),
  df_noAHT %>% select(Metabolite, Weight, Direction, Stability),
  by = "Metabolite", suffix = c("_all", "_noAHT")
)

pearson <- cor.test(df_merged$Weight_all, df_merged$Weight_noAHT, method = "pearson")

cat("\n=== Weight Correlation (n =", nrow(df_merged), "overlapping metabolites) ===\n")
cat(sprintf("Pearson r = %.4f, p %s\n", pearson$estimate,
            ifelse(pearson$p.value < 2.2e-16, "< 2.2e-16", sprintf("= %.2e", pearson$p.value))))

dir_match <- mean(df_merged$Direction_all == df_merged$Direction_noAHT)
cat(sprintf("Direction match: %d/%d = %.2f%%\n",
            sum(df_merged$Direction_all == df_merged$Direction_noAHT),
            nrow(df_merged), dir_match * 100))


# Pearson r by Weight


cat("\n=== Correlation by |Weight| bin ===\n")
bins   <- c(0, 0.01, 0.02, 0.03, 0.04, Inf)
labels <- c("<0.01", "0.01-0.02", "0.02-0.03", "0.03-0.04", ">0.04")

df_merged <- df_merged %>% mutate(abs_w_all = abs(Weight_all))

for (i in seq_along(labels)) {
  sub <- df_merged %>%
    filter(abs_w_all >= bins[i], abs_w_all < bins[i+1])
  if (nrow(sub) > 5) {
    r <- cor(sub$Weight_all, sub$Weight_noAHT)
    cat(sprintf("  %s: n=%4d, r=%.4f\n", labels[i], nrow(sub), r))
  }
}



df_all_only   <- df_all   %>% filter(Metabolite %in% only_all)
df_noAHT_only <- df_noAHT %>% filter(Metabolite %in% only_no)
df_all_ov     <- df_all   %>% filter(Metabolite %in% overlap)
df_noAHT_ov   <- df_noAHT %>% filter(Metabolite %in% overlap)

cat("\n=== Overlap Characterization ===\n")
cat(sprintf("Overlapping     : %d / %d (%.1f%%)\n",
            length(overlap), length(m_all), length(overlap)/length(m_all)*100))
cat(sprintf("Non-overlapping : %d / %d (%.1f%%)\n\n",
            length(only_all), length(m_all), length(only_all)/length(m_all)*100))

cat("--- Median |Weight| ---\n")
cat(sprintf("  Overlapping    (all)    : %.5f\n", median(abs(df_all_ov$Weight))))
cat(sprintf("  Non-overlapping (all)   : %.5f\n", median(abs(df_all_only$Weight))))
cat(sprintf("  Overlapping    (noAHT)  : %.5f\n", median(abs(df_noAHT_ov$Weight))))
cat(sprintf("  Non-overlapping (noAHT) : %.5f\n\n", median(abs(df_noAHT_only$Weight))))

cat("--- Mean Stability ---\n")
cat(sprintf("  Overlapping    (all)    : %.4f\n", mean(df_all_ov$Stability)))
cat(sprintf("  Non-overlapping (all)   : %.4f\n", mean(df_all_only$Stability)))
cat(sprintf("  Overlapping    (noAHT)  : %.4f\n", mean(df_noAHT_ov$Stability)))
cat(sprintf("  Non-overlapping (noAHT) : %.4f\n", mean(df_noAHT_only$Stability)))


# Each point = one metabolite; axes = weight from each analysis
# Dashed line = line of identity (y = x); tight clustering indicates concordance

lim <- max(abs(c(df_merged$Weight_all, df_merged$Weight_noAHT))) * 1.05

p <- ggplot(df_merged, aes(x = Weight_all, y = Weight_noAHT)) +
  geom_point(alpha = 0.25, size = 1.2, color = "#3266ad") +
  geom_abline(intercept = 0, slope = 1,
              linetype = "dashed", color = "#D85A30", linewidth = 0.8) +
  annotate("text", x = -lim * 0.95, y = lim * 0.95,
           label = sprintf("r = %.3f", pearson$estimate),
           hjust = 0, vjust = 1, size = 3.5, color = "black") +
  scale_x_continuous(limits = c(-lim, lim)) +
  scale_y_continuous(limits = c(-lim, lim)) +
  labs(
    x = "Weight (all)",
    y = "Weight (excluding antihypertensive users)") +
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(size = 15, face = "bold"),
    plot.subtitle = element_text(size = 12, color = "grey40"),
    axis.title    = element_text(size = 13)
  ) +
  coord_fixed()

ggsave("sensitivity_scatter.svg", plot = p,
       width = 5, height = 5)

print(p)

cat("\nFigure saved: sensitivity_scatter.svg\n")
