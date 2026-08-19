### Supplmenetary Table 1&2

library(dplyr)
library(tibble)
library(ggplot2)
library(gtsummary)

cat("vegan installed:", as.character(packageVersion("vegan")), "\n")

setwd("C:/Users/sirju/OneDrive/Desktop/AgingPaper2/For_github/Demographic_table_20260818")
df <- read.csv("0-other-res-df-group.csv", check.names = TRUE)

df$host_age             <- as.numeric(df$host_age)
df$host_body_mass_index <- as.numeric(df$host_body_mass_index)
df$physresilience       <- as.numeric(df$physresilience)

sex_vec <- factor(df$sex, levels = c("male", "female"))
grp_vec <- factor(df$physresilience_group, levels = c("Low", "High"))
fem <- sex_vec == "female"
hi  <- grp_vec == "High"


tab_df <- df %>%
  mutate(
    physresilience_group = factor(physresilience_group, levels = c("High", "Low")),
    sex       = factor(sex, levels = c("female", "male"),
                       labels = c("Female", "Male")),
    bmi_women = if_else(sex == "Female", host_body_mass_index, NA_real_),
    bmi_men   = if_else(sex == "Male",   host_body_mass_index, NA_real_)
  ) %>%
  select(physresilience_group, host_age, sex, host_body_mass_index,
         bmi_women, bmi_men, physresilience)

tblS1 <- tab_df %>%
  tbl_summary(
    by = physresilience_group,
    missing = "no",
    statistic = list(all_continuous()  ~ "{mean} ({sd})",
                     all_categorical() ~ "{n} ({p}%)"),
    digits = list(host_age ~ 2, host_body_mass_index ~ 2,
                  bmi_women ~ 2, bmi_men ~ 2, physresilience ~ 3),
    label = list(host_age             ~ "Age, years",
                 sex                  ~ "Sex",
                 host_body_mass_index ~ "BMI, kg/m\u00b2",
                 bmi_women            ~ "BMI, women only",
                 bmi_men              ~ "BMI, men only",
                 physresilience       ~ "Physical resilience score")
  ) %>%
  add_overall(col_label = "**Overall**  \nN = {N}") %>%
  add_p(test = list(all_continuous()  ~ "wilcox.test",
                    all_categorical() ~ "chisq.test")) %>%
  modify_table_styling(columns = label,
                       rows = variable %in% c("bmi_women", "bmi_men"),
                       text_format = "indent") %>%
  modify_footnote(
    p.value ~ paste(
      "Wilcoxon rank-sum test (age, BMI, physical resilience score);",
      "Pearson's Chi-squared test (sex). Indented BMI rows compare",
      "resilience groups within women and within men separately."
    )
  ) %>%
  bold_labels()

tblS1
tblS1 %>% as_gt() %>% gt::gtsave("TableS1_demographics.html")


rclr_manual <- function(m) {
  m  <- as.matrix(m)
  nz <- m > 0
  L  <- matrix(0, nrow(m), ncol(m), dimnames = dimnames(m))
  L[nz] <- log(m[nz])
  gm  <- rowSums(L) / rowSums(nz)
  out <- matrix(0, nrow(m), ncol(m), dimnames = dimnames(m))
  out[nz] <- (L - gm)[nz]
  out
}

metab_cols <- grep("^X[0-9]+$", colnames(df), value = TRUE)
cat("metabolite columns:", length(metab_cols), "\n")

df_rclr <- rclr_manual(df[, metab_cols])
df_rclr <- df_rclr[, apply(df_rclr, 2, var) > 1e-6]
cat("after variance filter:", ncol(df_rclr), "\n")

prop   <- 0.10
keep_n <- ceiling(prop * ncol(df_rclr))
X_sub  <- df_rclr[, order(apply(df_rclr, 2, var), decreasing = TRUE)[1:keep_n]]
Y      <- df$physresilience

cat(sprintf("retained for testing (top %.0f%% by variance): %d\n", prop * 100, keep_n))



for (s in c("female", "male")) {
  d  <- df[df$sex == s, ]
  wb <- wilcox.test(host_body_mass_index ~ physresilience_group, data = d)
  cat(sprintf("%-7s n = %3d | BMI High vs Low p = %.3f\n",
              s, nrow(d), wb$p.value))
}


cliff <- function(a, b) {
  r  <- rank(c(a, b))
  n1 <- length(a); n2 <- length(b)
  U  <- sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2
  2 * U / (n1 * n2) - 1
}

van_elteren <- function(x, g, strata) {
  num <- 0; vr <- 0
  for (s in unique(strata)) {
    m  <- strata == s
    xs <- x[m]; gs <- g[m]
    nh <- length(xs); n1 <- sum(gs)
    if (n1 == 0 || n1 == nh) next
    r <- rank(xs); w <- 1 / (nh + 1)
    num <- num + w * (sum(r[gs]) - n1 * (nh + 1) / 2)
    ties <- sum(sapply(table(xs), function(k) k^3 - k))
    vr <- vr + w^2 * n1 * (nh - n1) / (12 * nh * (nh - 1)) *
      ((nh^2 - 1) * nh - ties)
  }
  if (vr <= 0) return(1)
  2 * pnorm(-abs(num / sqrt(vr)))
}



p_unadj <- apply(X_sub, 2, function(x)
  suppressWarnings(wilcox.test(x[hi], x[!hi])$p.value))
q_unadj <- p.adjust(p_unadj, method = "BH")
sig <- which(q_unadj < 0.05)

cat(sprintf("\nSignificant at BH q < 0.05: %d of %d\n", length(sig), ncol(X_sub)))
if (length(sig) > 1000)
  warning("More than 1000 hits - check that rclr did not impute the zeros")

res <- lapply(sig, function(j) {
  x <- X_sub[, j]
  tibble(
    feature             = sub("^X", "", colnames(X_sub)[j]),
    cliff_delta_all     = cliff(x[hi], x[!hi]),
    cliff_delta_female  = cliff(x[fem & hi],  x[fem & !hi]),
    cliff_delta_male    = cliff(x[!fem & hi], x[!fem & !hi]),
    cliff_delta_sex     = cliff(x[fem], x[!fem]),   # size of the sex effect itself
    q_unadjusted        = q_unadj[j],
    p_stratified        = van_elteren(x, hi, as.integer(fem))
  )
}) %>%
  bind_rows() %>%
  mutate(q_stratified   = p.adjust(p_stratified, method = "BH"),
         same_direction = sign(cliff_delta_female) == sign(cliff_delta_male))

ct <- cor.test(res$cliff_delta_female, res$cliff_delta_male)

cat(sprintf("Cliff's delta, women vs men: r = %.3f (p = %.2e)\n",
            ct$estimate, ct$p.value))
cat(sprintf("Same direction in both sexes: %d of %d (%.0f%%)\n",
            sum(res$same_direction), nrow(res), 100 * mean(res$same_direction)))
cat(sprintf("Still significant after stratification: %d of %d (%.0f%%)\n",
            sum(res$q_stratified < 0.05), nrow(res), 100 * mean(res$q_stratified < 0.05)))
cat(sprintf("Median |delta|  overall %.3f | women %.3f | men %.3f\n",
            median(abs(res$cliff_delta_all)), median(abs(res$cliff_delta_female)),
            median(abs(res$cliff_delta_male))))

write.csv(res, "TableS2_sex_stratified.csv", row.names = FALSE)



lim <- max(abs(c(res$cliff_delta_female, res$cliff_delta_male))) * 1.15

figS1 <- ggplot(res, aes(cliff_delta_female, cliff_delta_male)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70", linewidth = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey85") +
  geom_smooth(method = "lm", formula = y ~ x,
              colour = "#0072B2", fill = "#0072B2", alpha = 0.18, linewidth = 1) +
  geom_point(aes(fill = cliff_delta_all > 0), size = 2.6, alpha = 0.75,
             shape = 21, colour = "white", stroke = 0.5) +
  scale_fill_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#D6A000"),
                    labels = c("TRUE"  = "Higher in high resilience",
                               "FALSE" = "Lower in high resilience"),
                    name = NULL) +
  annotate("label", x = -lim, y = lim, hjust = 0, vjust = 1, size = 5,
           fill = "white", parse = TRUE,
           label = sprintf("italic(r) == %.3f", ct$estimate)) +
  coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  labs(x = expression("Cliff's " * delta * " in women (" * italic(n) * " = 152)"),
       y = expression("Cliff's " * delta * " in men (" * italic(n) * " = 85)")) +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.line   = element_line(colour = "black"),
        axis.ticks  = element_line(colour = "black"),
        axis.ticks.length = unit(3, "pt"),
        legend.position = c(0.98, 0.02),
        legend.justification = c(1, 0),
        legend.background = element_blank())

figS1
ggsave("sex_sensitivity.svg", figS1, width = 6.5, height = 6.2, dpi = 300)
