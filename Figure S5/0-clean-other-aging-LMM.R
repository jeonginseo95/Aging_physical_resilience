#Packages
library(data.table)
library(dplyr)
library(tibble)
library(ggpubr)
library(stringr)
library(tidyr)
library(caret)
library(mixOmics)
library(tidyverse)
library(vegan)
library(ggrepel)
library(ggplot2)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(MuMIn)
library(performance)

#Loading data
df1 <- readRDS("clean_data_wellcome.rds")

metabs <- fread("0-ranked_metabolites-physres.csv")

carnitines <- c(
  7245, 9257, 14167, 28018, 33484,
  46322, 46997, 50582, 66650, 67549,
  71700, 75141, 77192, 78113, 80107,
  80586, 80681
)

glutamine_conjugates <- c(9454L,
                            12439L,
                            20833L,
                            40982L,
                            42850L,
                            43949L,
                            74560L)


phosphocholines <- c(
  89371L,
  100405L,
  102880L,
  103401L,
  104087L,
  106172L,
  107906L,
  136475L,
  141286L
)

BA <- c(62124,
        71371)

all_targets <- c(carnitines, glutamine_conjugates, phosphocholines, BA)

selected_metabs <- metabs[metabs$Metabolite %in% all_targets, ]

df1$host_age<-as.numeric(df1$host_age)
df1$host_body_mass_index<-as.numeric(df1$host_body_mass_index)
df1$physresilience_pctl<-as.numeric(df1$physresilience_pctl)
df1$physresilience<-as.numeric(df1$physresilience)

df <- df1[, colnames(df1) %in% c("filename", selected_metabs$Metabolite)]

#Selecting metabolites
df_metab <- df %>%
  column_to_rownames(var = "filename") %>%  
  dplyr::select(1:35) %>%                
  filter(rowSums(.) > 0)       


#Selecting prevelance
prev <- colMeans(df_metab > 0, na.rm = TRUE)

hist(
  prev,
  breaks = 50,
  col = "grey",
  border = "white",
  main = "Metabolite prevalence (non-zero proportion)",
  xlab = "Proportion of samples with non-zero values"
)

plot(
  sort(prev),
  seq_along(prev) / length(prev),
  type = "l",
  xlab = "Prevalence threshold",
  ylab = "Proportion of metabolites retained",
  main = "Cumulative metabolite prevalence"
) +
  abline(v = c(0.5, 0.6, 0.7), col = c("blue", "red", "darkgreen"), lty = 2)

sum(prev >= 0.7)
sum(prev >= 0.6)
sum(prev >= 0.5)
sum(prev >= 0.4)
sum(prev >= 0.3)
sum(prev >= 0.2)
sum(prev >= 0.1)
sum(prev >= 0.01)

#Cut-off
keep_metabs <- names(prev[prev >= 0.7])
df_filt <- df_metab[, keep_metabs]

#Pseudocount for downstream stats
pseudocount <- apply(df_filt, 2, function(x) min(x[x > 0], na.rm = TRUE) / 2)
df_pc <- sweep(df_filt, 2, pseudocount, "+")

#Log transform
df_log <- log2(df_pc)

#Prepping metadata and calculating within-person age
meta <- df1 %>%
  dplyr::select(filename, record_ID, visit, host_age, sex)

meta <- meta %>%
  group_by(record_ID) %>%
  mutate(age_within = host_age - mean(host_age)) %>%
  ungroup()

##Changing variable names to fit LMM package
# Save original names 
original_names <- colnames(df_log)

# Make valid R names
colnames(df_log) <- make.names(colnames(df_log))

# Create a mapping
name_map <- data.frame(original = original_names, valid = colnames(df_log))

#Combining with metadata
df_lmm <- cbind(meta[, c("record_ID", "visit", "age_within")], df_log)
df_lmm <- as.data.frame(df_lmm)

#LMM loop
results <- lapply(colnames(df_log), function(metab){
  formula <- as.formula(paste(metab, "~ age_within + (1|record_ID)"))
  fit <- lmer(formula, data = df_lmm)
  tidy(fit) %>%
    filter(term == "age_within") %>%
    mutate(metabolite = metab)  # <-- add metabolite name
})

results_df <- do.call(rbind, results)
results_df$p_adj <- p.adjust(results_df$p.value, method = "BH")

#Marking singular metabolites
results_df$singular <- sapply(results_df$metabolite, function(m){
  fit <- lmer(as.formula(paste(m, "~ age_within + (1|record_ID)")), data=df_lmm)
  isSingular(fit)
})
table(results_df$singular)


#Metabolite selection
final_metabs <- results_df %>%
  filter(
    #p_adj < 0.05,             # FDR threshold
    singular == FALSE          # non-singular random effect
  ) %>%
  arrange(desc(abs(estimate)))  

#Prepping plot df
final_metabs_clean <- final_metabs %>%
  mutate(metabolite = sub("^X", "", metabolite))

final_metabs_clean <- final_metabs_clean %>%
  rename(Metabolite = metabolite)

metabs$Metabolite <- as.character(metabs$Metabolite)

df_merged <- metabs %>%
  left_join(final_metabs_clean, by = "Metabolite")

#Stats
model <- lm(estimate ~ Weight, data = df_merged)
summary(model)
summary_model <- summary(model)


#stats
beta <- round(summary_model$coefficients["Weight", "Estimate"], 2)
p_val <- signif(summary_model$coefficients["Weight", "Pr(>|t|)"], 3)
r2 <- round(summary_model$r.squared, 2)

#position of stats annotation
x_pos <- max(df_merged$Weight, na.rm = TRUE) -
  0.05 * diff(range(df_merged$Weight, na.rm = TRUE))

y_pos <- min(df_merged$estimate, na.rm = TRUE) +
  0.25 * diff(range(df_merged$estimate, na.rm = TRUE))



table <- df_merged %>%
  dplyr::select(Metabolite, Weight, Direction, Stability, estimate, p_adj) %>%
  dplyr::filter(!is.na(estimate))

#write.csv(table, "0-other-aging-supl_table_2.csv", row.names = FALSE)

df1$visit <- factor(df1$visit,
                    levels = c("V5", "V8", "V10", "V12"))

#Plotting individual features

metab_age_33484 <- ggplot(df1, aes(x = host_age, y = `33484`, color = visit)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "#D95F0E", fill = "#FEE6CE") +
  scale_y_log10(labels = scales::label_log()) +
  theme_classic(base_size = 13) +
  theme(
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  ) +
  labs(
    x = "Host age",
    y = "Peak area (log10 scale)",
    color = "Visit"
  )


metab_age_77192 <- ggplot(df1, aes(x = host_age, y = `77192`, color = visit)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "#D95F0E", fill = "#FEE6CE") +
  scale_y_log10(labels = scales::label_log()) +
  theme_classic(base_size = 13) +
  theme(
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  ) +
  labs(
    x = "Host age",
    y = "Peak area (log10 scale)",
    color = "Visit"
  )

metab_age_77192
metab_age_33484

#ggsave("metab_age_77192.png", plot = metab_age_77192, width = 6, height = 4, dpi = 300)
#ggsave("metab_age_77192.svg", plot = metab_age_77192, width = 6, height = 4, dpi = 300)

#ggsave("metab_age_33484.png", plot = metab_age_33484, width = 6, height = 4, dpi = 300)
#ggsave("metab_age_33484.svg", plot = metab_age_33484, width = 6, height = 4, dpi = 300)

