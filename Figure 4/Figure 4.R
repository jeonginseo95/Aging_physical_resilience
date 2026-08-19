#Figure 4
#Assign each resilience according as four quartiles.
library(dplyr)
library(readr)
library(tidyverse)
library(ggplot2)
library(scales)
library(purrr)
library(broom)
library(ggpubr)
library(data.table)
library(pheatmap)
library(svglite)
library(conflicted)
library(effsize)
library(ggrepel)
conflicted::conflicts_prefer(dplyr::summarise)
conflicted::conflicts_prefer(dplyr::group_by)
conflicted::conflicts_prefer(dplyr::mutate)
conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::select)
conflicted::conflicts_prefer(dplyr::rename)

df <- read.csv("C:/Users/sirju/OneDrive/Desktop/AgingPaper2/For_github/Figure 4_20260811/0-other-res-df.csv")

df <- df %>%
  mutate(
    physresilience_quartile = ntile(physresilience, 4)
  )

write.csv(df, "C:/Users/sirju/OneDrive/Desktop/AgingPaper2/For_github/Figure 4_20260811/0-other-res-df-quartile.csv", row.names = FALSE)

df <- read.csv("C:/Users/sirju/OneDrive/Desktop/AgingPaper2/For_github/Figure 4_20260811/0-other-res-df-quartile.csv")

df <- df %>%
  mutate(
    physresilience_group = ifelse(physresilience_quartile %in% c(1, 2), "Low", "High")
  )

write.csv(df, "C:/Users/sirju/OneDrive/Desktop/AgingPaper2/For_github/Figure 4_20260811/0-other-res-df-group.csv", row.names = FALSE)



# Merging with library_compound_name column

# Set file path
base_path <- "C:/Users/sirju/OneDrive/Desktop/AgingPaper2/For_github/Figure 4_20260811/"

# Load files
metabolites <- read_csv(file.path(base_path, "0-ranked_metabolites-physres.csv"))
gnps <- read_tsv(file.path(base_path, "Library_results_with_gnps.tsv"))

# Match Metabolite with #Scan# and add Compound_Name and LibraryName as new columns
result <- metabolites %>%
  left_join(
    gnps %>% select(`#Scan#`, Compound_Name, LibraryName),
    by = c("Metabolite" = "#Scan#")
  )

# Save result
write_csv(result, file.path(base_path, "0-ranked_metabolites-physres_with_compound.csv"))

cat("Done! Matched rows:", sum(!is.na(result$Compound_Name)), "\n")
cat("Total rows:", nrow(result), "\n")


# Merging with classification and BA query_validation columns

# Set file path
base_path <- "C:/Users/sirju/OneDrive/Desktop/AgingPaper2/For_github/Figure 4_20260811/"

# Load files
metabolites <- read_csv(file.path(base_path, "0-ranked_metabolites-physres_with_compound.csv"))
ba <- read_csv(file.path(base_path, "BA_annotation_results.csv"))

# Match Metabolite with #Scan# and add classification, query_validation
result <- metabolites %>%
  left_join(
    ba %>% select(`#Scan#`, classification, query_validation),
    by = c("Metabolite" = "#Scan#")
  )

# Save result (overwrite the same file)
write_csv(result, file.path(base_path, "0-ranked_metabolites-physres_with_compound.csv"))

cat("Done! Matched rows (classification):", sum(!is.na(result$classification)), "\n")
cat("Done! Matched rows (query_validation):", sum(!is.na(result$query_validation)), "\n")
cat("Total rows:", nrow(result), "\n")




#Figure 4b
### Wilcoxon + Bar chart: Carnitine (filter: adjusted p < 0.05; FC > 1)

# Load data
base_path <- "C:/Users/sirju/OneDrive/Desktop/AgingPaper2/For_github/Figure 4_20260811/"

main_df   <- read.csv(file.path(base_path, "0-other-res-df-group.csv"))
phys_meta <- read.csv(file.path(base_path, "0-ranked_metabolites-physres_with_compound.csv"))

# Filter carnitine features
carnitine_meta <- phys_meta %>%
  filter(grepl("carnitine", Compound_Name, ignore.case = TRUE))

cat("Number of carnitine features:", nrow(carnitine_meta), "\n")

carnitine_ids <- carnitine_meta$Metabolite

# Row-sum normalization on ALL features
all_feature_cols <- grep("^X", colnames(main_df), value = TRUE)

main_norm <- main_df %>%
  mutate(row_sum = rowSums(select(., all_of(all_feature_cols)))) %>%
  mutate(across(all_of(all_feature_cols), ~ . / row_sum)) %>%
  select(-row_sum)

# Subset carnitine features only
carnitine_cols <- paste0("X", carnitine_ids)
carnitine_cols <- intersect(carnitine_cols, colnames(main_norm))

cat("Matched carnitine columns in main_df:", length(carnitine_cols), "\n")

phys <- main_norm %>%
  filter(physresilience_group %in% c("High", "Low")) %>%
  select(physresilience_group, any_of(carnitine_cols)) %>%
  rename_with(~ sub("^X", "", .), starts_with("X"))

# Wilcoxon rank-sum test + BH correction + Log2FC
results_pval <- phys %>%
  pivot_longer(-physresilience_group, names_to = "Feature", values_to = "Value") %>%
  group_by(Feature) %>%
  summarise(
    p_value = tryCatch(
      wilcox.test(Value ~ physresilience_group)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(adj_p_value = p.adjust(p_value, method = "BH"))

# Log2FC
results_log2fc <- phys %>%
  pivot_longer(-physresilience_group, names_to = "Feature", values_to = "Value") %>%
  group_by(physresilience_group, Feature) %>%
  summarise(mean_val = mean(Value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = physresilience_group, values_from = mean_val) %>%
  mutate(across(c("High", "Low"), ~ ifelse(. == 0, 1e-8, .))) %>%
  mutate(Log2FC = log2(High / Low))

results <- results_log2fc %>%
  left_join(results_pval, by = "Feature")

cat("Significant features (adj_p < 0.05):", sum(results$adj_p_value < 0.05, na.rm = TRUE), "\n")

# Filter significant & add Compound_Name
label_map <- carnitine_meta %>%
  select(Metabolite, Compound_Name) %>%
  mutate(Metabolite = as.character(Metabolite))

sig_results <- results %>%
  filter(adj_p_value < 0.05 & abs(Log2FC) > 1) %>%
  mutate(Feature = as.character(Feature)) %>%
  left_join(label_map, by = c("Feature" = "Metabolite")) %>%
  group_by(Compound_Name) %>%
  mutate(label = if_else(n() > 1,
                         paste0(Compound_Name, " [", Feature, "]"),
                         Compound_Name)) %>%
  ungroup() %>%
  arrange(Log2FC)


# Cliff's delta for significant carnitine features
cliff_results_carnitine <- purrr::map_dfr(unique(sig_results$Feature), function(feat) {
  vals <- phys %>%
    select(physresilience_group, all_of(feat)) %>%
    rename(Value = all_of(feat)) %>%
    mutate(physresilience_group = factor(physresilience_group, levels = c("High", "Low")))
  
  cd <- tryCatch(
    effsize::cliff.delta(Value ~ physresilience_group, data = vals),
    error = function(e) NULL
  )
  
  if (is.null(cd)) {
    tibble(Feature = feat, cliffs_delta = NA_real_, magnitude = NA_character_,
           ci_lower = NA_real_, ci_upper = NA_real_)
  } else {
    tibble(
      Feature      = feat,
      cliffs_delta = as.numeric(cd$estimate),
      magnitude    = as.character(cd$magnitude),
      ci_lower     = cd$conf.int[1],
      ci_upper     = cd$conf.int[2]
    )
  }
})

sig_results_carnitine_cliff <- sig_results %>%
  mutate(Feature = as.character(Feature)) %>%
  left_join(cliff_results_carnitine, by = "Feature")

write_csv(sig_results_carnitine_cliff, file.path(base_path, "cliffsdelta_phys_carnitine.csv"))
cat("Saved: cliffsdelta_phys_carnitine.csv\n")


# Bar chart
if (nrow(sig_results) == 0) {
  cat("No significant features (adj_p < 0.05) to plot.\n")
} else {
  p <- ggplot(sig_results, aes(x = reorder(label, Log2FC), y = Log2FC, fill = Log2FC > 0)) +
    geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.6) +
    scale_fill_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#D09200")) +    coord_flip() +
    labs(
      title = "Wilcoxon Rank-Sum Test: Physical Resilience\n(Carnitine Features, BH adj_p < 0.05)",
      x     = NULL,
      y     = "Log2FC (High / Low)"
    ) +
    theme_minimal() +
    theme(
      plot.title        = element_text(hjust = 0.5, size = 13, face = "bold"),
      axis.title        = element_text(size = 11),
      axis.text         = element_text(size = 9),
      axis.text.y       = element_text(margin = margin(r = 5)),
      panel.grid        = element_blank(),
      axis.line         = element_line(color = "black"),
      axis.ticks        = element_line(color = "black"),
      axis.ticks.length = unit(0.15, "cm"),
      legend.position   = "none",
      plot.margin       = margin(10, 10, 10, 120)
    )
  
  print(p)
  ggsave(file.path(base_path, "barplot_phys_carnitine_wilcoxon.svg"),
         plot = p, device = "svg", width = 20, height = max(3, nrow(sig_results) * 1.2), dpi = "retina")
  cat("Saved: barplot_phys_carnitine_wilcoxon.svg\n")
}





###Figure 4a
library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)

base_path <- "C:/Users/sirju/OneDrive/Desktop/AgingPaper2/For_github/Figure 4_20260811/"


df <- read_csv(file.path(base_path, "results_phys_carnitine_all_with_flag.csv"))

df <- df %>%
  mutate(
    chain_length = str_match(Compound_Name, "putative explanation:\\s*CAR\\s*C(\\d+)")[, 2],
    chain_length = as.numeric(chain_length),
    chain_length = if_else(
      is.na(chain_length) & str_detect(Compound_Name, regex("undecanoic", ignore_case = TRUE)),
      11, chain_length
    ),
    chain_category = case_when(
      is.na(chain_length)                     ~ "No putative explanation",
      chain_length >= 2  & chain_length <= 5   ~ "SHORT (C2-C5)",
      chain_length >= 6  & chain_length <= 12  ~ "MEDIUM (C6-C12)",
      chain_length >= 13 & chain_length <= 21  ~ "LONG (C13-C21)",
      chain_length >= 22                       ~ "VERY LONG (>=C22)",
      TRUE ~ "No putative explanation"
    ),
    chain_category = factor(chain_category,
                            levels = c("No putative explanation", "SHORT (C2-C5)",
                                       "MEDIUM (C6-C12)", "LONG (C13-C21)", "VERY LONG (>=C22)")
    ),
    neglog10p = -log10(adj_p_value),   # <- adjusted p-value, matches `significant`
    car_label = str_match(Compound_Name, "putative explanation:\\s*CAR\\s*([^)]+)")[, 2],
    car_label = str_trim(car_label),
    delta_mass = str_match(Compound_Name, "delta mass\\s*([0-9]+\\.[0-9]+)")[, 2],
    car_label = if_else(is.na(car_label) & !is.na(delta_mass),
                        paste0("\u0394m/z ", delta_mass), car_label),
    car_label = if_else(
      is.na(car_label) & str_detect(Compound_Name, regex("undecanoic", ignore_case = TRUE)),
      "C11:0", car_label
    ),
    label_this = significant & abs(Log2FC) > 1
  )


stopifnot(all(df$significant == (df$adj_p_value < 0.05)))
cat(sprintf("significant rows: %d (all adj_p_value < 0.05: %s)\n",
            sum(df$significant), all(df$adj_p_value[df$significant] < 0.05)))

labelled <- df %>% filter(label_this)
cat(sprintf("rows to be labelled: %d (%d with CAR name, %d with delta-mass fallback, %d unlabelled)\n",
            nrow(labelled),
            sum(str_detect(labelled$Compound_Name, "putative explanation:\\s*CAR")),
            sum(!str_detect(labelled$Compound_Name, "putative explanation:\\s*CAR") &
                  str_detect(labelled$Compound_Name, "delta mass")),
            sum(is.na(labelled$car_label))))

colnames(df)
stopifnot(all(c("neglog10p", "chain_category", "car_label", "label_this") %in% colnames(df)))


color_map <- c(
  "No putative explanation" = "grey70",
  "SHORT (C2-C5)"           = "#440154",
  "MEDIUM (C6-C12)"         = "#31688e",
  "LONG (C13-C21)"          = "#35b779",
  "VERY LONG (>=C22)"       = "#fde725"
)


p <- ggplot(df, aes(x = Log2FC, y = neglog10p, fill = chain_category)) +
  geom_point(shape = 21, color = "black", size = 2.5, alpha = 0.75, stroke = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "grey40") +
  geom_text_repel(
    data = df %>% filter(label_this & !is.na(car_label)),
    aes(label = car_label),
    size = 5, fontface = "bold", color = "black",
    segment.color = "black", segment.size = 0.4,
    min.segment.length = 0,
    max.overlaps = Inf, box.padding = 0.4,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = color_map, name = "Acyl chain length") +
  labs(
    title = "Volcano Plot: Carnitine Features by Acyl Chain Length",
    x = "Log2FC (High / Low)",
    y = "-log10(adjusted p-value)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title  = element_text(hjust = 0.5, face = "bold"),
    axis.title  = element_text(size = 16),   
    axis.text   = element_text(size = 14)    
  )

print(p)

ggsave(file.path(base_path, "volcano_carnitine_chainlength.svg"),
       plot = p, device = "svg", width = 10, height = 8, dpi = "retina")










#Figure 4c
### Wilcoxon + Bar chart: Glutamine (filter: adjusted p < 0.05)

# Load data
base_path <- "C:/Users/sirju/OneDrive/Desktop/AgingPaper2/For_github/Figure 4_20260811/"

main_df   <- read.csv(file.path(base_path, "0-other-res-df-group.csv"))
phys_meta <- read.csv(file.path(base_path, "0-ranked_metabolites-physres_with_compound.csv"))

# Filter glutamine features
glutamine_meta <- phys_meta %>%
  filter(grepl("glutamine", Compound_Name, ignore.case = TRUE))

cat("Number of glutamine features:", nrow(glutamine_meta), "\n")

glutamine_ids <- glutamine_meta$Metabolite

# Row-sum normalization on ALL features
all_feature_cols <- grep("^X", colnames(main_df), value = TRUE)

main_norm <- main_df %>%
  mutate(row_sum = rowSums(select(., all_of(all_feature_cols)))) %>%
  mutate(across(all_of(all_feature_cols), ~ . / row_sum)) %>%
  select(-row_sum)

# Subset glutamine features only
glutamine_cols <- paste0("X", glutamine_ids)
glutamine_cols <- intersect(glutamine_cols, colnames(main_norm))

cat("Matched glutamine columns in main_df:", length(glutamine_cols), "\n")

phys <- main_norm %>%
  filter(physresilience_group %in% c("High", "Low")) %>%
  select(physresilience_group, any_of(glutamine_cols)) %>%
  rename_with(~ sub("^X", "", .), starts_with("X"))

# Wilcoxon rank-sum test + BH correction + Log2FC
results_pval <- phys %>%
  pivot_longer(-physresilience_group, names_to = "Feature", values_to = "Value") %>%
  group_by(Feature) %>%
  summarise(
    p_value = tryCatch(
      wilcox.test(Value ~ physresilience_group)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(adj_p_value = p.adjust(p_value, method = "BH"))

# Log2FC
results_log2fc <- phys %>%
  pivot_longer(-physresilience_group, names_to = "Feature", values_to = "Value") %>%
  group_by(physresilience_group, Feature) %>%
  summarise(mean_val = mean(Value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = physresilience_group, values_from = mean_val) %>%
  mutate(across(c("High", "Low"), ~ ifelse(. == 0, 1e-8, .))) %>%
  mutate(Log2FC = log2(High / Low))

results <- results_log2fc %>%
  left_join(results_pval, by = "Feature")

cat("Significant features (adj_p < 0.05):", sum(results$adj_p_value < 0.05, na.rm = TRUE), "\n")

# Filter significant & add Compound_Name
label_map <- glutamine_meta %>%
  select(Metabolite, Compound_Name) %>%
  mutate(Metabolite = as.character(Metabolite))

sig_results <- results %>%
  filter(adj_p_value < 0.05) %>%
  mutate(Feature = as.character(Feature)) %>%
  left_join(label_map, by = c("Feature" = "Metabolite")) %>%
  group_by(Compound_Name) %>%
  mutate(label = paste0(Compound_Name, " [", Feature, "]")) %>%
  ungroup() %>%
  arrange(Log2FC)


# Cliff's delta for significant glutamine features
cliff_results_glutamine <- purrr::map_dfr(unique(sig_results$Feature), function(feat) {
  vals <- phys %>%
    select(physresilience_group, all_of(feat)) %>%
    rename(Value = all_of(feat)) %>%
    mutate(physresilience_group = factor(physresilience_group, levels = c("High", "Low")))
  
  cd <- tryCatch(
    effsize::cliff.delta(Value ~ physresilience_group, data = vals),
    error = function(e) NULL
  )
  
  if (is.null(cd)) {
    tibble(Feature = feat, cliffs_delta = NA_real_, magnitude = NA_character_,
           ci_lower = NA_real_, ci_upper = NA_real_)
  } else {
    tibble(
      Feature      = feat,
      cliffs_delta = as.numeric(cd$estimate),
      magnitude    = as.character(cd$magnitude),
      ci_lower     = cd$conf.int[1],
      ci_upper     = cd$conf.int[2]
    )
  }
})

sig_results_glutamine_cliff <- sig_results %>%
  mutate(Feature = as.character(Feature)) %>%
  left_join(cliff_results_glutamine, by = "Feature")

write_csv(sig_results_glutamine_cliff, file.path(base_path, "cliffsdelta_phys_glutamine.csv"))
cat("Saved: cliffsdelta_phys_glutamine.csv\n")



# Bar chart
if (nrow(sig_results) == 0) {
  cat("No significant features (adj_p < 0.05) to plot.\n")
} else {
  p <- ggplot(sig_results, aes(x = reorder(label, Log2FC), y = Log2FC, fill = Log2FC > 0)) +
    geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.6) +
    scale_fill_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#D09200")) +
    coord_flip() +
    labs(
      title = "Wilcoxon Rank-Sum Test: Physical Resilience\n(Glutamine Features, BH adj_p < 0.05)",
      x     = NULL,
      y     = "Log2FC (High / Low)"
    ) +
    theme_minimal() +
    theme(
      plot.title        = element_text(hjust = 0.5, size = 13, face = "bold"),
      axis.title        = element_text(size = 11),
      axis.text         = element_text(size = 9),
      axis.text.y       = element_text(margin = margin(r = 5)),
      panel.grid        = element_blank(),
      axis.line         = element_line(color = "black"),
      axis.ticks        = element_line(color = "black"),
      axis.ticks.length = unit(0.15, "cm"),
      legend.position   = "none",
      plot.margin       = margin(10, 10, 10, 120)
    )
  
  print(p)
  ggsave(file.path(base_path, "barplot_phys_glutamine_wilcoxon.svg"),
         plot = p, device = "svg", width = 12, height = max(3, nrow(sig_results) * 1.2), dpi = "retina")
  cat("Saved: barplot_phys_glutamine_wilcoxon.svg\n")
}





#Figure 4d
### Wilcoxon + Bar chart: Phosphocholine / PC-related (filter: adjusted p < 0.05)
# Load data
base_path <- "C:/Users/sirju/OneDrive/Desktop/AgingPaper2/For_github/Figure 4_20260811/"

main_df   <- read.csv(file.path(base_path, "0-other-res-df-group.csv"))
phys_meta <- read.csv(file.path(base_path, "0-ranked_metabolites-physres_with_compound.csv"))

# Filter phosphocholine-related features
pc_keywords <- c("phosphorylcholine", "phosphocholine", "Lyso PC", "PC\\(", "PAF", "PC-", "LPC", "SM\\(", "-PC")
pc_pattern  <- paste(pc_keywords, collapse = "|")

pc_meta <- phys_meta %>%
  filter(grepl(pc_pattern, Compound_Name, ignore.case = TRUE))

cat("Number of PC-related features:", nrow(pc_meta), "\n")

pc_ids <- pc_meta$Metabolite

# Row-sum normalization on ALL features
all_feature_cols <- grep("^X", colnames(main_df), value = TRUE)

main_norm <- main_df %>%
  mutate(row_sum = rowSums(select(., all_of(all_feature_cols)))) %>%
  mutate(across(all_of(all_feature_cols), ~ . / row_sum)) %>%
  select(-row_sum)

# Subset PC-related features only
pc_cols <- paste0("X", pc_ids)
pc_cols <- intersect(pc_cols, colnames(main_norm))

cat("Matched PC-related columns in main_df:", length(pc_cols), "\n")

phys <- main_norm %>%
  filter(physresilience_group %in% c("High", "Low")) %>%
  select(physresilience_group, any_of(pc_cols)) %>%
  rename_with(~ sub("^X", "", .), starts_with("X"))

# Wilcoxon rank-sum test + BH correction + Log2FC
results_pval <- phys %>%
  pivot_longer(-physresilience_group, names_to = "Feature", values_to = "Value") %>%
  group_by(Feature) %>%
  summarise(
    p_value = tryCatch(
      wilcox.test(Value ~ physresilience_group)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(adj_p_value = p.adjust(p_value, method = "BH"))

# Log2FC
results_log2fc <- phys %>%
  pivot_longer(-physresilience_group, names_to = "Feature", values_to = "Value") %>%
  group_by(physresilience_group, Feature) %>%
  summarise(mean_val = mean(Value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = physresilience_group, values_from = mean_val) %>%
  mutate(across(c("High", "Low"), ~ ifelse(. == 0, 1e-8, .))) %>%
  mutate(Log2FC = log2(High / Low))

results <- results_log2fc %>%
  left_join(results_pval, by = "Feature")

cat("Significant features (adj_p < 0.05):", sum(results$adj_p_value < 0.05, na.rm = TRUE), "\n")

# Filter significant & add Compound_Name
label_map <- pc_meta %>%
  select(Metabolite, Compound_Name) %>%
  mutate(Metabolite = as.character(Metabolite))

sig_results <- results %>%
  filter(adj_p_value < 0.05) %>%
  mutate(Feature = as.character(Feature)) %>%
  left_join(label_map, by = c("Feature" = "Metabolite")) %>%
  group_by(Compound_Name) %>%
  mutate(label = paste0(Compound_Name, " [", Feature, "]")) %>%
  ungroup() %>%
  arrange(Log2FC)


# Cliff's delta for significant PC-related features
cliff_results_pc <- purrr::map_dfr(unique(sig_results$Feature), function(feat) {
  vals <- phys %>%
    select(physresilience_group, all_of(feat)) %>%
    rename(Value = all_of(feat)) %>%
    mutate(physresilience_group = factor(physresilience_group, levels = c("High", "Low")))
  
  cd <- tryCatch(
    effsize::cliff.delta(Value ~ physresilience_group, data = vals),
    error = function(e) NULL
  )
  
  if (is.null(cd)) {
    tibble(Feature = feat, cliffs_delta = NA_real_, magnitude = NA_character_,
           ci_lower = NA_real_, ci_upper = NA_real_)
  } else {
    tibble(
      Feature      = feat,
      cliffs_delta = as.numeric(cd$estimate),
      magnitude    = as.character(cd$magnitude),
      ci_lower     = cd$conf.int[1],
      ci_upper     = cd$conf.int[2]
    )
  }
})

sig_results_pc_cliff <- sig_results %>%
  mutate(Feature = as.character(Feature)) %>%
  left_join(cliff_results_pc, by = "Feature")

write_csv(sig_results_pc_cliff, file.path(base_path, "cliffsdelta_phys_PC.csv"))
cat("Saved: cliffsdelta_phys_PC.csv\n")




# Bar chart
if (nrow(sig_results) == 0) {
  cat("No significant features (adj_p < 0.05) to plot.\n")
} else {
  p <- ggplot(sig_results, aes(x = reorder(label, Log2FC), y = Log2FC, fill = Log2FC > 0)) +
    geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.6) +
    scale_fill_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#D09200")) +
    coord_flip() +
    labs(
      title = "Wilcoxon Rank-Sum Test: Physical Resilience\n(PC-related Features, BH adj_p < 0.05)",
      x     = NULL,
      y     = "Log2FC (High / Low)"
    ) +
    theme_minimal() +
    theme(
      plot.title        = element_text(hjust = 0.5, size = 13, face = "bold"),
      axis.title        = element_text(size = 11),
      axis.text         = element_text(size = 9),
      axis.text.y       = element_text(margin = margin(r = 5)),
      panel.grid        = element_blank(),
      axis.line         = element_line(color = "black"),
      axis.ticks        = element_line(color = "black"),
      axis.ticks.length = unit(0.15, "cm"),
      legend.position   = "none",
      plot.margin       = margin(10, 10, 10, 120)
    )
  
  print(p)
  ggsave(file.path(base_path, "barplot_phys_PC_wilcoxon.svg"),
         plot = p, device = "svg", width = 12, height = max(3, nrow(sig_results) * 1.2), dpi = "retina")
  cat("Saved: barplot_phys_PC_wilcoxon.svg\n")
}






# Figure 4d-f
# Fisher's exact test
# Load data
df <- read.csv("C:/Users/sirju/OneDrive/Desktop/AgingPaper2/For_github/Figure 4_20260811/NumberOf_match_unmatch_20260412.csv",
               stringsAsFactors = FALSE)

features <- unique(df$Feature)

results_all <- data.frame()

for (feat in features) {
  cat("\n", strrep("=", 60), "\n")
  cat("Feature:", feat, "\n")
  cat(strrep("=", 60), "\n")
  
  df_feat <- df %>% filter(Feature == feat)
  healthy_row <- df_feat %>% filter(DOIDCommonName == "healthy")
  disease_rows <- df_feat %>% filter(DOIDCommonName != "healthy")
  
  if (nrow(healthy_row) == 0) {
    cat("No healthy row found, skipping\n")
    next
  }
  
  for (i in 1:nrow(disease_rows)) {
    disease <- disease_rows$DOIDCommonName[i]
    
    # Build 2x2 contingency table
    # Rows: Presence / Absence
    # Cols: Disease / Healthy
    mat <- matrix(
      c(disease_rows$Presence[i], disease_rows$Absence[i],
        healthy_row$Presence,     healthy_row$Absence),
      nrow = 2, byrow = FALSE,
      dimnames = list(
        c("Presence", "Absence"),
        c(disease, "healthy")
      )
    )
    
    cat("\n[", disease, "vs healthy ]\n")
    print(mat)
    
    # Fisher's Exact Test
    result <- fisher.test(mat)
    cat("p-value =", format(result$p.value, digits = 4),
        "| Odds Ratio =", round(result$estimate, 4), "\n")
    
    # Store result
    results_all <- rbind(results_all, data.frame(
      Feature         = feat,
      Disease         = disease,
      Disease_Presence = disease_rows$Presence[i],
      Disease_Absence  = disease_rows$Absence[i],
      Healthy_Presence = healthy_row$Presence,
      Healthy_Absence  = healthy_row$Absence,
      OddsRatio       = round(result$estimate, 4),
      P_value         = result$p.value
    ))
  }
}

# Save to CSV
output_path <- "C:/Users/sirju/OneDrive/Desktop/AgingPaper2/For_github/Figure 4_20260811/Fisher_results.csv"
write.csv(results_all, output_path, row.names = FALSE)
cat("\nResults saved to:", output_path, "\n")

