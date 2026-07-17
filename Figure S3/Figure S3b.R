#Figure S3b
# Quadrant Plot: Cognitive vs Physical Resilience Metabolites
library(dplyr)
library(ggplot2)
library(ggrepel)

# Load data
setwd("...")

cog <- read.csv("ranked_metabolites_cognitive_resilience_updated.csv",
                stringsAsFactors = FALSE)

phy <- read.csv("0-ranked_metabolites-physres_with_compound.csv",
                stringsAsFactors = FALSE)

# Filter by top 10% AbsWeight
cog_cutoff <- quantile(cog$AbsWeight,          0.90)
phy_cutoff <- quantile(phy$AbsWeight_physical, 0.90)

cat(sprintf("Cognitive cutoff (90th pct) : %.5f\n", cog_cutoff))
cat(sprintf("Physical  cutoff (90th pct) : %.5f\n", phy_cutoff))

cog_filtered <- cog %>%
  filter(AbsWeight >= cog_cutoff) %>%
  select(Metabolite,
         Weight_cog    = Weight,
         AbsWeight_cog = AbsWeight,
         Label_cog     = Compound_Name)

phy_filtered <- phy %>%
  filter(AbsWeight_physical >= phy_cutoff) %>%
  select(Metabolite    = Metabolite_physical,
         Weight_phy    = Weight_physical,
         AbsWeight_phy = AbsWeight_physical,
         Label_phy     = Compound_Name_phy)

merged <- full_join(
  cog_filtered %>% mutate(Metabolite = as.character(Metabolite)),
  phy_filtered  %>% mutate(Metabolite = as.character(Metabolite)),
  by = "Metabolite"
) %>%
  mutate(
    Weight_cog    = ifelse(is.na(Weight_cog),    0, Weight_cog),
    Weight_phy    = ifelse(is.na(Weight_phy),    0, Weight_phy),
    AbsWeight_cog = ifelse(is.na(AbsWeight_cog), 0, AbsWeight_cog),
    AbsWeight_phy = ifelse(is.na(AbsWeight_phy), 0, AbsWeight_phy)
  )

merged <- merged %>%
  mutate(
    Source = case_when(
      AbsWeight_cog > 0 & AbsWeight_phy > 0 ~ "Both",
      AbsWeight_cog > 0                      ~ "Cognitive only",
      TRUE                                   ~ "Physical only"
    )
  )

cat(sprintf("Total metabolites plotted : %d\n", nrow(merged)))
cat(sprintf("  Shared (both)           : %d\n", sum(merged$Source == "Both")))
cat(sprintf("  Cognitive only          : %d\n", sum(merged$Source == "Cognitive only")))
cat(sprintf("  Physical only           : %d\n", sum(merged$Source == "Physical only")))

# Assign quadrant labels
merged <- merged %>%
  mutate(
    Quadrant = case_when(
      Weight_cog >= 0 & Weight_phy >= 0 ~ "Q1: Phy+ / Cog+",
      Weight_cog <  0 & Weight_phy >= 0 ~ "Q2: Phy+ / Cog-",
      Weight_cog <  0 & Weight_phy <  0 ~ "Q3: Phy- / Cog-",
      Weight_cog >= 0 & Weight_phy <  0 ~ "Q4: Phy- / Cog+"
    )
  )

print(table(merged$Quadrant))

label_shared <- merged %>%
  filter(Source == "Both") %>%
  filter(!is.na(Label_cog) & Label_cog != "NA" & nchar(trimws(Label_cog)) > 0) %>%
  group_by(Quadrant) %>%
  mutate(AvgAbs = (AbsWeight_cog + AbsWeight_phy) / 2) %>%
  slice_max(AvgAbs, n = 5) %>%
  ungroup()

label_phy_only <- merged %>%
  filter(Source == "Physical only") %>%
  filter(!is.na(Label_phy) & Label_phy != "NA" & nchar(trimws(Label_phy)) > 0) %>%
  slice_max(AbsWeight_phy, n = 5)

label_df <- bind_rows(label_shared, label_phy_only)

# Assign colors
merged <- merged %>%
  mutate(
    PlotColor = case_when(
      Source == "Cognitive only"                  ~ "grey75",
      Source == "Both"                            ~ "#7B6BB5",
      Source == "Physical only" & Weight_phy >= 0 ~ "#3690C0",
      Source == "Physical only" & Weight_phy <  0 ~ "#D09200"
    )
  )

label_df <- label_df %>%
  mutate(
    PlotColor = case_when(
      Source == "Both"                            ~ "#7B6BB5",
      Source == "Physical only" & Weight_phy >= 0 ~ "#3690C0",
      Source == "Physical only" & Weight_phy <  0 ~ "#D09200"
    )
  )

merged_cog  <- merged %>% filter(Source == "Cognitive only")
merged_phy  <- merged %>% filter(Source == "Physical only")
merged_both <- merged %>% filter(Source == "Both")

size_min <- min(c(
  merged_cog$AbsWeight_cog,
  merged_phy$AbsWeight_phy,
  (merged_both$AbsWeight_cog + merged_both$AbsWeight_phy) / 2
))
size_max <- max(c(
  merged_cog$AbsWeight_cog,
  merged_phy$AbsWeight_phy,
  (merged_both$AbsWeight_cog + merged_both$AbsWeight_phy) / 2
))

# Build plot
p <- ggplot(merged, aes(x = Weight_cog, y = Weight_phy)) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  
  geom_point(data  = merged_cog,
             aes(size = AbsWeight_cog),
             color = "grey75", shape = 16, alpha = 0.5) +
  
  geom_point(data  = merged_phy,
             aes(size = AbsWeight_phy, color = PlotColor),
             shape = 16, alpha = 0.55) +
  
  geom_point(data  = merged_both,
             aes(size = (AbsWeight_cog + AbsWeight_phy) / 2, color = PlotColor),
             shape = 16, alpha = 0.60) +
  
  geom_text_repel(
    data               = label_df,
    aes(label          = Label_phy),
    color              = "black",
    size               = 2.4,
    max.overlaps       = 30,
    box.padding        = 0.6,
    point.padding      = 0,        # line goes into the dot
    min.segment.length = 0,        # always draw segment
    segment.color      = "black",
    segment.size       = 0.4,
    segment.linetype   = 1,
    arrow              = NULL,     # no arrow
    force              = 3,
    force_pull         = 0.5,
    show.legend        = FALSE
  ) +
  
  scale_color_identity(
    guide  = "legend",
    name   = "Group",
    breaks = c("#7B6BB5", "#3690C0", "#D09200", "grey75"),
    labels = c("Both", "Physical resilience only (Positive)",
               "Physical resilience only (Negative)", "Cognitive only")
  ) +
  
  scale_size_continuous(
    range  = c(1.2, 5),
    name   = "Weight",
    breaks = c(size_min, size_max),
    labels = c(sprintf("%.3f", size_min),
               sprintf("%.3f", size_max))
  ) +
  
  annotate("text", x =  Inf, y =  Inf,
           label = "Physical resilience +\nCognitive resilience +",
           hjust = 1.1, vjust = 1.4,  size = 4.0, color = "black", fontface = "bold") +
  annotate("text", x = -Inf, y =  Inf,
           label = "Physical resilience +\nCognitive resilience \u2212",
           hjust = -0.1, vjust = 1.4, size = 4.0, color = "black", fontface = "bold") +
  annotate("text", x = -Inf, y = -Inf,
           label = "Physical resilience \u2212\nCognitive resilience \u2212",
           hjust = -0.1, vjust = -0.5, size = 4.0, color = "black", fontface = "bold") +
  annotate("text", x =  Inf, y = -Inf,
           label = "Physical resilience \u2212\nCognitive resilience +",
           hjust = 1.1, vjust = -0.5, size = 4.0, color = "black", fontface = "bold") +
  
  labs(
    title    = "Cognitive vs Physical Resilience \u2014 Top 10% by Absolute Weight",
    subtitle = sprintf("Total: %d  |  Shared: %d  |  Cognitive only: %d  |  Physical only: %d",
                       nrow(merged),
                       sum(merged$Source == "Both"),
                       sum(merged$Source == "Cognitive only"),
                       sum(merged$Source == "Physical only")),
    x = "Cognitive resilience (Weight)",
    y = "Physical resilience (Weight)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(color = "grey40", size = 10),
    legend.position  = "right",
    axis.title.x     = element_text(size = 15, face = "bold"),
    axis.title.y     = element_text(size = 15, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Save and preview
ggsave("quadrant_plot_cog_vs_phy.svg", plot = p,
       width = 9, height = 6)

cat("Saved: quadrant_plot_cog_vs_phy.svg\n")

print(p)













