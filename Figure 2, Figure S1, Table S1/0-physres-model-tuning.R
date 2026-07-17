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
library(reshape2)
library(ppcor)
library(lme4)
library(lmerTest)       
library(broom.mixed)   
library(future.apply) 
library(patchwork)
library(broom)
library(forcats)
library(flextable)
library(gtsummary)
library(gt)      
library(zCompositions) # for zero replacement
library(compositions) # for clr transform
library(pls) # for PLS regression (plsr)
library(pheatmap)
library(MASS)


df <- readRDS("clean_data_wellcome.rds")

df$host_age<-as.numeric(df$host_age)
df$host_body_mass_index<-as.numeric(df$host_body_mass_index)
df$physresilience_pctl<-as.numeric(df$physresilience_pctl)
df$physresilience<-as.numeric(df$physresilience)


df <- df %>%
  group_by(record_ID) %>%
  slice_min(order_by = host_age, with_ties = FALSE) %>%
  ungroup()


#Selction of metabolites only
df_metab <- df %>%
  column_to_rownames(var = "filename") %>%  
  select(3:22663) %>%                
  filter(rowSums(.) > 0)           

#RCLR conversion
df_rclr <- decostand(df_metab, method = "rclr")
df_rclr_matrix <- as.matrix(df_rclr)

#Removing zero variance metabolites
var_resid <- apply(df_rclr_matrix, 2, var)
metab_resid_filtered <- df_rclr_matrix[, var_resid > 1e-6]


X <- metab_resid_filtered
Y <- df$physresilience

var_resid <- apply(metab_resid_filtered, 2, var)
top_idx <- order(var_resid, decreasing = TRUE)[1:5000]
X_sub <- metab_resid_filtered[, top_idx]

set.seed(42)

tune <- tune.spls(
  X = X_sub,
  Y = Y,
  ncomp = 5,
  test.keepX = c(50, 100, 200, 500),
  measure = "MSE",        # <--- specify explicitly
  validation = "Mfold",
  folds = 5,
  nrepeat = 5,
  progressBar = TRUE
)

tune$choice.ncomp


#Running model for later inflection
spls_final <- spls(X_sub, Y, ncomp = 1, keepX = 5000) 

# Extract absolute weights from sPLS model
loadings <- abs(spls_final$loadings$X[,1])


# Rank them
df_ranked <- data.frame(
  Metabolite = names(loadings),
  Weight = loadings
) %>%
  arrange(desc(Weight)) %>%
  mutate(Rank = row_number())


x <- 1:length(df_ranked$Weight) 
y <- df_ranked$Weight

elbows <- inflection::findiplist(x, y, index = 1)
elbow <- elbows[length(elbows)]

cutoff <- y[elbow]
selected <- sum(y > cutoff)

cat("Elbow at metabolite rank:", elbow, "\n")
cat("Cutoff weight:", cutoff, "\n")
cat("Number of metabolites above cutoff:", selected, "\n")


summary(y)
any(is.na(y))
length(unique(y))


