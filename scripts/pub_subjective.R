library(readxl)
library(tidyverse)
library(RColorBrewer)
library(see)
library(mdthemes)
library(ggtext)
library(ggrepel)
library(ggplot2)
library(gridExtra)
library(grid)
library(ggpubr)
library(ggsci)
library(cowplot)
library(ggpp)
library(patchwork)
library(Hmisc)
library(writexl)
library(ggstatsplot)
library(emmeans)
library(nlme)
library(lmerTest)
library(janitor)
library(dplyr)
library(plyr, include.only = "ddply")

# ---- Load and prepare data ----
#master_df <- read_excel("/Users/administrator/Documents/MATLAB/project_cb_tom/data/behaviour/peq_data_linear_model.xlsx")
#master_df <- read_excel("/Users/administrator/Documents/MATLAB/project_cb_tom/data/behaviour/significant_behaviour.xlsx")
master_df <- read_excel("/Users/administrator/Documents/MATLAB/project_cb_tom/data/behaviour/paper_outcomes.xlsx")

# Exclude participants 150 and 170
master_df <- master_df[!(master_df$participant %in% c(15, 17)), ]

# Add grouped drug variable
master_df$drug_grouped <- ifelse(master_df$drug == "cPlacebo", "Placebo", "Psychedelic")
master_df$drug_grouped <- factor(master_df$drug_grouped, levels = c("Placebo", "Psychedelic"))

#outcome_measures <- names(master_df)[4:(length(master_df) - 11)] #for peq
outcome_measures <- names(master_df)[4:(length(master_df)-1)] 

# Initialize outputs
outputs_model <- data.frame(matrix(nrow = 1, ncol = 11))
outputs_emmeans <- data.frame(matrix(nrow = 1, ncol = 7))
outputs_pairwise <- data.frame(matrix(nrow = 1, ncol = 12))

p_threshold <- 0.05

for (i in outcome_measures) {
  print(paste0("Running grouped model for: ", i))
  tmp_name <- t(as.data.frame(i, 'variable'))
  
  # --- Grouped model (Placebo vs Psychedelic) ---
  m1_grouped <- lme(
    as.formula(paste(i, "~ drug_grouped")),
    random = ~1 | participant,
    data = master_df,
    method = "REML",
    na.action = na.omit
  )
  
  tmp_bic <- as.data.frame(summary(m1_grouped)$BIC)
  tmp_resid <- t(as.data.frame(summary(m1_grouped)$residuals))
  tmp_anova <- anova(m1_grouped)[2, ]
  
  tmp_model_results <- cbind(tmp_name, tmp_anova, tmp_bic, tmp_resid)
  rownames(tmp_model_results) <- NULL
  
  # Save model outputs
  names(outputs_model) <- names(tmp_model_results)
  outputs_model <- rbind(outputs_model, tmp_model_results)
  
  # If significant main effect, run detailed EMMeans using full drug factor
  if (tmp_model_results$`p-value` < p_threshold) {
    print(paste0("Main effect p = ", round(tmp_model_results$`p-value`, 4), " found! (Grouped model)"))
    
    # --- Full model with 3-level drug factor ---
    m1_full <- lme(
      as.formula(paste(i, "~ drug")),
      random = ~1 | participant,
      data = master_df,
      method = "REML",
      na.action = na.omit
    )
    
    tmp_posthoc <- emmeans(m1_full, list(pairwise ~ drug), adjust = "Tukey")
    
    # Estimated marginal means
    tmp_emeans <- summary(tmp_posthoc)$emmeans
    tmp_emeans_results <- cbind(as.data.frame(rep(i, 3)), tmp_emeans)
    names(tmp_emeans_results)[1] <- "variable"
    
    # Pairwise comparisons
    tmp_pairwise <- summary(tmp_posthoc)$pairwise
    
    # Effect sizes (Cohen's d)
    tmp_cohensd <- eff_size(tmp_posthoc, sigma = sigma(m1_full), edf = anova(m1_full)$denDF[2])
    tmp_cohensd_df <- as.data.frame(tmp_cohensd)
    
    # Uncorrected p-values
    tmp_uncorrected_p <- summary(
      emmeans(m1_full, list(pairwise ~ drug), adjust = "none")
    )$pairwise[6]
    
    # Combine outputs
    tmp_pairwise_final <- cbind(
      as.data.frame(rep(i, 3)),
      tmp_pairwise,
      tmp_cohensd_df[c("effect.size", "SE", "lower.CL", "upper.CL")],
      tmp_uncorrected_p
    )
    names(tmp_pairwise_final)[1] <- "variable"
    names(tmp_pairwise_final)[12] <- "p.value_uncorrected"
    tmp_pairwise_final <- tmp_pairwise_final %>% relocate(p.value_uncorrected, .after = p.value)
    
    # Format p-values (no scientific notation)
    tmp_pairwise_final$p.value <- format(as.numeric(gsub("<", "", tmp_pairwise_final$p.value)), scientific = FALSE, digits = 4)
    tmp_pairwise_final$p.value_uncorrected <- format(as.numeric(gsub("<", "", tmp_pairwise_final$p.value_uncorrected)), scientific = FALSE, digits = 4)
    
    # Append results
    names(outputs_emmeans) <- names(tmp_emeans_results)
    outputs_emmeans <- rbind(outputs_emmeans, tmp_emeans_results)
    
    names(outputs_pairwise) <- names(tmp_pairwise_final)
    outputs_pairwise <- rbind(outputs_pairwise, tmp_pairwise_final)
  }
}

# --- CLEAN OUTPUTS ---
outputs_model <- outputs_model[-1, ]
outputs_emmeans <- outputs_emmeans[-1, ]
outputs_pairwise <- outputs_pairwise[-1, ]

# === Format Means with Asterisks and Hashtags ===
all_means <- master_df %>%
  select(drug, outcome_measures) %>%
  pivot_longer(cols = 2:length(.), names_to = "variable", values_to = "score")

all_means_data <- ddply(all_means, c("drug", "variable"), summarise,
                        N = sum(!is.na(score)),
                        mean = mean(score, na.rm = TRUE),
                        se = sd(score, na.rm = TRUE) / sqrt(N))
all_means_data$mean_se <- paste0(round(all_means_data$mean, 2), "±", round(all_means_data$se, 2))
all_means_data <- subset(all_means_data, select = -c(mean, se, N))
all_means_data <- spread(all_means_data, key = drug, value = mean_se)
all_means_data <- merge(all_means_data, outputs_model, by = 'variable')

# Add significance flags (* for active vs placebo; # for Psilo vs 2C-B)
p_index <- outputs_pairwise[, c("variable", "1", "p.value")]
p_index$flag <- ""

for (i in 1:nrow(p_index)) {
  # Clean p-values
  p_val_clean <- suppressWarnings(as.numeric(gsub("<", "", p_index[i, "p.value"])))
  if (is.na(p_val_clean)) next
  
  # Flag type
  flag <- ifelse(p_index[i, "1"] != "aPsilocybin - (b2C-B)", "*", "#")
  
  # Threshold-based flags
  if (p_val_clean < 0.05 & p_val_clean > 0.01) {
    p_index[i, "flag"] <- flag
  } else if (p_val_clean < 0.01 & p_val_clean > 0.001) {
    p_index[i, "flag"] <- paste0(flag, flag)
  } else if (p_val_clean < 0.001) {
    p_index[i, "flag"] <- paste0(flag, flag, flag)
  } else {
    p_index[i, "flag"] <- NA
  }
}

# Apply flags to means
p_index_input <- p_index[p_index$`1` != "aPsilocybin - (b2C-B)", ]
for (i in 1:nrow(p_index_input)) {
  if (grepl("Psilocybin", p_index_input[i, "1"])) {
    all_means_data[all_means_data$variable == p_index_input[i, "variable"], "aPsilocybin"] <-
      paste0(all_means_data[all_means_data$variable == p_index_input[i, "variable"], "aPsilocybin"], p_index_input[i, "flag"])
  } else {
    all_means_data[all_means_data$variable == p_index_input[i, "variable"], "b2C-B"] <-
      paste0(all_means_data[all_means_data$variable == p_index_input[i, "variable"], "b2C-B"], p_index_input[i, "flag"])
  }
}

# Add # flags (Psilo vs 2C-B)
p_index_input <- p_index[p_index$`1` == "aPsilocybin - (b2C-B)", ]
for (i in 1:nrow(p_index_input)) {
  all_means_data[all_means_data$variable == p_index_input[i, "variable"], "aPsilocybin"] <-
    paste0(all_means_data[all_means_data$variable == p_index_input[i, "variable"], "aPsilocybin"], p_index_input[i, "flag"])
}

# Final cleaning
all_means_data[] <- lapply(all_means_data, gsub, pattern = 'NA', replacement = '')

# Export (if desired)
out_names = list('model_tidy' = all_means_data,'model_pairwise' = outputs_pairwise,'model_emmeans' = outputs_emmeans,'input_data' = master_df)
#write_xlsx(out_names, "/Users/administrator/Documents/MATLAB/project_cb_tom/data/behaviour/single_stats.xlsx")


