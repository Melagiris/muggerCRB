## Load data for GLM spatial CV
glm_data <- read.csv('/Users/argos/Desktop/Croc Modeling Cauvery River Basin/R/Revision3analysis/BiasWeighted_BGpoints_5k_env_latestwithspThin.CSV')

## required libraries
library(blockCV) 
library(pROC) 
library(caret)
library(dplyr) 
library(purrr)

## Convert to sf object
glm_sf <- st_as_sf(glm_data, coords = c("Longitude", "Latitude"), crs = 4326)

## Create spatial blocks (~20 km range)
set.seed(42)
sb <- spatialBlock(
  speciesData = glm_sf,
  species = "presence",
  theRange = 20000,
  k = 5,
  selection = "random",
  iteration = 100,
  showBlocks = TRUE,
  progress = TRUE
)

folds <- sb$foldID

## Partial ROC function
partial_roc_manual <- function(pred, obs, proportion = 0.95, reps = 1000) {
  roc_obj <- pROC::roc(obs, pred)
  auc_full <- pROC::auc(roc_obj)
  tpr_vals <- roc_obj$sensitivities
  fpr_vals <- 1 - roc_obj$specificities
  auc_partial <- pROC::auc(fpr_vals[fpr_vals <= proportion],
                           tpr_vals[fpr_vals <= proportion])
  return(as.numeric(auc_partial / auc_full))
}

## Spatial CV results
spatial_results <- list()

for (i in 1:max(folds, na.rm = TRUE)) {
  train <- glm_sf[folds != i, ]
  test  <- glm_sf[folds == i, ]
  
  if (length(unique(test$presence)) < 2) {
    cat("Fold", i, "skipped (single class)\n")
    next
  }
  
  model <- glm(
    presence ~ FlowAcc_FM_8000 + LULC_Base +
      TWI_SD_500 + FlowAcc_Base + DistWater_Base +
      DistSettlement_FM_4000 + LULC_CWED_1000 +
      DistWater_SD_8000 + bio3_Base + SolarRad_SD_2000 +
      DistRoad_Base + Wetland_Base + Wetland_GYRATE_AM_1000,
    family = "binomial",
    data = train
  )
  
  preds <- predict(model, newdata = test, type = "response")
  
  auc_val <- tryCatch(as.numeric(pROC::auc(test$presence, preds)), error = function(e) NA)
  partial_roc_val <- tryCatch(partial_roc_manual(preds, test$presence), error = function(e) NA)
  omission_rate <- mean(preds[test$presence == 1] < 0.5, na.rm = TRUE)
  
  if (length(unique(test$presence)) == 2) {
    preds_class <- ifelse(preds >= 0.5, 1, 0)
    tab <- table(factor(preds_class, levels = c(0,1)),
                 factor(test$presence, levels = c(0,1)))
    TP <- tab[2,2]; TN <- tab[1,1]; FP <- tab[2,1]; FN <- tab[1,2]
    sens <- TP / (TP + FN)
    spec <- TN / (TN + FP)
    tss_val <- sens + spec - 1
  } else {
    tss_val <- NA
  }
  
  dev_exp <- 1 - (model$deviance / model$null.deviance)
  
  spatial_results[[i]] <- data.frame(
    Fold = i,
    AUC = auc_val,
    PartialROC = partial_roc_val,
    OmissionRate = omission_rate,
    TSS = as.numeric(tss_val),
    DevianceExplained = dev_exp
  )
  
  cat("Spatial fold", i, "complete\n")
}

spatial_cv <- bind_rows(spatial_results)

spatial_summary <- spatial_cv %>%
  summarise(across(AUC:DevianceExplained, mean, na.rm = TRUE))

spatial_summary

## Naive CV (random)
set.seed(123)
nfolds <- 5
folds_random <- sample(rep(1:nfolds, length.out = nrow(glm_sf)))

naive_results <- list()

for (i in 1:nfolds) {
  train <- glm_sf[folds_random != i, ]
  test  <- glm_sf[folds_random == i, ]
  
  model <- glm(
    presence ~ FlowAcc_FM_8000 + LULC_Base + DistWater_Base +
      TWI_SD_500 + TWI_FM_8000 + Wetland_Base + Wetland_GYRATE_AM_4000 +
      Wetland_GYRATE_AM_2000 + bio3_Base + LULC_CWED_500 + HydroConditionedDEM_Base +
      DistRoad_Base + LULC_CWED_2000 + TWI_SD_8000 + bio7_FM_8000 + bio4_FM_2000,
    family = "binomial",
    data = train
  )
  
  preds <- predict(model, newdata = test, type = "response")
  
  auc_val <- tryCatch(as.numeric(pROC::auc(test$presence, preds)), error = function(e) NA)
  partial_roc_val <- tryCatch(partial_roc_manual(preds, test$presence), error = function(e) NA)
  omission_rate <- mean(preds[test$presence == 1] < 0.5, na.rm = TRUE)
  
  if (length(unique(test$presence)) == 2) {
    preds_class <- ifelse(preds >= 0.5, 1, 0)
    tab <- table(factor(preds_class, levels = c(0,1)),
                 factor(test$presence, levels = c(0,1)))
    TP <- tab[2,2]; TN <- tab[1,1]; FP <- tab[2,1]; FN <- tab[1,2]
    sens <- TP / (TP + FN)
    spec <- TN / (TN + FP)
    tss_val <- sens + spec - 1
  } else {
    tss_val <- NA
  }
  
  dev_exp <- 1 - (model$deviance / model$null.deviance)
  
  naive_results[[i]] <- data.frame(
    Fold = i,
    AUC = auc_val,
    PartialROC = partial_roc_val,
    OmissionRate = omission_rate,
    TSS = tss_val,
    DevianceExplained = dev_exp
  )
  
  cat("Naive fold", i, "complete\n")
}

## Combine results
spatial_df <- do.call(rbind, spatial_results) %>% mutate(Type = "Spatial")
naive_df   <- do.call(rbind, naive_results) %>% mutate(Type = "Naive")
comparison_df <- rbind(spatial_df, naive_df)

## Summary
comparison_summary <- comparison_df %>%
  group_by(Type) %>%
  summarise(across(AUC:DevianceExplained,
                   list(mean = mean, sd = sd), na.rm = TRUE),
            .groups = "drop")

comparison_summary

comparison_summary_fmt <- comparison_summary %>%
  mutate(across(ends_with("_mean"), round, 3),
         across(ends_with("_sd"), round, 3)) %>%
  mutate(across(AUC_mean:DevianceExplained_mean,
                ~ paste0(.x, " ± ", get(sub("_mean$", "_sd", cur_column())))))

comparison_summary_fmt %>%
  dplyr::select(Type, dplyr::ends_with("_mean")) %>%
  dplyr::rename_with(~ gsub("_mean", "", .x))

## Prepare data for violin plots
library(tidyr)
library(ggplot2)

comparison_long <- comparison_df %>%
  dplyr::select(Type, PartialROC, TSS, OmissionRate, DevianceExplained) %>%
  tidyr::pivot_longer(cols = -Type, names_to = "Metric", values_to = "Value")

summary_labels <- comparison_long %>%
  group_by(Type, Metric) %>%
  summarise(
    mean_val = mean(Value, na.rm = TRUE),
    sd_val = sd(Value, na.rm = TRUE),
    ymax = max(Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("\n%.2f ± %.2f", mean_val, sd_val),
         y = ymax * 1.006)

vPlot <- ggplot(comparison_long, aes(x = Type, y = Value)) +
  geom_violin(fill = NA, color = "black", linewidth = 0.3, width = 0.9, trim = FALSE) +
  geom_boxplot(width = 0.12, color = "black", fill = NA,
               outlier.shape = 21, outlier.size = 1.3, outlier.stroke = 0.3, linewidth = 0.3) +
  geom_errorbar(
    data = summary_labels,
    aes(x = Type, ymin = mean_val - sd_val, ymax = mean_val + sd_val),
    width = 0.1, linewidth = 0.3, color = "black", inherit.aes = FALSE
  ) +
  geom_point(
    data = summary_labels,
    aes(x = Type, y = mean_val),
    shape = 21, fill = "white", color = "black", size = 2, stroke = 0.3,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = summary_labels,
    aes(x = Type, y = y, label = label),
    size = 3.2, fontface = "italic", color = "black",
    inherit.aes = FALSE
  ) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 2) +
  labs(x = "Cross-validation type", y = "Metric value") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 11),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    panel.spacing = unit(1, "lines"),
    plot.margin = margin(15, 15, 10, 10)
  )

vPlot

ggsave("PerformanceTestforspatialCVvsNaiveCVPlain.jpeg",
       plot = vPlot, width = 12, height = 7, dpi = 300)