### Empty memory and workspace ----------------------
gc()
rm(list = ls())

### Load libraries ----------------------
library(dplyr)  
library(pROC)
library(ecospat)
library(ggplot2)
library(caret)

### Load datasets ----------------------
data_uniform <- read.csv("croc_datav1.1.csv")
data_bias5k  <- read.csv("croc_data_biased5k.csv")
data_bias10k <- read.csv("croc_data_biased10k.csv")
datasets <- list(Uniform = data_uniform, Bias_5k = data_bias5k, Bias_10k = data_bias10k)

### Partial ROC function ----------------------
partial_roc_manual <- function(pred, obs, proportion = 0.95, iterations = 1000) {
  auc_full <- pROC::auc(obs, pred)
  auc_partial <- pROC::auc(
    obs, pred,
    partial.auc = c(1, proportion),
    partial.auc.focus = "specificity",
    partial.auc.correct = TRUE
  )
  auc_ratio <- as.numeric(auc_partial / auc_full)
  return(auc_ratio)
}

### Helper function: Continuous Boyce Index
robust_boyce <- function(fit, obs, res = 200, window.w = "default", plot = FALSE) {
  fit <- as.numeric(fit)
  obs <- as.numeric(obs)
  fit <- fit[!is.na(fit)]
  obs <- obs[!is.na(obs)]
  if (length(unique(obs)) < 2) return(NA)
  if (window.w == "default") window.w <- (max(fit) - min(fit)) / 20
  mini <- min(fit)
  maxi <- max(fit)
  vec.mov <- seq(from = mini, to = maxi - window.w, by = (maxi - mini - window.w) / res)
  vec.mov[res + 1] <- vec.mov[res + 1] + 1
  interval <- cbind(vec.mov, vec.mov + window.w)
  f <- apply(interval, 1, function(interval) {
    fit.bin <- fit
    obs.bin <- obs
    fit.bin[fit >= interval[1] & fit <= interval[2]] <- "i"; fit.bin[fit.bin != "i"] <- 0
    obs.bin[obs >= interval[1] & obs <= interval[2]] <- "i"; obs.bin[obs.bin != "i"] <- 0
    pi <- length(which(obs.bin == "i")) / length(obs)
    ei <- length(which(fit.bin == "i")) / length(fit.bin)
    fi <- ifelse(ei == 0, NA, pi / ei)
    return(fi)
  })
  f <- f[!is.na(f) & f != "NaN" & f > 0]
  if (length(f) < 3) return(NA)
  HS <- apply(interval, 1, mean)[1:length(f)]
  b <- cor(f, HS[1:length(f)], method = "spearman")
  return(round(b, 3))
}


### Main evaluation loop ----------------------
results <- list()

for (bg in names(datasets)) {
  df <- datasets[[bg]] %>% filter(complete.cases(.))
  df$presence <- as.numeric(df$presence)
  ms <- glm(
    formula = presence ~ FlowAcc_FM_8000 + LULC_Base +
      TWI_SD_500 + FlowAcc_Base + DistWater_Base +
      DistSettlement_FM_4000 + LULC_CWED_1000 +
      DistWater_SD_8000 + bio3_Base + SolarRad_SD_2000 +
      DistRoad_Base + Wetland_Base + Wetland_GYRATE_AM_1000,
    family = "binomial",
    data = df
  )
  preds <- predict(ms, type = "response")
  auc_val <- pROC::auc(df$presence, preds)
  partial_roc_val <- partial_roc_manual(preds, df$presence, proportion = 0.95)
  omission_rate <- mean(preds[df$presence == 1] < 0.5)
  pred_binary <- ifelse(preds > 0.5, 1, 0)
  conf <- tryCatch({
    caret::confusionMatrix(
      factor(pred_binary, levels = c(0, 1)),
      factor(df$presence, levels = c(0, 1)),
      positive = "1"
    )
  }, error = function(e) NULL)
  tss_val <- if (!is.null(conf)) {
    as.numeric(conf$byClass["Sensitivity"] + conf$byClass["Specificity"] - 1)
  } else {
    NA
  }
  boyce_val <- tryCatch({
    robust_boyce(fit = preds, obs = preds[df$presence == 1])
  }, error = function(e) NA)
  metrics <- data.frame(
    AIC = AIC(ms),
    AUC = as.numeric(auc_val),
    PartialROC = partial_roc_val,
    OmissionRate = omission_rate,
    TSS = as.numeric(tss_val),
    Boyce = as.numeric(boyce_val),
    DevianceExplained = 1 - (ms$deviance / ms$null.deviance),
    Background = bg
  )
  results[[bg]] <- metrics
  cat("Model run complete for:", bg, "\n")
}


### Combine results ----------------------
results_df <- do.call(rbind, results)
print(results_df)


### Normalization and composite scoring ----------------------
results_df <- results_df %>%
  mutate(
    AUC_scaled = (AUC - min(AUC)) / (max(AUC) - min(AUC)),
    Boyce_scaled = (Boyce - min(Boyce)) / (max(Boyce) - min(Boyce)),
    TSS_scaled = (TSS - min(TSS)) / (max(TSS) - min(TSS)),
    PartialROC_scaled = (PartialROC - min(PartialROC)) / (max(PartialROC) - min(PartialROC)),
    AIC_scaled = 1 - (AIC - min(AIC)) / (max(AIC) - min(AIC)),
    predictive_score = (AUC_scaled + Boyce_scaled + TSS_scaled + PartialROC_scaled + AIC_scaled) / 5,
    realism_score = ifelse(grepl("Uniform", Background), 0, 1)
  )


### Sensitivity analysis ----------------------
realism_weights <- seq(0, 0.5, by = 0.05)
sensitivity_df <- expand.grid(Background = results_df$Background, realism_weight = realism_weights) %>%
  left_join(results_df, by = "Background") %>%
  mutate(Composite = predictive_score * (1 - realism_weight) + realism_score * realism_weight) %>%
  group_by(realism_weight) %>%
  mutate(Composite_norm = (Composite - min(Composite)) / (max(Composite) - min(Composite))) %>%
  ungroup()


### Plot sensitivity curves ----------------------
sensitity_uniform_5k_10k <- ggplot(sensitivity_df, aes(x = realism_weight, y = Composite_norm, color = Background)) +
  geom_line(size = 0.8) +
  geom_point(size = 2) +
  labs(
    title = "Sensitivity of Model Selection to Ecological Realism Weight",
    subtitle = "Trade-off between predictive accuracy and ecological realism",
    x = "Ecological realism weight (0 = purely predictive, 1 = purely realistic)",
    y = "Normalized composite score",
    color = "Background type"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.3),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 12, color = "gray25"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    legend.text = element_text(size = 11)
  )
ggsave("sensitity_uniform_5k_10k.jpeg", sensitity_uniform_5k_10k, width = 8, height = 5, dpi = 300)