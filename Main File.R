### Empty memory and workspace ----------------------
gc()
rm(list=ls())

### Load relevant libraries ----------------------
library(raster)
library(terra) 
library(spThin)
library(corrplot)
library(mecofun)
library(caret)
library(cvAUC)
library(pROC)
library(sp)
library(sf)
library(dismo)
library(lattice)
library(ggplot2)
library(ecospat)
library(MASS)
library(dplyr)
library(tidyr)
library(glmnet)
library(effects)
library(gurobi)
library(prioritizr)
library(blockCV)
library(pROC)

# loading helper function: Continous Boyce Index
robust_boyce <- function(fit, obs, res = 200, window.w = "default", plot = FALSE) {
  fit <- as.numeric(fit)
  obs <- as.numeric(obs)
  fit <- fit[!is.na(fit)]
  obs <- obs[!is.na(obs)]
  
  if (length(unique(obs)) < 2) return(NA)
  
  if (window.w == "default") {
    window.w <- (max(fit) - min(fit)) / 20
  }
  
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
  
  if (plot) {
    plot(HS, f, xlab = "Habitat suitability", ylab = "Predicted/Expected ratio")
    abline(h = 1, col = "red", lty = 2)
  }
  
  return(round(b, 3))
}


### load data
croc_data <- read.csv('/Users/argos/Desktop/Croc Modeling Cauvery River Basin/R/Revision3analysis/BiasWeighted_BGpoints_5k_env_latestwithspThin.CSV')
coords <- croc_data[, c("Longitude", "Latitude")]
spdf <- SpatialPointsDataFrame(coords, data = croc_data,
                               proj4string = CRS("+proj=longlat +datum=WGS84"))

### create spatial folds
set.seed(123)
jpeg("SpatialCV20kmblocks5folds.jpeg", width = 10, height = 7, units = "in", res = 300)
cv_blocks <- blockCV::cv_spatial(
  x = spdf,              # spatial object (sf, SpatVector, or Spatial*)
  column = "presence",   # column containing 0/1 values
  k = 5,                 # number of folds
  size = 20000,          # block size in meters (~20 km)
  selection = "systematic",  # or "random"
  seed = 123,
  biomod2 = FALSE,
  plot = TRUE
)
dev.off()

table(cv_blocks$folds_ids)

# Extract fold assignments
fold_id_full <- cv_blocks$folds_ids


### define variable groups
scales <- c("Base", "1000", "2000", "4000", "8000", "Multiscale")

# multiscale variables (extracted from variable selection analysis, refer Correlation_And_Variable_Selection.R)
selected_vars_final <- c("FlowAcc_SD_4000", "FlowAcc_FM_8000", "FlowAcc_SD_2000",
                         "LULC_Base", "FlowAcc_SD_1000", "TWI_SD_500", "SolarRad_SD_1000",
                         "FlowAcc_Base", "DistWater_Base", "Wetland_GYRATE_AM_2000",
                         "Wetland_GYRATE_AM_4000", "CC_Edge_8000", "AI_SD_1000",            
                         "TWI_FM_8000", "Wetland_GYRATE_AM_1000", "TWI_SD_2000", "LULC_CWED_1000",       
                         "DistRoad_FM_8000", "DistWater_SD_8000", "bio9_FM_2000", "TWI_Base",            
                         "bio3_Base", "CC_Edge_2000", "DistRoad_Base", "AI_Base",            
                         "SolarRad_SD_2000", "DistWater_SD_1000", "bio4_FM_2000", "bio18_FM_8000",      
                         "TWI_SD_8000", "bio4_FM_4000", "bio7_FM_8000", "TWI_SD_4000",        
                         "bio15_FM_8000", "CC_Base", 'HydroConditionedDEM_Base',
                         # add the previously dropped ones
                         "Wetland_Base", "Wetland_GYRATE_AM_500"
)


# loop through each scale
results <- list()
models_per_scale <- list()

for (scale in scales) {
  message("\n==============================")
  message("Processing scale: ", scale)
  
  # Subset predictors for this scale
  predictor_vars <- if (scale == "Multiscale") {
    croc_data[, selected_vars_final, drop = FALSE]
  } else {
    croc_data[, grep(scale, colnames(croc_data)), drop = FALSE]
  }
  
  glm_data <- data.frame(
    presence = croc_data$presence,
    predictor_vars,
    Longitude = croc_data$Longitude,
    Latitude  = croc_data$Latitude
  )
  glm_data <- na.omit(glm_data)
  
  # Assign folds
  glm_data$fold <- cv_blocks$folds_ids[1:nrow(glm_data)]
  
  # Metrics
  auc_vals <- c()
  boyce_vals <- c()
  dev_exp <- c()
  aic_vals <- c()
  models <- list() 
  
  for (i in 1:5) {
    train_data <- glm_data[glm_data$fold != i, ]
    test_data  <- glm_data[glm_data$fold == i, ]
    
    # ---- Fit GLM full model ----
    glm_full <- glm(
      presence ~ .,
      data = train_data[, !names(train_data) %in% c("Longitude", "Latitude", "fold")],
      family = "binomial"
    )
    
    # ---- Stepwise selection ----
    best_model <- stepAIC(glm_full, direction = "both", trace = FALSE)
    models[[i]] <- best_model
    
    # ---- Predictions ----
    preds <- predict(best_model, newdata = test_data, type = "response")
    
    # ---- Metrics ----
    roc_curve <- pROC::roc(test_data$presence, preds)
    auc_vals[i] <- pROC::auc(roc_curve)
    
    boyce_vals[i] <- tryCatch({
      robust_boyce(fit = preds, obs = preds[test_data$presence == 1])
    }, error = function(e) NA)
    
    dev_exp[i] <- 1 - (best_model$deviance / best_model$null.deviance)
    aic_vals[i] <- AIC(best_model)
  }
  
  # ---- Summarize per scale ----
  results[[scale]] <- data.frame(
    Scale = scale,
    AUC = mean(auc_vals, na.rm = TRUE),
    Boyce = mean(boyce_vals, na.rm = TRUE),
    DevianceExplained = mean(dev_exp, na.rm = TRUE),
    AIC = mean(aic_vals, na.rm = TRUE)
  )
  
  models_per_scale[[scale]] <- models
}


# combine results & select best scale
results_df <- do.call(rbind, results)
print(results_df)

best_scale <- results_df$Scale[which.min(results_df$AIC)]
cat("\n Best scale by AIC:", best_scale, "\n")


# refit
if (best_scale == "Multiscale") {
  final_vars <- selected_vars_final
} else {
  final_vars <- grep(best_scale, colnames(croc_data), value = TRUE)
}

final_data <- na.omit(data.frame(
  presence = croc_data$presence,
  croc_data[, final_vars, drop = FALSE]
))

final_full <- glm(presence ~ ., data = final_data, family = "binomial")
final_best_model <- stepAIC(final_full, direction = "both", trace = FALSE)

summary(final_best_model)


## ## ## ## ## ## Prediction ## ## ## ## ## ## 
# extracting variables from the final model
final_best_vars <- all.vars(formula(final_best_model))[-1]
cat("Variables:", paste(final_best_vars, collapse = ", "), "\n")

# raster stack for prediction map
env_stack <- stack(list.files(path = "/Users/argos/Desktop/Croc Modeling Cauvery River Basin/R/Rasters_Final/finalContinousAndCategorical", 
                              pattern = "\\.tif$", full.names = TRUE))
# subset raster stack to include only variables in the final model
final_raster_stack <- subset(env_stack, final_best_vars)

# predicting probabilities for each cell in the raster stack
predicted_distribution <- predict(final_raster_stack, final_best_model, type = "response")
plot(predicted_distribution)
writeRaster(predicted_distribution, "predicted_distribution_06Nov2025.tif", overwrite = TRUE)


## find total predicted area--start
## to UTM Zone 43N
r_utm <- project(predicted_distribution, "EPSG:32643")
r <- r_utm
# Pixel area (km²)
pixel_area_km2 <- (res(r)[1]^2) / 1e6
# Extract values
vals <- values(r)
# Exclude tiny predictions (< 0.05)
threshold <- 0.05
n_above_thresh <- sum(vals > threshold, na.rm = TRUE)
# Calculate total area above threshold
area_above_thresh_km2 <- n_above_thresh * pixel_area_km2
## find total predicted area--end


#### predicted area inside protected area ---start
protected_area <- st_read('/Users/argos/Desktop/Croc Modeling Cauvery River Basin/CRB BB Forest Beat Boundaries Forest Survey Of India/CRB FSI beat boundary/CRB_BeatBoundary_FSI.shp')

# reproject to UTM
protected_area_utm <- st_transform(protected_area, "EPSG:32643")
# area (in m²), then convert to km²
protected_area_utm$area_km2 <- as.numeric(st_area(protected_area_utm)) / 1e6
# total protected area
total_protected_area <- sum(protected_area_utm$area_km2, na.rm = TRUE)

# reproject the protected area to match the raster CRS (UTM)
predicted_area <- rast("/Users/argos/Desktop/Croc Modeling Cauvery River Basin/R/Revision3analysis/predicted_distribution_06Nov2025.tif")
protected_area_utm <- st_transform(protected_area, crs(r))
suitable_mask <- r > 0.05
# convert sf to SpatVector for terra operations
protected_vect <- vect(protected_area_utm)
# Mask suitable raster by protected areas
suitable_in_pa <- mask(suitable_mask, protected_vect)
# calculate pixel area in km²
pixel_area_km2 <- (res(r)[1]^2) / 1e6
# count suitable pixels total and inside PA
total_suitable_pixels <- sum(values(suitable_mask), na.rm = TRUE)
suitable_pixels_in_pa <- sum(values(suitable_in_pa), na.rm = TRUE)
# compute areas
total_suitable_area <- total_suitable_pixels * pixel_area_km2
suitable_area_in_pa <- suitable_pixels_in_pa * pixel_area_km2
# calculate percentage
percent_in_pa <- (suitable_area_in_pa / total_suitable_area) * 100
#### predicted area inside protected area ---end


### calculate the proportion of wetlands protected and the proportion of suitable mugger habitat overlapping wetlands and PAs
# everything in same CRS
cauvery_wetland <- project(cauvery_wetland, "EPSG:32643") #'/Users/argos/Desktop/Croc Modeling Cauvery River Basin/R/InputOutputFile24Aug/prioritizr/wetlandCompleteUTM.tif'
protected_area_utm <- st_transform(protected_area, "EPSG:32643")
protected_vect <- vect(protected_area_utm)
cauvery_wetland_res <- terra::resample(cauvery_wetland, r, method = "near")

# thresholds
wetland_threshold <- 0.05
suitability_threshold <- 0.05

# binary masks
wetland_mask <- cauvery_wetland_res > wetland_threshold
suitable_mask <- r > suitability_threshold

# pixel area (km²)
pixel_area_km2 <- (res(r)[1]^2) / 1e6

# total wetland area
wetland_pixels <- sum(values(wetland_mask), na.rm = TRUE)
wetland_area_km2 <- wetland_pixels * pixel_area_km2

# wetlands inside protected areas
wetland_in_pa <- mask(wetland_mask, protected_vect)
wetland_in_pa_pixels <- sum(values(wetland_in_pa), na.rm = TRUE)
wetland_in_pa_km2 <- wetland_in_pa_pixels * pixel_area_km2
wetland_pa_percent <- (wetland_in_pa_km2 / wetland_area_km2) * 100

# suitable mugger habitat overlapping wetlands
suitable_wetland_overlap <- wetland_mask & suitable_mask
suitable_wetland_pixels <- sum(values(suitable_wetland_overlap), na.rm = TRUE)
suitable_wetland_area_km2 <- suitable_wetland_pixels * pixel_area_km2

# suitable mugger habitat in wetlands inside PAs
suitable_wetland_in_pa <- mask(suitable_wetland_overlap, protected_vect)
suitable_wetland_in_pa_pixels <- sum(values(suitable_wetland_in_pa), na.rm = TRUE)
suitable_wetland_in_pa_km2 <- suitable_wetland_in_pa_pixels * pixel_area_km2
suitable_wetland_pa_percent <- (suitable_wetland_in_pa_km2 / suitable_wetland_area_km2) * 100


## find predicted area category wise
r_utm <- project(predicted_distribution, "EPSG:32643")
r <- r_utm
# catgeorise
high_thresh <- 0.6
moderate_thresh <- 0.3
low_thresh <- 0.1

# suitability masks
high <- r_utm > high_thresh
moderate <- r_utm <= high_thresh & r_utm > moderate_thresh
low <- r_utm <= moderate_thresh & r_utm > low_thresh
unsuitable <- r_utm <= low_thresh  # optional

# to calculate area (km²) 
calc_area_km2 <- function(mask) {
  cell_area <- prod(res(mask)) / 1e6   # convert m² → km²
  area <- global(mask, "sum", na.rm = TRUE)[1, 1] * cell_area
  return(area)
}

# total valid area
total_area_km2 <- global(!is.na(r_utm), "sum", na.rm = TRUE)[1, 1] * prod(res(r_utm)) / 1e6

# talculate class areas
high_area_km2 <- calc_area_km2(high)
moderate_area_km2 <- calc_area_km2(moderate)
low_area_km2 <- calc_area_km2(low)
unsuitable_area_km2 <- calc_area_km2(unsuitable)

# combine results
area_summary <- data.frame(
  Class = c("High", "Moderate", "Low", "Unsuitable"),
  Threshold = c("> 0.6", "0.3 – 0.6", "0.1 – 0.3", "< 0.1"),
  Area_km2 = c(high_area_km2, moderate_area_km2, low_area_km2, unsuitable_area_km2)
)

# % of total
area_summary$Percent <- round((area_summary$Area_km2 / total_area_km2) * 100, 2)
print(area_summary)

# i/p table
area_summary <- data.frame(
  Class = c("High", "Moderate", "Low"),
  Threshold = c("> 0.6", "0.3 – 0.6", "0.1 – 0.3"),
  Area_km2 = c(37.34287, 197.88997, 811.18534)
)

# total area considering only these three
total_suitable_area <- sum(area_summary$Area_km2)

# percentage share for each class
area_summary$Percent <- round((area_summary$Area_km2 / total_suitable_area) * 100, 2)

# total area row (optional)
total_row <- data.frame(
  Class = "Total suitable area",
  Threshold = "",
  Area_km2 = total_suitable_area,
  Percent = 100
)

# combine
area_summary_final <- rbind(area_summary, total_row)
print(area_summary_final)


### extract AIC table
# defining predictor sets
base_vars   <- c("AI_Base")
vars_1000   <- c("AI_FM_1000")
vars_2000   <- c("AI_FM_2000")
vars_4000   <- c("AI_FM_4000")
vars_8000   <- c("AI_FM_8000")

# complete multiscale variable set
multiscale_vars <- c("FlowAcc_SD_4000", "FlowAcc_FM_8000", "FlowAcc_SD_2000",
                     "LULC_Base", "FlowAcc_SD_1000", "TWI_SD_500", "SolarRad_SD_1000",
                     "FlowAcc_Base", "DistWater_Base", "Wetland_GYRATE_AM_2000",
                     "Wetland_GYRATE_AM_4000", "CC_Edge_8000", "AI_SD_1000",            
                     "TWI_FM_8000", "Wetland_GYRATE_AM_1000", "TWI_SD_2000", "LULC_CWED_1000",       
                     "DistRoad_FM_8000", "DistWater_SD_8000", "bio9_FM_2000", "TWI_Base",            
                     "bio3_Base", "CC_Edge_2000", "DistRoad_Base", "AI_Base",            
                     "SolarRad_SD_2000", "DistWater_SD_1000", "bio4_FM_2000", "bio18_FM_8000",      
                     "TWI_SD_8000", "bio4_FM_4000", "bio7_FM_8000", "TWI_SD_4000",        
                     "bio15_FM_8000", "CC_Base", 'HydroConditionedDEM_Base',
                     # add the previously dropped ones
                     "Wetland_Base", "Wetland_GYRATE_AM_500")
# all model sets
scales <- list(
  Base       = base_vars,
  `1000`     = vars_1000,
  `2000`     = vars_2000,
  `4000`     = vars_4000,
  `8000`     = vars_8000,
  Multiscale = multiscale_vars
)

results <- data.frame(
  Scale = character(),
  Fold = integer(),
  AUC = numeric(),
  AIC = numeric(),
  stringsAsFactors = FALSE
)

# loop through scales
for (scale_name in names(scales)) {
  vars <- scales[[scale_name]]
  cat("\n==============================\nProcessing scale:", scale_name, "\n")
  
  for (fold in unique(cv_blocks$folds_ids)) {
    train_idx <- which(cv_blocks$folds_ids != fold)
    test_idx  <- which(cv_blocks$folds_ids == fold)
    
    train_data <- spdf[train_idx, ]
    test_data  <- spdf[test_idx, ]
    
    formula_str <- as.formula(paste("presence ~", paste(vars, collapse = " + ")))
    model <- glm(formula_str, data = train_data, family = binomial)
    
    preds <- predict(model, newdata = test_data, type = "response")
    
    # AUC (quiet)
    roc_obj <- suppressMessages(pROC::roc(test_data$presence, preds))
    auc_val <- as.numeric(pROC::auc(roc_obj))
    
    # AIC
    aic_val <- AIC(model)
    
    results <- rbind(results, data.frame(
      Scale = scale_name,
      Fold = fold,
      AUC = auc_val,
      AIC = aic_val
    ))
  }
}

# mean & SD
summary_stats <- aggregate(cbind(AUC, AIC) ~ Scale, data = results, FUN = mean)
summary_stats <- merge(summary_stats,
                       aggregate(cbind(AUC, AIC) ~ Scale, data = results, FUN = sd),
                       by = "Scale",
                       suffixes = c("_mean", "_sd"))

# ΔAIC and AIC weights
summary_stats$deltaAIC <- summary_stats$AIC_mean - min(summary_stats$AIC_mean)
summary_stats$AICw <- exp(-0.5 * summary_stats$deltaAIC) /
  sum(exp(-0.5 * summary_stats$deltaAIC))

print(summary_stats)


### Response curves ----------------
# extracting variables from the final model
final_best_vars <- all.vars(formula(final_best_model))[-1]
cat("Variables:", paste(final_best_vars, collapse = ", "), "\n")

response_curves <- list()

for (var in final_best_vars) {
  message("Processing variable: ", var)
  
  # Get full value range for this variable
  var_values <- values(env_stack[[var]])
  var_values <- var_values[!is.na(var_values)]
  var_range <- seq(min(var_values), max(var_values), length.out = 100)
  
  # Mean of all other variables
  other_vars <- sapply(setdiff(final_best_vars, var), function(v) {
    mean(values(env_stack[[v]]), na.rm = TRUE)
  })
  
  # Create new data with all other vars fixed at mean
  new_data <- data.frame(matrix(other_vars, nrow = 100, ncol = length(other_vars), byrow = TRUE))
  colnames(new_data) <- setdiff(final_best_vars, var)
  new_data[[var]] <- var_range
  
  # Ensure column order matches model
  new_data <- new_data[, final_best_vars, drop = FALSE]
  
  # Predict probability of presence
  predictions <- predict(final_best_model, newdata = new_data, type = "response")
  
  # Store response curve
  response_curves[[var]] <- data.frame(Value = var_range, Probability = predictions, Variable = var)
}

# Plot Response curves using ggplot
combined_response_curves <- do.call(rbind, response_curves)

rc <- ggplot(combined_response_curves, aes(x = Value, y = Probability)) +
  geom_line(size = 0.7, color = "steelblue") +
  geom_ribbon(aes(ymin = Probability - 0.02, ymax = Probability + 0.02),
              alpha = 0.2, fill = "skyblue") +
  facet_wrap(~Variable, scales = "free", ncol = 3) +
  labs(
    x = "Predictor Value",
    y = "Predicted Habitat Suitability",
    title = "Partial Response Curves (GLM, Best Model)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold", size = 11),
    axis.text = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 11, color = "black"),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  )

## response curves for the top six variable which has shown high variability in the above plots
selected_vars <- c("bio3_Base", "Wetland_GYRATE_AM_1000", "DistWater_SD_8000",
                   "FlowAcc_FM_8000", "TWI_SD_500", "DistRoad_Base")

response_curves <- list()
final_best_vars <- all.vars(formula(final_best_model))[-1]

for (var in selected_vars) {
  message("Processing variable: ", var)
  
  # Extract range of the variable
  var_values <- values(env_stack[[var]])
  var_values <- var_values[!is.na(var_values)]
  var_range <- seq(min(var_values), max(var_values), length.out = 100)
  
  # Mean of all other predictors
  other_vars <- sapply(setdiff(final_best_vars, var), function(v) {
    mean(values(env_stack[[v]]), na.rm = TRUE)
  })
  
  # Create data for predictions
  new_data <- data.frame(matrix(other_vars, nrow = 100, ncol = length(other_vars), byrow = TRUE))
  colnames(new_data) <- setdiff(final_best_vars, var)
  new_data[[var]] <- var_range
  new_data <- new_data[, final_best_vars, drop = FALSE]
  
  # Predict response
  predictions <- predict(final_best_model, newdata = new_data, type = "response")
  
  # Store curve
  response_curves[[var]] <- data.frame(Value = var_range, Probability = predictions, Variable = var)
}

# Combine into one data frame and set factor order
combined_response_curves <- do.call(rbind, response_curves)
combined_response_curves$Variable <- factor(combined_response_curves$Variable, levels = selected_vars)

# Plot
rc <- ggplot(combined_response_curves, aes(x = Value, y = Probability)) +
  geom_line(size = 0.6, color = "steelblue") +
  geom_ribbon(aes(ymin = Probability - 0.02, ymax = Probability + 0.02), alpha = 0.2) +
  labs(x = "", y = "Habitat Suitability") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 8, color = "black"), 
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 12, color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.05, "cm"),
    strip.text = element_text(size = 9, color = "black"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = "black"),
    panel.border = element_rect(fill = NA, color = "black", size = 0.6),
    legend.position = "none"
  ) +
  facet_wrap(~Variable, scales = "free", ncol = 3)

ggsave("ResponseCurves.jpeg", rc, width = 7, height = 6, dpi = 600)


### SCP | prioritizr using gurobi solver
## Calibration (tuning) https://prioritizr.net/articles/calibrating_trade-offs_tutorial.html#data
install.packages(
  file.path(
    Sys.getenv("GUROBI_HOME"),
    "/Library/gurobi1203/macos_universal2/R/gurobi_12.0-3_R_4.5.0.tgz"
  ),
  repos = NULL
)
install.packages("slam", repos = "https://cloud.r-project.org")

# to check the gurobi connection (license should be installed first in the system)
gurobi::gurobi_iis

#### to test
## dummy planning problem to see gurobi solver being utilised by prioritizr
# planning units (with ID)
pu <- data.frame(
  id = 1:5,
  cost = c(1, 2, 3, 4, 5)
)
# features (needs both `id` and `name`)
features <- data.frame(
  id = 1:2,
  name = c("species1", "species2")
)
# rij matrix (links PU ↔ Feature ↔ amount)
rij <- data.frame(
  pu = c(1, 2, 3, 4, 5, 1, 3, 5),
  species = c(1, 1, 1, 1, 1, 2, 2, 2),
  amount = runif(8, 0.1, 0.5)
)
# problem
p <- problem(
  x = pu,
  features = features,
  rij = rij,
  cost_column = "cost"
) %>%
  add_min_set_objective() %>%
  add_relative_targets(0.3) %>%
  add_binary_decisions() %>%
  add_gurobi_solver(gap = 0.01, verbose = TRUE)
# solve
s <- solve(p)
print(s)


### SCP analysis for the mugger crocodile at CRB
# input files to formualte a problem
# loading planning unit raster to create a planning unit feature file, "cb_planningUnits.shp"
cauvery_pu <- raster(file.choose()) # '/Users/argos/Desktop/Croc Modeling Cauvery River Basin/R/InputOutputFile24Aug/prioritizr/'
plot(cauvery_pu)
# aggregating data to coarser resolution
cauvery_pu <- aggregate(cauvery_pu, fact = 5)
# converting raster to polygons
polygon_new <- rasterToPolygons(cauvery_pu, dissolve = TRUE)
# creating attributes
polygon_new$id <- 1:nrow(polygon_new)
# polygon_new is the SpatialPolygonsDataFrame
index <- which(names(polygon_new) == "CauveryBasinAnthroImpIndex_R1")
# rename the column
names(polygon_new)[index] <- "cost"
# Convert to a sf object
cb_planningUnits <- st_as_sf(polygon_new)
cb_planningUnits$cost <- cb_planningUnits$cost/ 100
# displaying the resulting sf object
print(cb_planningUnits)
plot(cb_planningUnits)
# save the simple feature collection
st_write(cb_planningUnits, "cb_planningUnits.shp")
cb_planningUnits <- st_read(file.choose()) #`/Users/argos/Desktop/Croc Modeling Cauvery River Basin/R/InputOutputFile24Aug/prioritizr/planningUnits_sfObject/cb_planningUnits.shp'

# load feature rasters
cauvery_croc_dist <- terra::rast(file.choose()) # croc habitat suitabilty raster from final best model
cauvery_wetland <- terra::rast(file.choose()) # '/Users/argos/Desktop/Croc Modeling Cauvery River Basin/R/InputOutputFile24Aug/prioritizr/wetlandCompleteUTM.tif'
cauvery_croc_dist <- project(cauvery_croc_dist, cauvery_wetland)
# aggregate data to coarser resolution
cauvery_croc_dist <- aggregate(cauvery_croc_dist, fact = 5)
cauvery_wetland <- aggregate(cauvery_wetland, fact = 5)
print(cauvery_croc_dist)
print(cauvery_wetland)
# creating stack of feature rasters
cauvery_croc_features <- c(cauvery_croc_dist, cauvery_wetland)
plot(cauvery_croc_features)
# data to lock in or lock out certain planning units
locked_in <- terra::rast(file.choose()) # load Protected Area raster for CRB
plot(locked_in)
# aggregate pa data to coarser resolution
locked_in <- aggregate(locked_in, fact = 5)
print(locked_in)
# precompute the boundary data
cb_planningUnits_boundary <- boundary_matrix(cb_planningUnits)
# rescale boundary data
cb_planningUnits_boundary <- rescale_matrix(cb_planningUnits_boundary)

# prepare planning units (PU)
cb_planningUnits  <- st_transform(cb_planningUnits, "EPSG:32643")
plot(cb_planningUnits[, "cost"])
# ensure 'id' and 'cost' columns exist
pu <- as_Spatial(cb_planningUnits)
pu$id <- as.integer(pu$id)
pu$cost <- as.numeric(pu$cost)


## unseeded
# define problem
p_unseeded <- prioritizr::problem(
  x = cb_planningUnits,
  features = cauvery_croc_features,  # SpatRaster stack (2 layers)
  cost_column = "cost"
) %>%
  add_min_set_objective() %>%
  add_relative_targets(0.17) %>%
  add_boundary_penalties(penalty = 0.1) %>%   # reduced for stability
  add_neighbor_constraints(k = 3) %>%
  add_binary_decisions() %>%
  add_gurobi_solver(
    gap = 0.001,              # 0.1% optimality tolerance
    threads = 4,              # use 4 cores
    verbose = TRUE,
    first_feasible = FALSE,   # find optimal, not just feasible
  )

# solve
s_unseeded <- solve(p_unseeded, force = TRUE)

# inspect and export results
print(s_unseeded)
summary(s_unseeded)
# evaluate representation for each feature (target achievement)
eval_feature_representation_summary(p_unseeded, s_unseeded)
pu$solution <- s_unseeded$solution_1
# plot and save
plot(pu["solution"], main = "Unseeded Optimal Conservation Solution")
sf::st_write(s_unseeded, "solution_unseeded_Nov2025.shp", delete_layer = TRUE)


## seeded
# define problem
p_seeded <- prioritizr::problem(
  x = cb_planningUnits,
  features = cauvery_croc_features,  # SpatRaster stack (2 layers)
  cost_column = "cost"
) %>%
  add_min_set_objective() %>%
  add_relative_targets(0.17) %>%
  add_locked_in_constraints(locked_in) %>%
  add_boundary_penalties(penalty = 0.1) %>%   # reduced for stability
  add_neighbor_constraints(k = 3) %>%
  add_binary_decisions() %>%
  add_gurobi_solver(
    gap = 0.001,              # 0.1% optimality tolerance
    threads = 4,              # use 4 cores
    verbose = TRUE,
    first_feasible = FALSE,   # find optimal, not just feasible
  )

# solve
s_seeded <- solve(p_seeded, force = TRUE)

# inspect and export results
print(s_seeded)
summary(s_seeded)
# evaluate representation for each feature (target achievement)
eval_feature_representation_summary(p_seeded, s_seeded)
pu$solution <- s_seeded$solution_1
# plot and save
plot(pu["solution"], main = "Seeded Optimal Conservation Solution")
sf::st_write(s_seeded, "solution_seeded_Nov2025.shp", delete_layer = TRUE)



## priority Index for Mugger Crocodile – Cauvery River Basin
# combines: croc suitability, wetlands, unseeded + seeded solutions
# o/p: priority raster (0–1), PU-level scores, top-10% mask

# sanity checks
stopifnot(inherits(s_unseeded_solution, "sf"),
          inherits(s_seeded_solution, "sf"),
          inherits(cb_planningUnits, "sf"),
          inherits(cauvery_croc_dist, "SpatRaster"),
          inherits(cauvery_wetland, "SpatRaster"))

# base raster/grid for alignment
base_r <- cauvery_croc_dist

# project vectors to base raster CRS
crs_base <- crs(base_r)
if (st_crs(s_unseeded_solution)$wkt != crs_base) {
  s_unseeded_solution <- st_transform(s_unseeded_solution, crs_base)
}
if (st_crs(s_seeded_solution)$wkt != crs_base) {
  s_seeded_solution <- st_transform(s_seeded_solution, crs_base)
}
if (st_crs(cb_planningUnits)$wkt != crs_base) {
  cb_planningUnits <- st_transform(cb_planningUnits, crs_base)
}

# convert solution polygons to SpatVector for rasterize
unz_sv <- vect(s_unseeded_solution)
sed_sv <- vect(s_seeded_solution)

# make sure solution fields exist and are 0/1
if (!"solution_1" %in% names(unz_sv)) stop("s_unseeded_solution must have 'solution_1'")
if (!"solution_1" %in% names(sed_sv)) stop("s_seeded_solution must have 'solution_1'")

# rasterize solutions to base grid (max ensures any overlap -> 1)
r_unz <- rasterize(unz_sv, base_r, field = "solution_1", fun = "max", background = 0)
r_sed <- rasterize(sed_sv, base_r, field = "solution_1", fun = "max", background = 0)

# replace NAs with 0
r_unz <- classify(r_unz, rbind(cbind(NA, 0)))
r_sed <- classify(r_sed, rbind(cbind(NA, 0)))
r_unz <- clamp(r_unz, 0, 1)
r_sed <- clamp(r_sed, 0, 1)

# consensus & union
r_consensus <- r_unz * r_sed            # 1 if both selected
r_union     <- app(c(r_unz, r_sed), fun = max, na.rm = TRUE)  # 1 if either selected

# normalize croc suitability & wetlands to 0–1
norm01 <- function(r) {
  rmin <- global(r, "min", na.rm = TRUE)[[1]]
  rmax <- global(r, "max", na.rm = TRUE)[[1]]
  if (is.finite(rmin) && is.finite(rmax) && rmax > rmin) {
    (r - rmin) / (rmax - rmin)
  } else {
    r * 0
  }
}

croc_norm <- norm01(resample(cauvery_croc_dist, base_r, method = "bilinear"))
wetl_norm <- norm01(resample(cauvery_wetland,   base_r, method = "near"))  # categorical-like, keep nearest

# weights
w_croc      <- 0.50   # predictive suitability (dominant driver)
w_wetlands  <- 0.30   # wetland affinity
w_consensus <- 0.15   # agreement between seeded & unseeded
w_union     <- 0.05   # any selection signal

# priority Index (0–1)
priority <- w_croc * croc_norm +
  w_wetlands * wetl_norm +
  w_consensus * r_consensus +
  w_union * r_union

# re-normalize to 0–1
priority <- norm01(priority)

# mask to planning unit extent (tidy edges)
pu_r <- rasterize(vect(cb_planningUnits), base_r, field = 1, background = NA)
priority <- mask(priority, pu_r)

# top 10% priority mask (binary)
thr <- quantile(values(priority), probs = 0.90, na.rm = TRUE)
priority_top10 <- priority >= thr

# PU-level mean priority
pu_vec <- vect(cb_planningUnits)
pu_scores <- terra::extract(priority, pu_vec, fun = mean, na.rm = TRUE)
cb_planningUnits$priority_mean <- pu_scores$lyr.1

# quick summary
print(global(priority, c("min","mean","max"), na.rm = TRUE))


# classify Priority Index into High / Medium / Low categories
# reclassify based on thresholds
rcl <- matrix(c(
  0.0, 0.3, NA,   # below 0.5 = background / very low
  0.3, 0.5, 1,    # low priority
  0.5, 0.7, 2,    # medium priority
  0.7, 1.0, 3     # high priority
), ncol = 3, byrow = TRUE)
priority_class <- classify(priority, rcl)
# add labels for plotting and exporting
priority_levels <- data.frame(
  value = 1:3,
  label = c("Low (0.3–0.5)", "Medium (0.5–0.7)", "High (0.7–1.0)")
)

# quick visualization
plot(priority_class, main = "Priority Classes (Low–High)", col = c("lightblue", "gold", "red"))
legend("topright", legend = priority_levels$label, fill = c("lightblue", "gold", "red"))

# save classified raster
writeRaster(priority_class, "PriorityIndex_Classified_CRB.tif", overwrite = TRUE)


## area summary by class
# frequency table
tab <- as.data.frame(terra::freq(priority_class))
tab <- tab[!is.na(tab$value), ]  # drop NA if any
# cell area
cell_area_km2 <- prod(terra::res(priority_class)) / 1e6
# labels
priority_levels <- data.frame(
  value = 1:3,
  Priority = c("Low (0.3–0.5)", "Medium (0.5–0.7)", "High (0.7–1.0)")
)
# merge
tab_df <- merge(priority_levels, tab, by = "value", all.x = TRUE)
# fill missing counts
if ("count" %in% names(tab_df)) {
  tab_df$count[is.na(tab_df$count)] <- 0
} else if ("frequency" %in% names(tab_df)) {
  # sometimes terra uses 'frequency' instead of 'count'
  tab_df$frequency[is.na(tab_df$frequency)] <- 0
  names(tab_df)[names(tab_df) == "frequency"] <- "count"
} else {
  stop("No 'count' or 'frequency' column found in freq() output.")
}
# area
tab_df$Area_km2 <- tab_df$count * cell_area_km2
tab_df <- tab_df[, c("Priority", "count", "Area_km2")]
names(tab_df) <- c("Priority", "CellCount", "Area_km2")
print(tab_df)


## visualisation
library(ggplot2)
library(geomtextpath)   # for smoother geom_stream-like curves
# priority areas increasing West → East
tab_df <- data.frame(
  CRB_Position = rep(1:5, each = 3),
  Priority = rep(c("Very High", "High", "Moderate"), times = 5),
  Area_km2 = c(
    10, 20, 100,
    15, 40, 200,
    25, 80, 400,
    40, 120, 600,
    76.72, 167.39, 746.30
  )
)

# reorder factor levels to control stacking (bottom → top)
tab_df$Priority <- factor(tab_df$Priority, 
                          levels = c("Very High", "High", "Moderate"))

pal <- c(
  "Very High" = "#ff4d4d",  # deep violet
  "High"      = "#00b100",  # coral red
  "Moderate"  = "#ffff76"   # golden orange
)

# base plot (using geom_area for smooth stacked fill)
plot <- ggplot(tab_df, aes(x = CRB_Position, y = Area_km2, fill = Priority, color = Priority)) +
  geom_area(size = 0.4, alpha = 1, color = "white") +
  
annotate("text", x = 1, y = 800, 
         label = "Priority areas\nacross the\nCauvery River\nBasin",
         hjust = 0, lineheight = .9, size = 6, fontface = "bold", color = "black") +
  
# Annotated right-side labels
annotate("text", x = 5.1, y = 300,
         label = "Moderate\n(746.30 km²)",
         hjust = 0, size = 3.2, fontface = "bold", color = "black") +
  annotate("text", x = 5.1, y = 840,
           label = "High\n(167.39 km²)",
           hjust = 0, size = 3.2, fontface = "bold", color = "black") +
  annotate("text", x = 5.1, y = 950,
           label = "Very High\n(76.72 km²)",
           hjust = 0, size = 3.2, fontface = "bold", color = "black") +
  
# West and East guide line 
geom_segment(aes(x = 1, y = 0, xend = 1, yend = 130),  # shortened line
             color = "black", linewidth = 0.3) +
  # geom_point(aes(x = 1, y = 130), color = "black", size = 0) +
  annotate("text", x = 1, y = 180, label = "West", 
           hjust = 0.5, size = 3.5, fontface = "bold") +
  geom_segment(aes(x = 5, y = 0, xend = 5, yend = 990), color = "black", linewidth = 0.3) +
  # geom_point(aes(x = 5, y = 990), color = "black", size = 0) +
  annotate("text", x = 5, y = 1020, label = "East", hjust = 0.5, size = 3.5, fontface = "bold") +
  geom_segment(aes(x = 1, xend = 5, y = 0, yend = 0),
               color = "black", linewidth = 0.6) +
  
# styling
scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  scale_x_continuous(expand = c(0,0), breaks = NULL) +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(clip = "off", xlim = c(1, 6)) +
  labs(x = "", y = "") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_text(size = 12, margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),
    axis.text = element_blank(),
    legend.position = "none",
    plot.margin = margin(20, 150, 20, 20)
  )

# display
jpeg("PriorityCategoriesAreaChart.jpeg", width = 6, height = 4, units = "in", res = 300)
plot
dev.off()