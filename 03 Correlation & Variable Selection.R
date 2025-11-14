## Correlation and variable selection

library(corrplot)
library(caret)
library(mecofun)
library(ggplot2)
library(dplyr) 

## Load and prepare data
data <- read.csv("/Users/argos/Desktop/Croc Modeling Cauvery River Basin/R/Revision3analysis/BiasWeighted_BGpoints_5k_env_latestwithspThin.CSV")
data <- data %>% dplyr::select(-Longitude, -Latitude, -presence)
num_data <- data %>% dplyr::select(where(is.numeric))

full_data <- read.csv("/Users/argos/Desktop/Croc Modeling Cauvery River Basin/R/Revision3analysis/BiasWeighted_BGpoints_5k_env_latestwithspThin.CSV")
response <- as.numeric(full_data$presence)

num_data <- data.frame(lapply(num_data, as.numeric))
num_data <- num_data[, colSums(is.na(num_data)) < nrow(num_data)]

clean_data <- cbind(num_data, response)
clean_data <- na.omit(clean_data)

response <- clean_data$response
num_data <- clean_data[, !(names(clean_data) == "response")]

num_data <- num_data[, sapply(num_data, function(x) length(unique(na.omit(x))) > 2)]

removed_vars <- setdiff(names(data), names(num_data))
cat("Dropped variables with <3 unique values:\n")
print(removed_vars)

## Run select07
sel07 <- mecofun::select07(X = num_data, y = response, threshold = 0.7)
selected_vars <- sel07$pred_sel
print(selected_vars)

var_status <- data.frame(
  Variable = colnames(num_data),
  Status = ifelse(colnames(num_data) %in% selected_vars, "Retained", "Removed")
)
write.csv(var_status, "TableS2_select07_status.csv", row.names = FALSE)

## Correlation heatmap (|r| < 0.7)

library(dplyr)
library(corrplot)

selected_vars_final <- c(
  "FlowAcc_SD_4000", "FlowAcc_FM_8000", "FlowAcc_SD_2000",
  "LULC_Base", "FlowAcc_SD_1000", "TWI_SD_500", "SolarRad_SD_1000",
  "FlowAcc_Base", "DistWater_Base", "Wetland_GYRATE_AM_2000",
  "Wetland_GYRATE_AM_4000", "CC_Edge_8000", "AI_SD_1000",
  "TWI_FM_8000", "Wetland_GYRATE_AM_1000", "TWI_SD_2000",
  "LULC_CWED_1000", "DistRoad_FM_8000", "DistWater_SD_8000",
  "bio9_FM_2000", "TWI_Base", "bio3_Base", "CC_Edge_2000",
  "DistRoad_Base", "AI_Base", "SolarRad_SD_2000", "DistWater_SD_1000",
  "bio4_FM_2000", "bio18_FM_8000", "TWI_SD_8000", "bio4_FM_4000",
  "bio7_FM_8000", "TWI_SD_4000", "bio15_FM_8000", "CC_Base",
  "HydroConditionedDEM_Base", "Wetland_Base", "Wetland_GYRATE_AM_500"
)

corr_data <- dplyr::select(data, dplyr::any_of(selected_vars_final))
corr_data[] <- lapply(corr_data, as.numeric)
cor_matrix <- cor(corr_data, use = "pairwise.complete.obs", method = "pearson")
cor_matrix[is.na(cor_matrix)] <- 0

cor_matrix_filtered <- cor_matrix
cor_matrix_filtered[abs(cor_matrix) >= 0.7] <- NA

jpeg("Supp_Fig_CorrelationMatrix_LowCorr.jpeg",
     width = 10, height = 8, units = "in", res = 600)

corrplot(
  cor_matrix_filtered,
  method = "color",
  type = "upper",
  order = "original",
  tl.cex = 0.7,
  tl.col = "grey20",
  tl.srt = 60,
  na.label = " ",
  addCoef.col = "grey20",
  number.cex = 0.45,
  col = colorRampPalette(c("#a4c8ed", "#F7F7F7", "#db4241"))(200),
  mar = c(0, 0, 1, 0)
)

dev.off()