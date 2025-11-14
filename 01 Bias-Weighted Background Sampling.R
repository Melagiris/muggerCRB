## Bias-weighted background sampling (wetlands, roads, settlements)
library(terra)
library(sf)

## Load input rasters and basin boundary
wet <- rast("wetLandCompleteUTM.tif")
dist_road <- rast("DistRoad_Base.tif")
dist_settlement <- rast("DistSettlement_Base.tif")
boundary_v <- vect("CauevryBasinBoundary.shp")

## Project boundary and distance layers to match wetland CRS
boundary_v_utm <- project(boundary_v, crs(wet))
dist_road_utm <- project(dist_road, wet)
dist_settle_utm <- project(dist_settlement, wet)

## Smooth wetland raster and mask to basin extent
wet_fill <- classify(wet, cbind(NA, 0))
wet_smooth <- focal(wet_fill, w = 23, fun = mean, na.policy = "omit")
wet_smooth <- crop(wet_smooth, boundary_v_utm)
wet_smooth <- mask(wet_smooth, boundary_v_utm)

## Road and settlement bias layers (reversed distance)
bias_road <- 1 - (dist_road_utm / global(dist_road_utm, "max", na.rm = TRUE)[[1]])
bias_settle <- 1 - (dist_settle_utm / global(dist_settle_utm, "max", na.rm = TRUE)[[1]])
bias_road_resamp <- resample(bias_road, wet_smooth, method = "bilinear")
bias_settle_resamp <- resample(bias_settle, wet_smooth, method = "bilinear")
bias_road_resamp <- mask(bias_road_resamp, boundary_v_utm)
bias_settle_resamp <- mask(bias_settle_resamp, boundary_v_utm)

## Combine layers (wetlands dominate bias) and normalize to 0–1
bias_combined <- (0.6 * wet_smooth) + (0.2 * bias_road_resamp) + (0.2 * bias_settle_resamp)
bias_min <- global(bias_combined, "min", na.rm = TRUE)[[1]]
bias_max <- global(bias_combined, "max", na.rm = TRUE)[[1]]
bias_combined <- (bias_combined - bias_min) / (bias_max - bias_min)

## Final mask and trimming
bias_final <- mask(bias_combined, boundary_v_utm)
bias_final <- trim(bias_final)
valid_mask <- !is.na(bias_final)
bias_final_masked <- mask(bias_final, valid_mask)
writeRaster(bias_final_masked, "bias_final.tif", overwrite = TRUE)

## Sample 5000 bias-weighted background points
set.seed(123)
bg_points <- spatSample(
  bias_final_masked,
  size = 5000,
  method = "random",
  prob = TRUE,
  as.points = TRUE,
  na.rm = TRUE,
  ext = ext(boundary_v_utm)
)

## Verify all points fall within basin
inside_check <- relate(bg_points, boundary_v_utm, "intersects")
cat("Points outside boundary:", sum(!inside_check), "\n")

## Export background points
writeVector(bg_points, "bias_weighted_background_points_5k.shp", overwrite = TRUE)
write.csv(bg_points, "bias_bg_5k_latlong.csv", row.names = FALSE)


## Extract environmental values for 5k points
library(raster)

## Load bias raster and points
bias_final <- rast("bias_final.tif")
bg_points <- vect("bias_weighted_background_points_5k.shp")

## Load environmental stack
directory_mean <- "Rasters_Final/finalContinousAndCategorical"
files_mean <- list.files(directory_mean, pattern = "\\.tif$", full.names = TRUE)
env_stack <- rast(files_mean)

## Extract environmental values
bg_env_data <- extract(env_stack, bg_points)

## Add geographic coordinates (WGS84)
bg_points_wgs <- project(bg_points, "EPSG:4326")
coords_wgs <- crds(bg_points_wgs, df = TRUE)
colnames(coords_wgs) <- c("Longitude", "Latitude")
bg_env_data <- cbind(coords_wgs, bg_env_data)
bg_env_data$Presence <- 0

## Save dataset
write.csv(bg_env_data, "BiasWeighted_BGpoints_5k_env.csv", row.names = FALSE)


## Generate 10k bias-weighted background points
valid_mask <- !is.na(bias_final)
set.seed(124)
bg_points_10k <- spatSample(
  valid_mask,
  size = 10000,
  method = "random",
  as.points = TRUE,
  na.rm = TRUE
)

## Extract environmental values
bg_env_data_10k <- extract(env_stack, bg_points_10k)

## Add WGS84 coordinates
bg_points_10k_wgs <- project(bg_points_10k, "EPSG:4326")
coords_10k_wgs <- crds(bg_points_10k_wgs, df = TRUE)
colnames(coords_10k_wgs) <- c("Longitude", "Latitude")
bg_env_data_10k <- cbind(coords_10k_wgs, bg_env_data_10k)
bg_env_data_10k$Presence <- 0

## Remove rows containing NA
cat("Total rows before cleaning:", nrow(bg_env_data_10k), "\n")
bg_env_data_10k_clean <- bg_env_data_10k[complete.cases(bg_env_data_10k), ]
cat("Number of complete records remaining:", nrow(bg_env_data_10k_clean), "\n")

## Save both raw and cleaned datasets
write.csv(bg_env_data_10k_clean, "BiasWeighted_BGpoints_10k_env_CLEAN_2.csv", row.names = FALSE)
write.csv(bg_env_data_10k, "BiasWeighted_BGpoints_10k_env.csv", row.names = FALSE)

## Column-wise NA count
na_summary <- colSums(is.na(bg_env_data_10k))
print(na_summary)