library(terra)
library(randomForest)
library(dplyr)
library(sf)
library(raster)
library(caTools)
library(corrplot)
library(e1071)
library(RStoolbox)
library(rsample)
library(pROC)
library(PRROC)

# load rasters
nlcd <- rast("NLCD_masked_raster.tif")
dts <- rast("Buncombe_DistanceToStream_Reproject.tif")
twi <- rast("Buncombe_TWI_Reproject.tif")
ndvi <- rast("NDVI_Sep_2024 1.tif")
svi <- rast("Buncombe_SocialVuln_Reproject.tif")
buffer <- rast("Buffer_Raster.tif")
flood <- rast("FloodHazard_Raster.tif")

# reproject NLCD
nlcd <- project(nlcd, ndvi)

# check coordinate systems and resolution
st_crs(twi) == st_crs(dts)
crs(nlcd) == crs(ndvi)
st_crs(flood) == st_crs(twi)

# list of rasters
rasters <- list(nlcd, dts, twi, ndvi, svi, buffer, flood)

# initialize common extent with the extent of the first raster
common_extent <- ext(ndvi)

# iteratively find the intersection of all extents
for (i in 2:length(rasters)) {
  common_extent <- terra::intersect(common_extent, ext(rasters[[i]]))
}

for (raster in rasters){
  cropped_raster <- crop(raster, common_extent)
}

# manually crop
nlcd_cropped <- crop(nlcd, common_extent)
dts_cropped <- crop(dts, common_extent)
twi_cropped <- crop(twi, common_extent)
ndvi_cropped <- crop(ndvi, common_extent)
svi_cropped <- crop(svi, common_extent)
buffer_cropped <- crop(buffer, common_extent)
flood_cropped <- crop(flood, common_extent)

# resample
nlcd_resampled <- resample(nlcd_cropped, ndvi_cropped, method = "near")
dts_resampled <- resample(dts_cropped, ndvi_cropped, method = "near")
twi_resampled <- resample(twi_cropped, ndvi_cropped, method = "near")
svi_resampled <- resample(svi_cropped, ndvi_cropped, method = "near")
buffer_resampled <- resample(buffer_cropped, ndvi_cropped, method = "near")
flood_resampled <- resample(flood_cropped, ndvi_cropped, method = "near")

# create raster stack
raster_stack <- c(nlcd_resampled, dts_resampled, twi_resampled, ndvi_cropped,
                  svi_resampled, buffer_resampled, flood_resampled)
names(raster_stack) <- c("NLCD", "DTS", "TWI", "NDVI", "SVI", "Buffer", "Flood")

print(raster_stack)
plot(raster_stack, stretch = "hist")

# load in county shapefile to crop rasters
buncombe_county <- st_read("BucombeCounty.shp")
plot(buncombe_county)
buncombe_county <- st_transform(buncombe_county, crs(ndvi))

# turn raster stack into a data frame
flood_df <- as.data.frame(raster_stack, xy = TRUE, na.rm = TRUE)
names(flood_df) <- c("x", "y", "NLCD", "DTS",
                     "TWI", "NDVI", "SVI", "Buffer", "Flood")
flood_df$Flood <- as.factor(flood_df$Flood)


yes_flood <- flood_df%>%filter(Flood == 1)
no_flood <- flood_df%>%filter(Flood == 0)

set.seed(123)
no_flood_sample <- no_flood%>%sample_n(nrow(yes_flood))
balanced_data <- bind_rows(yes_flood, no_flood_sample)


# random forest

#remove coordinate columns
predictors <- names(balanced_data)[!names(balanced_data)
                                   %in% c("x", "y", "balanced_data", "Flood")]
rf_test <- as.formula(paste("Flood~", paste(predictors, collapse = "+")))

#subset sample of data
set.seed(123)
sample_size <- floor(0.7 * nrow(balanced_data))
train_indices <- sample(seq_len(nrow(balanced_data)),
                        size = sample_size, replace = FALSE)

train_data <- balanced_data[train_indices, ]
test_data <- balanced_data[-train_indices, ]

#data are too big, so need to sample
sample_rows <- sample(nrow(train_data), 10000)
train_sample <- train_data[sample_rows, ]

#use this if the testing is taking too long
sample_rows <- sample(nrow(test_data), 10000)
test_sample <- test_data[sample_rows, ]

#train
rf_model <- randomForest(rf_test, data = train_sample,
                         ntree = 500, importance = TRUE)

print(rf_model)

# predict classification
rf_predict <- predict(rf_model, newdata = test_sample)
rf_pred_raster <- predict(raster_stack, rf_model, type = "response")
plot(rf_pred_raster, main = "RF Predicted Flood Risk Binary")

# writeRaster(rf_pred_raster,
#             filename = "rf_flood_risk_binary.tif",
#             datatype = "INT1U", overwrite = TRUE)


# predict probabilities
rf_pred_prob <- predict(raster_stack, rf_model, type = "prob")

names(rf_pred_prob)
rf_pred_probability <- rf_pred_prob[[2]]

reclass_m <- matrix(c(
  0.0, 0.3, 1, # low
  0.3, 0.6, 2, # medium
  0.6, 1.0, 3 # high
), ncol = 3, byrow = TRUE)

flood_risk_classification <- classify(rf_pred_probability, reclass_m)
plot(flood_risk_classification)

writeRaster(flood_risk_classification,
            filename = "rf_flood_risk_classification.tif",
            datatype = "INT1U", overwrite = TRUE)



#use this confusionMatrix(rf_predict, test_data$Flood)
varImpPlot(rf_model)

# -------------------------------------------------------------------------
# svm model
svm_model <- svm(rf_test, data = train_sample,
                 kernel = "radial", cost = 10, gamma = 0.1)

#predict on test set
svm_predict <- predict(svm_model, newdata = test_sample)

#confusion matrix
svm_actual <- test_sample$Flood
svm_predicted <- factor(svm_predict, levels = levels(svm_actual))
svm_conf_matrix <- table(Predicted = svm_predicted, Actual = svm_actual)
print(svm_conf_matrix)

#accuracy
svm_accuracy <- mean(svm_predicted == svm_actual)
cat("SVM Accuracy:", svm_accuracy, "\n")

#hyperparameter tuning
tuned_svm <- tune(svm, rf_test, data = train_sample,
                  kernel = "radial",
                  ranges = list(cost = c(0.1, 1, 10),
                                gamma = c(0.01, 0.1, 1)))

summary(tuned_svm)

#best model
svm_best <- tuned_svm$best.model

#predict using best model
svm_best_pred <- predict(svm_best, newdata = test_sample)
svm_best_accuracy <- mean(svm_best_pred == svm_actual)
cat("Tuned SVM Accuracy:", svm_best_accuracy, "\n")

# -----------------------------------------

# Train SVM with probability output enabled
svm_model <- svm(rf_test, data = train_sample,
                 kernel = "radial", cost = 10, gamma = 0.1,
                 probability = TRUE)

# Predict on test set with probabilities
svm_predict <- predict(svm_model, newdata = test_sample, probability = TRUE)
svm_probs <- attr(svm_predict, "probabilities")[, "1"]

# Confusion matrix
svm_actual <- test_sample$Flood
svm_predicted <- factor(svm_predict, levels = levels(svm_actual))
svm_conf_matrix <- table(Predicted = svm_predicted, Actual = svm_actual)
print(svm_conf_matrix)

# Accuracy
svm_accuracy <- mean(svm_predicted == svm_actual)
cat("SVM Accuracy:", svm_accuracy, "\n")

# ROC Curve
flood_true <- as.numeric(as.character(test_sample$Flood))
roc_curve_svm <- roc(flood_true, svm_probs)
plot(roc_curve_svm, col = "darkred", main = "SVM ROC Curve", print.auc = TRUE)

# Predict raster layer using SVM (probability)
svm_pred_stack <- predict(raster_stack, svm_model, type = "response", na.rm=TRUE)

#Custom function to return probability of flood (class 1)
svm_prob_fun <- function(model, data) {
  p <- predict(model, data, probability = TRUE)
  attr(p, "probabilities")[, "1"]  # Return only class 1 (flood) probability
}

#Apply to raster stack using terra::predict with custom function
svm_pred_probability <- predict(raster_stack, svm_model, fun = svm_prob_fun, na.rm = TRUE)

# Reclassify probability raster into risk categories
reclass_m <- matrix(c(
  0.0, 0.3, 1,
  0.3, 0.6, 2,
  0.6, 1.0, 3
), ncol = 3, byrow = TRUE)
svm_risk_classification <- classify(svm_pred_probability, reclass_m)

# Plot maps
plot(svm_pred_stack, main = "SVM Predicted Flood Risk")
plot(svm_pred_probability, main = "SVM Predicted Flood Probability")
plot(svm_risk_classification, main = "SVM Flood Risk Classification")

# Save classified raster 
writeRaster(svm_pred_stack,
            filename = "svm_flood_risk_binary.tif",
            datatype = "INT1U", overwrite = TRUE)
