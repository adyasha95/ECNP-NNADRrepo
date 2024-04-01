setwd('/volume/projects/JF_ECNP/Analysis/Analysis/R_models_150224/Regression')
load('/volume/mitnvp1_scratch/AK_Testing/ECNP/regression_age_allsites_all_atlas_meanoffset_UniDistr_normHC_withgroup_offset.RData')

res$aggr
age_differences <- numeric()

# Initialize vectors to store the average predictions and original ages
average_predictions <- numeric()
original_ages <- numeric()

# Loop through each individual
for (individual_id in unique(res$pred$data$id)) {
  # Extract the predictions for the current individual
  predictions <- res$pred$data$response[res$pred$data$id == individual_id]
  
  # Extract the original age for the current individual
  original_age <- res$pred$data$truth[res$pred$data$id == individual_id][1]
  
  # Calculate the average prediction for the current individual
  avg_prediction <- mean(predictions)
  
  # Append the average prediction and original age to their respective vectors
  average_predictions <- c(average_predictions, avg_prediction)
  original_ages <- c(original_ages, original_age)
}

# Calculate the age differences outside the loop
age_differences <- average_predictions - original_ages
mean_age_difference<-mean(age_differences)
mean_age_difference
sd_age_difference<-sd(age_differences)
sd_age_difference
correlation <- cor(average_predictions, original_ages)
correlation
data_norm<-data.frame(average_predictions,original_ages)
write.csv(data_norm, "allatlas_HCnorm_new.csv", row.names = FALSE)

model <- lm(original_ages ~ average_predictions)

# Obtain summary of the model
summary_model <- summary(model)

# Extract the R-squared value
rsquared <- summary_model$r.squared
# Print the R-squared value
print(rsquared)
absolute_differences <- abs(average_predictions - original_ages)

# Calculate the Mean Absolute Error (MAE)
mae <- mean(absolute_differences)
mae


# application to HC rem
load('/volume/mitnvp1_scratch/AK_Testing/ECNP/regression_age_all_atlas_meanoffset_UniDistr_appliedremHC_withgroup_offset.RData')
# Extract the predicted and actual values

# Extract the predicted and actual values
# Get the predicted values from the resample results

average_predictions_HC_rem <- rowMeans(sapply(predictions_patients_regr_rem, function(pred) pred$data$response))
original_age_HC_rem <- rowMeans(sapply(predictions_patients_regr_rem, function(pred) pred$data$truth))
data_HCrem_new<-data.frame(average_predictions_HC_rem,original_age_HC_rem)
write.csv(data_HCrem_new, "allatlas_HCrem_new.csv", row.names = FALSE)

model <- lm(original_age_HC_rem ~ average_predictions_HC_rem)

# Obtain summary of the model
summary_model <- summary(model)

# Extract the R-squared value
rsquared <- summary_model$r.squared


# Print the R-squared value
print(rsquared)
correlation <- cor(average_predictions_HC_rem, original_age_HC_rem)
correlation
# Calculate the differences between predicted age and original age for the patient group
age_difference_patients <- average_predictions_HC_rem - original_age_HC_rem

# Calculate the mean of the differences for the patient group
mean_age_difference_patients <- mean(age_difference_patients)
mean_age_difference_patients
sd_age_difference_patients<-sd(age_difference_patients)
sd_age_difference_patients
# Calculate the absolute differences
absolute_differences <- abs(average_predictions_HC_rem - original_age_HC_rem)

# Calculate the Mean Absolute Error (MAE)
mae <- mean(absolute_differences)
mae

# Display the mean_age_difference for the patient group
cat("Mean of the differences between predicted age and original age:", mean_age_difference)
cat("Mean of the differences between predicted age and original age for patients:", mean_age_difference_patients)




# application to SCZ all
load('/volume/mitnvp1_scratch/AK_Testing/ECNP/regression_age_all_atlas_meanoffset_UniDistr_appliedsczall_withgroup_offset.RData')
# Extract the predicted and actual values

# Extract the predicted and actual values
# Get the predicted values from the resample results

average_predictions_patients_scz <- rowMeans(sapply(predictions_patients_scz, function(pred) pred$data$response))
original_age_scz <- rowMeans(sapply(predictions_patients_scz, function(pred) pred$data$truth))
data_scz<-data.frame(average_predictions_patients_scz,original_age_scz)
write.csv(data_scz, "allatlas_scz_new.csv", row.names = FALSE)

model <- lm(original_age_scz ~ average_predictions_patients_scz)

# Obtain summary of the model
summary_model <- summary(model)

# Extract the R-squared value
rsquared <- summary_model$r.squared

# Print the R-squared value
print(rsquared)
correlation <- cor(average_predictions_patients_scz, original_age_scz)
correlation
# Calculate the differences between predicted age and original age for the patient group
age_difference_patients <- average_predictions_patients_scz - original_age_scz

# Calculate the mean of the differences for the patient group
mean_age_difference_patients <- mean(age_difference_patients)
mean_age_difference_patients
sd_age_difference_patients<-sd(age_difference_patients)
sd_age_difference_patients
# Calculate the absolute differences
absolute_differences <- abs(average_predictions_patients_scz - original_age_scz)

# Calculate the Mean Absolute Error (MAE)
mae <- mean(absolute_differences)
mae

# Display the mean_age_difference for the patient group
cat("Mean of the differences between predicted age and original age:", mean_age_difference)
cat("Mean of the differences between predicted age and original age for patients:", mean_age_difference_patients)

plot(original_age_scz, average_predictions_patients_scz, xlab = "Actual age", ylab = "Predicted age", main = "Predicted vs Actual Age")
abline(a = 0, b = 1, col = "red")



## schaefers
load('/volume/mitnvp1_scratch/AK_Testing/ECNP/regression_age_allsites_schaefer200_gmv_meanoffset_UniDistr_normHC_withgroup_offset.RData')

res$aggr
age_differences <- numeric()

# Initialize vectors to store the average predictions and original ages
average_predictions <- numeric()
original_ages <- numeric()

# Loop through each individual
for (individual_id in unique(res$pred$data$id)) {
  # Extract the predictions for the current individual
  predictions <- res$pred$data$response[res$pred$data$id == individual_id]
  
  # Extract the original age for the current individual
  original_age <- res$pred$data$truth[res$pred$data$id == individual_id][1]
  
  # Calculate the average prediction for the current individual
  avg_prediction <- mean(predictions)
  
  # Append the average prediction and original age to their respective vectors
  average_predictions <- c(average_predictions, avg_prediction)
  original_ages <- c(original_ages, original_age)
}

# Calculate the age differences outside the loop
age_differences <- average_predictions - original_ages
mean_age_difference<-mean(age_differences)
mean_age_difference
sd_age_difference<-sd(age_differences)
sd_age_difference
correlation <- cor(average_predictions, original_ages)
correlation
data_norm<-data.frame(average_predictions,original_ages)
write.csv(data_norm, "schaefer200_gmv_HCnorm_new.csv", row.names = FALSE)
model <- lm(original_ages ~ average_predictions)

# Obtain summary of the model
summary_model <- summary(model)

# Extract the R-squared value
rsquared <- summary_model$r.squared
# Print the R-squared value
print(rsquared)
absolute_differences <- abs(average_predictions - original_ages)

# Calculate the Mean Absolute Error (MAE)
mae <- mean(absolute_differences)
mae

# application to HC rem
load('/volume/mitnvp1_scratch/AK_Testing/ECNP/regression_age_schaefer200_gmv_meanoffset_UniDistr_appliedremHC_withgroup_offset.RData')
# Extract the predicted and actual values

# Extract the predicted and actual values
# Get the predicted values from the resample results

average_predictions_HC_rem <- rowMeans(sapply(predictions_patients_regr_rem, function(pred) pred$data$response))
original_age_HC_rem <- rowMeans(sapply(predictions_patients_regr_rem, function(pred) pred$data$truth))
data_HCrem_new<-data.frame(average_predictions_HC_rem,original_age_HC_rem)
write.csv(data_HCrem_new, "schaefer200_gmv_HCrem_new.csv", row.names = FALSE)

model <- lm(original_age_HC_rem ~ average_predictions_HC_rem)

# Obtain summary of the model
summary_model <- summary(model)

# Extract the R-squared value
rsquared <- summary_model$r.squared

# Print the R-squared value
print(rsquared)
correlation <- cor(average_predictions_HC_rem, original_age_HC_rem)
correlation
# Calculate the differences between predicted age and original age for the patient group
age_difference_patients <- average_predictions_HC_rem - original_age_HC_rem

# Calculate the mean of the differences for the patient group
mean_age_difference_patients <- mean(age_difference_patients)
mean_age_difference_patients
sd_age_difference_patients<-sd(age_difference_patients)
sd_age_difference_patients
# Calculate the absolute differences
absolute_differences <- abs(average_predictions_HC_rem - original_age_HC_rem)

# Calculate the Mean Absolute Error (MAE)
mae <- mean(absolute_differences)
mae

##SCZ
load('/volume/mitnvp1_scratch/AK_Testing/ECNP/regression_age_schaefer200_gmv_meanoffset_UniDistr_appliedsczall_withgroup_offset.RData')
# Extract the predicted and actual values

# Extract the predicted and actual values
# Get the predicted values from the resample results

average_predictions_patients_scz <- rowMeans(sapply(predictions_patients_scz, function(pred) pred$data$response))
original_age_scz <- rowMeans(sapply(predictions_patients_scz, function(pred) pred$data$truth))
data_scz<-data.frame(average_predictions_patients_scz,original_age_scz)
write.csv(data_scz, "schaefer200_gmv_scz_new.csv", row.names = FALSE)

model <- lm(original_age_scz ~ average_predictions_patients_scz)

# Obtain summary of the model
summary_model <- summary(model)

# Extract the R-squared value
rsquared <- summary_model$r.squared

# Print the R-squared value
print(rsquared)
correlation <- cor(average_predictions_patients_scz, original_age_scz)
correlation
# Calculate the differences between predicted age and original age for the patient group
age_difference_patients <- average_predictions_patients_scz - original_age_scz

# Calculate the mean of the differences for the patient group
mean_age_difference_patients <- mean(age_difference_patients)
mean_age_difference_patients
sd_age_difference_patients<-sd(age_difference_patients)
sd_age_difference_patients
# Calculate the absolute differences
absolute_differences <- abs(average_predictions_patients_scz - original_age_scz)

# Calculate the Mean Absolute Error (MAE)
mae <- mean(absolute_differences)
mae

## hammers
load('/volume/mitnvp1_scratch/AK_Testing/ECNP/regression_age_allsites_hammers_gmv_meanoffset_UniDistr_normHC_withgroup_offset.RData')

res$aggr
age_differences <- numeric()

# Initialize vectors to store the average predictions and original ages
average_predictions <- numeric()
original_ages <- numeric()

# Loop through each individual
for (individual_id in unique(res$pred$data$id)) {
  # Extract the predictions for the current individual
  predictions <- res$pred$data$response[res$pred$data$id == individual_id]
  
  # Extract the original age for the current individual
  original_age <- res$pred$data$truth[res$pred$data$id == individual_id][1]
  
  # Calculate the average prediction for the current individual
  avg_prediction <- mean(predictions)
  
  # Append the average prediction and original age to their respective vectors
  average_predictions <- c(average_predictions, avg_prediction)
  original_ages <- c(original_ages, original_age)
}

# Calculate the age differences outside the loop
age_differences <- average_predictions - original_ages
mean_age_difference<-mean(age_differences)
mean_age_difference
sd_age_difference<-sd(age_differences)
sd_age_difference
correlation <- cor(average_predictions, original_ages)
correlation
data_norm<-data.frame(average_predictions,original_ages)
write.csv(data_norm, "hammers_gmv_HCnorm_new.csv", row.names = FALSE)
model <- lm(original_ages ~ average_predictions)

# Obtain summary of the model
summary_model <- summary(model)

# Extract the R-squared value
rsquared <- summary_model$r.squared
# Print the R-squared value
print(rsquared)
absolute_differences <- abs(average_predictions - original_ages)

# Calculate the Mean Absolute Error (MAE)
mae <- mean(absolute_differences)
mae

# application to HC rem
load('/volume/mitnvp1_scratch/AK_Testing/ECNP/regression_age_hammers_gmv_meanoffset_UniDistr_appliedremHC_withgroup_offset.RData')
# Extract the predicted and actual values

# Extract the predicted and actual values
# Get the predicted values from the resample results

average_predictions_HC_rem <- rowMeans(sapply(predictions_patients_regr_rem, function(pred) pred$data$response))
original_age_HC_rem <- rowMeans(sapply(predictions_patients_regr_rem, function(pred) pred$data$truth))
data_HCrem_new<-data.frame(average_predictions_HC_rem,original_age_HC_rem)
write.csv(data_HCrem_new, "hammers_gmv_HCrem_new.csv", row.names = FALSE)

model <- lm(original_age_HC_rem ~ average_predictions_HC_rem)

# Obtain summary of the model
summary_model <- summary(model)

# Extract the R-squared value
rsquared <- summary_model$r.squared

# Print the R-squared value
print(rsquared)
correlation <- cor(average_predictions_HC_rem, original_age_HC_rem)
correlation
# Calculate the differences between predicted age and original age for the patient group
age_difference_patients <- average_predictions_HC_rem - original_age_HC_rem

# Calculate the mean of the differences for the patient group
mean_age_difference_patients <- mean(age_difference_patients)
mean_age_difference_patients
sd_age_difference_patients<-sd(age_difference_patients)
sd_age_difference_patients
# Calculate the absolute differences
absolute_differences <- abs(average_predictions_HC_rem - original_age_HC_rem)

# Calculate the Mean Absolute Error (MAE)
mae <- mean(absolute_differences)
mae

##SCZ
load('/volume/mitnvp1_scratch/AK_Testing/ECNP/regression_age_hammers_gmv_meanoffset_UniDistr_appliedsczall_withgroup_offset.RData')
# Extract the predicted and actual values

# Extract the predicted and actual values
# Get the predicted values from the resample results

average_predictions_patients_scz <- rowMeans(sapply(predictions_patients_scz, function(pred) pred$data$response))
original_age_scz <- rowMeans(sapply(predictions_patients_scz, function(pred) pred$data$truth))
data_scz<-data.frame(average_predictions_patients_scz,original_age_scz)
write.csv(data_scz, "hammers_gmv_scz_new.csv", row.names = FALSE)

model <- lm(original_age_scz ~ average_predictions_patients_scz)

# Obtain summary of the model
summary_model <- summary(model)

# Extract the R-squared value
rsquared <- summary_model$r.squared

# Print the R-squared value
print(rsquared)
correlation <- cor(average_predictions_patients_scz, original_age_scz)
correlation
# Calculate the differences between predicted age and original age for the patient group
age_difference_patients <- average_predictions_patients_scz - original_age_scz

# Calculate the mean of the differences for the patient group
mean_age_difference_patients <- mean(age_difference_patients)
mean_age_difference_patients
sd_age_difference_patients<-sd(age_difference_patients)
sd_age_difference_patients
# Calculate the absolute differences
absolute_differences <- abs(average_predictions_patients_scz - original_age_scz)

# Calculate the Mean Absolute Error (MAE)
mae <- mean(absolute_differences)
mae

## aal3
load('/volume/mitnvp1_scratch/AK_Testing/ECNP/regression_age_allsites_aal3_meanoffset_UniDistr_normHC_withgroup_offset.RData')

res$aggr
age_differences <- numeric()

# Initialize vectors to store the average predictions and original ages
average_predictions <- numeric()
original_ages <- numeric()

# Loop through each individual
for (individual_id in unique(res$pred$data$id)) {
  # Extract the predictions for the current individual
  predictions <- res$pred$data$response[res$pred$data$id == individual_id]
  
  # Extract the original age for the current individual
  original_age <- res$pred$data$truth[res$pred$data$id == individual_id][1]
  
  # Calculate the average prediction for the current individual
  avg_prediction <- mean(predictions)
  
  # Append the average prediction and original age to their respective vectors
  average_predictions <- c(average_predictions, avg_prediction)
  original_ages <- c(original_ages, original_age)
}

# Calculate the age differences outside the loop
age_differences <- average_predictions - original_ages
mean_age_difference<-mean(age_differences)
mean_age_difference
sd_age_difference<-sd(age_differences)
sd_age_difference
correlation <- cor(average_predictions, original_ages)
correlation
data_norm<-data.frame(average_predictions,original_ages)
write.csv(data_norm, "aal3_HCnorm_new.csv", row.names = FALSE)
model <- lm(original_ages ~ average_predictions)

# Obtain summary of the model
summary_model <- summary(model)

# Extract the R-squared value
rsquared <- summary_model$r.squared
# Print the R-squared value
print(rsquared)
absolute_differences <- abs(average_predictions - original_ages)

# Calculate the Mean Absolute Error (MAE)
mae <- mean(absolute_differences)
mae
# application to HC rem
load('/volume/mitnvp1_scratch/AK_Testing/ECNP/regression_age_aal3_meanoffset_UniDistr_appliedremHC_withgroup_offset.RData')
# Extract the predicted and actual values

# Extract the predicted and actual values
# Get the predicted values from the resample results

average_predictions_HC_rem <- rowMeans(sapply(predictions_patients_regr_rem, function(pred) pred$data$response))
original_age_HC_rem <- rowMeans(sapply(predictions_patients_regr_rem, function(pred) pred$data$truth))
data_HCrem_new<-data.frame(average_predictions_HC_rem,original_age_HC_rem)
write.csv(data_HCrem_new, "aal3_HCrem_new.csv", row.names = FALSE)

model <- lm(original_age_HC_rem ~ average_predictions_HC_rem)

# Obtain summary of the model
summary_model <- summary(model)

# Extract the R-squared value
rsquared <- summary_model$r.squared

# Print the R-squared value
print(rsquared)
correlation <- cor(average_predictions_HC_rem, original_age_HC_rem)
correlation
# Calculate the differences between predicted age and original age for the patient group
age_difference_patients <- average_predictions_HC_rem - original_age_HC_rem

# Calculate the mean of the differences for the patient group
mean_age_difference_patients <- mean(age_difference_patients)
mean_age_difference_patients
sd_age_difference_patients<-sd(age_difference_patients)
sd_age_difference_patients
# Calculate the absolute differences
absolute_differences <- abs(average_predictions_HC_rem - original_age_HC_rem)

# Calculate the Mean Absolute Error (MAE)
mae <- mean(absolute_differences)
mae

##SCZ
load('/volume/mitnvp1_scratch/AK_Testing/ECNP/regression_age_aal3_meanoffset_UniDistr_appliedsczall_withgroup_offset.RData')
# Extract the predicted and actual values

# Extract the predicted and actual values
# Get the predicted values from the resample results

average_predictions_patients_scz <- rowMeans(sapply(predictions_patients_scz, function(pred) pred$data$response))
original_age_scz <- rowMeans(sapply(predictions_patients_scz, function(pred) pred$data$truth))
data_scz<-data.frame(average_predictions_patients_scz,original_age_scz)
write.csv(data_scz, "aal3_scz_new.csv", row.names = FALSE)

model <- lm(original_age_scz ~ average_predictions_patients_scz)

# Obtain summary of the model
summary_model <- summary(model)

# Extract the R-squared value
rsquared <- summary_model$r.squared

# Print the R-squared value
print(rsquared)
correlation <- cor(average_predictions_patients_scz, original_age_scz)
correlation
# Calculate the differences between predicted age and original age for the patient group
age_difference_patients <- average_predictions_patients_scz - original_age_scz

# Calculate the mean of the differences for the patient group
mean_age_difference_patients <- mean(age_difference_patients)
mean_age_difference_patients
sd_age_difference_patients<-sd(age_difference_patients)
sd_age_difference_patients
# Calculate the absolute differences
absolute_differences <- abs(average_predictions_patients_scz - original_age_scz)

# Calculate the Mean Absolute Error (MAE)
mae <- mean(absolute_differences)
mae

