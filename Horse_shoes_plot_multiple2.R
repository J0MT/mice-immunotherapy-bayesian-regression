# Load required libraries
library(rstanarm)
library(tidybayes)
library(dplyr)
library(ggplot2)
library(tidyr)
library(readxl)
library(stringr)
library(modelr)
library(glmnet)  # For Lasso regression
library(knitr)  # For displaying MSE and coefficient tables

# Load and prepare the dataset
dataset1 <- read_excel("C:/Users/jomkr/OneDrive/Documents/Imperial_Y3/UROP/Data/CBD-IL-12 myeloid.xlsx", 
                       sheet = "Myeloid", range = "A15:Y81")
names(dataset1)[1] <- "mice"

# Process the data: Separate the 'mice' column and assign days
puredata <- dataset1 %>%
  separate(col = mice, into = c('analysis', 'number'), sep = 2) %>%
  mutate(number = as.integer(number)) %>%
  mutate(days = case_when(
    number >= 1 & number <= 4 ~ 0,
    number >= 5 & number <= 9 ~ 1,
    number >= 10 & number <= 14 ~ 2,
    number >= 15 & number <= 19 ~ 4,
    number >= 20 & number <= 24 ~ 7,
    number >= 25 & number <= 28 ~ 10,
    number >= 29 & number <= 34 ~ 13,
    number >= 35 ~ 14
  ))

# Reshape the data to long format
puredata_long <- puredata %>%
  pivot_longer(cols = c(3:26), names_to = 'metric', values_to = 'value')

# Filter for TM analysis and rename the metrics
TM_data <- puredata_long %>%
  filter(analysis == "TM") %>%
  mutate(metric = factor(metric, levels = c('CD45+', 'CD11b+CD45+', 'M-MDSC', 'G-MDSC', 'Macrophages', 'M2 macrophages', 
                                            'M1 macrophages', 'B cells', 'Endothelial cells', 'Dendritic cells', 
                                            '% M-MDSC/CD45', '% G-MDSC/CD45', '%macrophage/CD45', '%M2 macrophage/CD45', 
                                            '%M1 macrophage/CD45', '% B cells/CD45', '% Endothelial cells/CD45', 
                                            '%Dendritic cells/CD45', '# macrophages/tumor', '#M2 macrophages/tumor mg', 
                                            '# M1 macrophage/tumor mg', '#M-MDSC/tumor mg', '# B cells/tumor', 
                                            'tumor mg'))) 

# Reshape to wide format
TM_data_wide <- TM_data %>%
  pivot_wider(names_from = metric, values_from = value)

# Clean and standardize the column names
cleaned_column_names <- names(TM_data_wide) %>%
  as.character() %>%
  str_replace_all("[^a-zA-Z0-9]", "") %>%
  str_replace_all(" ", "_") %>%
  str_replace_all("\\+", "_plus") %>%
  str_replace_all("-", "_")

# Apply cleaned column names to the data frame
names(TM_data_wide) <- cleaned_column_names

# Adjust the covariate names based on existing columns
covariate_names <- c('CD45_plus', 'CD11bCD45', 'MMDSC', 'GMDSC', 'Macrophages', 'M2macrophages', 
                     'M1macrophages', 'Bcells', 'Endothelialcells', 'Dendriticcells', 'tumormg')

# Check covariates are actually present in the data
existing_covariates <- intersect(covariate_names, names(TM_data_wide))

# Standardize the covariates that exist
standardized_data <- TM_data_wide %>%
  mutate(across(all_of(existing_covariates), ~ scale(.))) %>%
  drop_na()

# Define days for training and testing
train_test_days <- list(
  "0" = 1,
  "1" = 2,
  "2" = 4,
  "4" = 7,
  "7" = 10,
  "10" = 13
)

# Initialize lists to store results, MSE values, and coefficients for each model
results <- list()
mse_results <- list()
coef_results <- list()

# Utility function to pad shorter coefficient vectors with NAs
pad_with_na <- function(coefs, target_length) {
  c(coefs, rep(NA, target_length - length(coefs)))
}

# Perform one-step cross-validation for each model
for (train_day in names(train_test_days)) {
  test_day <- train_test_days[[train_day]]
  
  # Split the data into training and testing
  train_data <- standardized_data %>% filter(days == as.numeric(train_day))
  test_data <- standardized_data %>% filter(days == test_day)
  
  # Check if the test set is empty
  if (nrow(test_data) == 0) next
  
  # Model 1: Simple Linear Model (Tumor Weight ~ Time)
  fit_lm_simple <- lm(tumormg ~ days, data = train_data)
  coef_lm_simple <- coef(fit_lm_simple)  # Extract coefficients
  predictions_lm_simple <- predict(fit_lm_simple, newdata = test_data)
  mse_lm_simple <- mean((test_data$tumormg - predictions_lm_simple)^2)
  
  # Model 2: Linear model with covariates
  model_formula_covariates <- as.formula(paste("tumormg ~", paste(existing_covariates[-length(existing_covariates)], collapse = " + ")))
  fit_lm_covariates <- lm(model_formula_covariates, data = train_data)
  coef_lm_covariates <- coef(fit_lm_covariates)  # Extract coefficients
  predictions_lm_covariates <- predict(fit_lm_covariates, newdata = test_data)
  mse_lm_covariates <- mean((test_data$tumormg - predictions_lm_covariates)^2)
  
  # Model 3: Lasso regression with covariates
  X_train <- model.matrix(~ ., data = train_data[, existing_covariates])
  X_test <- model.matrix(~ ., data = test_data[, existing_covariates])
  fit_lasso <- cv.glmnet(X_train, train_data$tumormg, alpha = 1, nfolds = 3)  # Lasso with fewer folds
  coef_lasso <- as.matrix(coef(fit_lasso, s = "lambda.min"))  # Extract coefficients for the optimal lambda
  predictions_lasso <- as.numeric(predict(fit_lasso, newx = X_test, s = "lambda.min"))  # Coerce to numeric vector
  mse_lasso <- mean((test_data$tumormg - predictions_lasso)^2)
  
  # Model 4: Horseshoe regression model
  fit_horseshoe <- stan_glm(model_formula_covariates,
                            data = train_data, 
                            family = gaussian(),
                            prior = hs(global_scale = 1, slab_scale = 2, slab_df = 4),  # Regularized horseshoe prior
                            chains = 4, iter = 2000, warmup = 1000, seed = 123,
                            control = list(adapt_delta = 0.999, max_treedepth = 15))  # Adjusted adapt_delta
  coef_horseshoe <- coef(fit_horseshoe)  # Extract coefficients
  predictions_horseshoe <- posterior_predict(fit_horseshoe, newdata = test_data)
  mean_predictions_horseshoe <- apply(predictions_horseshoe, 2, mean)
  mse_horseshoe <- mean((test_data$tumormg - mean_predictions_horseshoe)^2)
  
  # Calculate the maximum length of coefficients among all models
  max_coef_length <- max(length(coef_lm_simple), length(coef_lm_covariates), length(coef_lasso), length(coef_horseshoe))
  
  # Pad each model's coefficients with NA to ensure they have the same length
  coef_lm_simple_padded <- pad_with_na(coef_lm_simple, max_coef_length)
  coef_lm_covariates_padded <- pad_with_na(coef_lm_covariates, max_coef_length)
  coef_lasso_padded <- pad_with_na(coef_lasso, max_coef_length)
  coef_horseshoe_padded <- pad_with_na(coef_horseshoe, max_coef_length)
  
  # Store MSE results for this fold
  mse_results[[paste0("Train_Day_", train_day, "_Test_Day_", test_day)]] <- data.frame(
    train_day = as.numeric(train_day),
    test_day = test_day,
    mse_lm_simple = mse_lm_simple,
    mse_lm_covariates = mse_lm_covariates,
    mse_lasso = mse_lasso,
    mse_horseshoe = mse_horseshoe
  )
  
  # Store coefficients for each model
  coef_results[[paste0("Train_Day_", train_day, "_Test_Day_", test_day)]] <- data.frame(
    Variable = rownames(as.matrix(coef_lasso)),
    Simple_Linear_Model = coef_lm_simple_padded,
    Linear_Model_with_Covariates = coef_lm_covariates_padded,
    Lasso_Regression = coef_lasso_padded,
    Horseshoe_Regression = coef_horseshoe_padded
  )
  
  # Store prediction results for this fold (for plotting later)
  results[[paste0("Train_Day_", train_day, "_Test_Day_", test_day)]] <- data.frame(
    train_day = as.numeric(train_day),
    test_day = test_day,
    actual_tumor_weight = test_data$tumormg,
    pred_lm_simple = predictions_lm_simple,
    pred_lm_covariates = predictions_lm_covariates,
    pred_lasso = predictions_lasso,
    pred_horseshoe = mean_predictions_horseshoe
  )
}

# Combine all the results into a single data frame
results_df <- bind_rows(results)

# Combine all the MSE results into a single data frame
mse_df <- bind_rows(mse_results)

# Combine all the coefficient results into a single data frame
coef_df <- bind_rows(coef_results)

# Display the MSE table for each fold
kable(mse_df, caption = "MSE for Each Model Across Folds")

# Display the coefficient table for each model
kable(coef_df, caption = "Coefficient Values for Each Model")

# Add prediction intervals (50%, 80%, and 95% for all models) to results_df
results_df <- results_df %>%
  mutate(
    # 50% credible interval for Simple Linear Model
    lower_50_lm_simple = pred_lm_simple - 0.674 * sqrt(mse_lm_simple),
    upper_50_lm_simple = pred_lm_simple + 0.674 * sqrt(mse_lm_simple),
    # 80% credible interval for Simple Linear Model
    lower_80_lm_simple = pred_lm_simple - 1.282 * sqrt(mse_lm_simple),
    upper_80_lm_simple = pred_lm_simple + 1.282 * sqrt(mse_lm_simple),
    # 95% credible interval for Simple Linear Model
    lower_95_lm_simple = pred_lm_simple - 1.96 * sqrt(mse_lm_simple),
    upper_95_lm_simple = pred_lm_simple + 1.96 * sqrt(mse_lm_simple),
    
    # 50% credible interval for Linear Model with Covariates
    lower_50_lm_covariates = pred_lm_covariates - 0.674 * sqrt(mse_lm_covariates),
    upper_50_lm_covariates = pred_lm_covariates + 0.674 * sqrt(mse_lm_covariates),
    # 80% credible interval for Linear Model with Covariates
    lower_80_lm_covariates = pred_lm_covariates - 1.282 * sqrt(mse_lm_covariates),
    upper_80_lm_covariates = pred_lm_covariates + 1.282 * sqrt(mse_lm_covariates),
    # 95% credible interval for Linear Model with Covariates
    lower_95_lm_covariates = pred_lm_covariates - 1.96 * sqrt(mse_lm_covariates),
    upper_95_lm_covariates = pred_lm_covariates + 1.96 * sqrt(mse_lm_covariates),
    
    # 50% credible interval for Lasso Regression
    lower_50_lasso = pred_lasso - 0.674 * sqrt(mse_lasso),
    upper_50_lasso = pred_lasso + 0.674 * sqrt(mse_lasso),
    # 80% credible interval for Lasso Regression
    lower_80_lasso = pred_lasso - 1.282 * sqrt(mse_lasso),
    upper_80_lasso = pred_lasso + 1.282 * sqrt(mse_lasso),
    # 95% credible interval for Lasso Regression
    lower_95_lasso = pred_lasso - 1.96 * sqrt(mse_lasso),
    upper_95_lasso = pred_lasso + 1.96 * sqrt(mse_lasso),
    
    # 50% credible interval for Horseshoe Regression
    lower_50_hs = pred_horseshoe - 0.674 * sqrt(mse_horseshoe),
    upper_50_hs = pred_horseshoe + 0.674 * sqrt(mse_horseshoe),
    # 80% credible interval for Horseshoe Regression
    lower_80_hs = pred_horseshoe - 1.282 * sqrt(mse_horseshoe),
    upper_80_hs = pred_horseshoe + 1.282 * sqrt(mse_horseshoe),
    # 95% credible interval for Horseshoe Regression
    lower_95_hs = pred_horseshoe - 1.96 * sqrt(mse_horseshoe),
    upper_95_hs = pred_horseshoe + 1.96 * sqrt(mse_horseshoe)
  )

# Plotting the results for each model

# Adjust the results_df to ensure all predictions are aligned
results_df <- results_df %>%
  filter(!is.na(actual_tumor_weight))  # Remove any rows with missing actual tumor weights

# Add confidence intervals for each model
results_df <- results_df %>%
  mutate(
    # 50%, 80%, and 95% credible intervals for Simple Linear Model
    lower_50_lm_simple = pred_lm_simple - 0.674 * sqrt(mse_lm_simple),
    upper_50_lm_simple = pred_lm_simple + 0.674 * sqrt(mse_lm_simple),
    lower_80_lm_simple = pred_lm_simple - 1.282 * sqrt(mse_lm_simple),
    upper_80_lm_simple = pred_lm_simple + 1.282 * sqrt(mse_lm_simple),
    lower_95_lm_simple = pred_lm_simple - 1.96 * sqrt(mse_lm_simple),
    upper_95_lm_simple = pred_lm_simple + 1.96 * sqrt(mse_lm_simple),
    
    # 50%, 80%, and 95% credible intervals for Linear Model with Covariates
    lower_50_lm_covariates = pred_lm_covariates - 0.674 * sqrt(mse_lm_covariates),
    upper_50_lm_covariates = pred_lm_covariates + 0.674 * sqrt(mse_lm_covariates),
    lower_80_lm_covariates = pred_lm_covariates - 1.282 * sqrt(mse_lm_covariates),
    upper_80_lm_covariates = pred_lm_covariates + 1.282 * sqrt(mse_lm_covariates),
    lower_95_lm_covariates = pred_lm_covariates - 1.96 * sqrt(mse_lm_covariates),
    upper_95_lm_covariates = pred_lm_covariates + 1.96 * sqrt(mse_lm_covariates),
    
    # 50%, 80%, and 95% credible intervals for Lasso Regression
    lower_50_lasso = pred_lasso - 0.674 * sqrt(mse_lasso),
    upper_50_lasso = pred_lasso + 0.674 * sqrt(mse_lasso),
    lower_80_lasso = pred_lasso - 1.282 * sqrt(mse_lasso),
    upper_80_lasso = pred_lasso + 1.282 * sqrt(mse_lasso),
    lower_95_lasso = pred_lasso - 1.96 * sqrt(mse_lasso),
    upper_95_lasso = pred_lasso + 1.96 * sqrt(mse_lasso),
    
    # 50%, 80%, and 95% credible intervals for Horseshoe Regression
    lower_50_hs = pred_horseshoe - 0.674 * sqrt(mse_horseshoe),
    upper_50_hs = pred_horseshoe + 0.674 * sqrt(mse_horseshoe),
    lower_80_hs = pred_horseshoe - 1.282 * sqrt(mse_horseshoe),
    upper_80_hs = pred_horseshoe + 1.282 * sqrt(mse_horseshoe),
    lower_95_hs = pred_horseshoe - 1.96 * sqrt(mse_horseshoe),
    upper_95_hs = pred_horseshoe + 1.96 * sqrt(mse_horseshoe)
  )

# Plot 1: Simple Linear Model
ggplot(results_df, aes(x = test_day)) +
  geom_rect(aes(xmin = test_day - 0.1, xmax = test_day + 0.1, ymin = lower_95_lm_simple, ymax = upper_95_lm_simple, fill = "95% Credible Interval"), alpha = 0.2) +
  geom_rect(aes(xmin = test_day - 0.1, xmax = test_day + 0.1, ymin = lower_80_lm_simple, ymax = upper_80_lm_simple, fill = "80% Credible Interval"), alpha = 0.15) +
  geom_rect(aes(xmin = test_day - 0.1, xmax = test_day + 0.1, ymin = lower_50_lm_simple, ymax = upper_50_lm_simple, fill = "50% Credible Interval"), alpha = 0.2) +
  geom_point(aes(y = pred_lm_simple, color = "Predicted Tumor Weight"), size = 4, shape = 17) +
  geom_point(aes(y = actual_tumor_weight, color = "Actual Tumor Weight"), size = 3, shape = 16) +
  scale_fill_manual(name = "Credible Intervals", values = c("95% Credible Interval" = "#D1E5F0", "80% Credible Interval" = "#7ab5fa", "50% Credible Interval" = "#2786f5")) +
  scale_color_manual(name = "", values = c("Actual Tumor Weight" = "black", "Predicted Tumor Weight" = "red")) +
  labs(title = "Simple Linear Model: Predicted Tumor Weights with 50%, 80%, and 95% Credible Intervals", 
       x = "Test Day", y = "Tumor Weight") +
  ylim(-5, 15) +
  theme_minimal()

# Plot 2: Linear Model with Covariates
ggplot(results_df, aes(x = test_day)) +
  geom_rect(aes(xmin = test_day - 0.1, xmax = test_day + 0.1, ymin = lower_95_lm_covariates, ymax = upper_95_lm_covariates, fill = "95% Credible Interval"), alpha = 0.2) +
  geom_rect(aes(xmin = test_day - 0.1, xmax = test_day + 0.1, ymin = lower_80_lm_covariates, ymax = upper_80_lm_covariates, fill = "80% Credible Interval"), alpha = 0.15) +
  geom_rect(aes(xmin = test_day - 0.1, xmax = test_day + 0.1, ymin = lower_50_lm_covariates, ymax = upper_50_lm_covariates, fill = "50% Credible Interval"), alpha = 0.2) +
  geom_point(aes(y = pred_lm_covariates, color = "Predicted Tumor Weight (LM with Covariates)"), size = 4, shape = 17) +
  geom_point(aes(y = actual_tumor_weight, color = "Actual Tumor Weight"), size = 3, shape = 16) +
  scale_fill_manual(name = "Credible Intervals", values = c("95% Credible Interval" = "#D1E5F0", "80% Credible Interval" = "#7ab5fa", "50% Credible Interval" = "#2786f5")) +
  scale_color_manual(name = "", values = c("Actual Tumor Weight" = "black", "Predicted Tumor Weight (LM with Covariates)" = "red")) +
  labs(title = "Linear Model with Covariates: Predicted Tumor Weights with 50%, 80%, and 95% Credible Intervals", 
       x = "Test Day", y = "Tumor Weight") +
  ylim(-5, 12) +
  theme_minimal()

# Plot 3: Lasso Regression
ggplot(results_df, aes(x = test_day)) +
  geom_rect(aes(xmin = test_day - 0.1, xmax = test_day + 0.1, ymin = lower_95_lasso, ymax = upper_95_lasso, fill = "95% Credible Interval"), alpha = 0.2) +
  geom_rect(aes(xmin = test_day - 0.1, xmax = test_day + 0.1, ymin = lower_80_lasso, ymax = upper_80_lasso, fill = "80% Credible Interval"), alpha = 0.15) +
  geom_rect(aes(xmin = test_day - 0.1, xmax = test_day + 0.1, ymin = lower_50_lasso, ymax = upper_50_lasso, fill = "50% Credible Interval"), alpha = 0.2) +
  geom_point(aes(y = pred_lasso, color = "Predicted Tumor Weight (Lasso)"), size = 4, shape = 17) +
  geom_point(aes(y = actual_tumor_weight, color = "Actual Tumor Weight"), size = 3, shape = 16) +
  scale_fill_manual(name = "Credible Intervals", values = c("95% Credible Interval" = "#D1E5F0", "80% Credible Interval" = "#7ab5fa", "50% Credible Interval" = "#2786f5")) +
  scale_color_manual(name = "", values = c("Actual Tumor Weight" = "black", "Predicted Tumor Weight (Lasso)" = "red")) +
  labs(title = "Lasso Regression: Predicted Tumor Weights with 50%, 80%, and 95% Credible Intervals", 
       x = "Test Day", y = "Tumor Weight") +
  ylim(-5, 15) +
  theme_minimal()

# Plot 4: Horseshoe Regression
ggplot(results_df, aes(x = test_day)) +
  geom_rect(aes(xmin = test_day - 0.1, xmax = test_day + 0.1, ymin = lower_95_hs, ymax = upper_95_hs, fill = "95% Credible Interval"), alpha = 0.2) +
  geom_rect(aes(xmin = test_day - 0.1, xmax = test_day + 0.1, ymin = lower_80_hs, ymax = upper_80_hs, fill = "80% Credible Interval"), alpha = 0.15) +
  geom_rect(aes(xmin = test_day - 0.1, xmax = test_day + 0.1, ymin = lower_50_hs, ymax = upper_50_hs, fill = "50% Credible Interval"), alpha = 0.2) +
  geom_point(aes(y = pred_horseshoe, color = "Predicted Tumor Weight (Horseshoe)"), size = 4, shape = 17) +
  geom_point(aes(y = actual_tumor_weight, color = "Actual Tumor Weight"), size = 3, shape = 16) +
  scale_fill_manual(name = "Credible Intervals", values = c("95% Credible Interval" = "#D1E5F0", "80% Credible Interval" = "#7ab5fa", "50% Credible Interval" = "#2786f5")) +
  scale_color_manual(name = "", values = c("Actual Tumor Weight" = "black", "Predicted Tumor Weight (Horseshoe)" = "red")) +
  labs(title = "Horseshoe Regression Model: Predicted Tumor Weights with 50%, 80%, and 95% Credible Intervals", 
       x = "Test Day", y = "Tumor Weight") +
  ylim(-5, 15) +
  theme_minimal()
