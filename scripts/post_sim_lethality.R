library(mgcv)


# Simulate posterior for 

# Step 1: Assuming `newdata` is the data frame of new observations; if you want to use the original data, replace `newdata` with your original dataset
newdata <- data.frame(yse = seq(min(reached5$yse), max(reached5$yse), length.out = 100))
lpmatrix <- predict(m1, newdata, type = "lpmatrix")

# Step 2: Obtain the covariance matrix of the coefficients
cov_mat <- vcov(m1)

# Step 3: Simulate new coefficients and predictions
n_sim <- 100000 # Number of simulations
simulated_coefs <- MASS::mvrnorm(n = n_sim, mu = coef(m1), Sigma = cov_mat)
simulated_preds <- lpmatrix %*% t(simulated_coefs) # Each column is one set of simulated predictions

# Step 4: Apply the inverse link function for binomial family (logit link)
simulated_probs <- plogis(simulated_preds)

# Calculate summary statistics (e.g., mean, quantiles) from simulations
mean_probs <- apply(simulated_probs, 1, median)
lower_ci <- apply(simulated_probs, 1, quantile, probs = 0.025)
upper_ci <- apply(simulated_probs, 1, quantile, probs = 0.975)

# Combine with `newdata` for plotting or further analysis
results_df <- data.frame(yse = newdata$yse, mean_probs, lower_ci, upper_ci)

# Plotting (assuming `ggplot2` is loaded)
library(ggplot2)
ggplot(results_df, aes(x = yse)) +
  geom_line(aes(y = mean_probs), color = "blue") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2) +
  lims(x = c(0,20), y = c(0,1)) +
  theme(aspect.ratio = 1)


library(mgcv)
library(MASS) # For mvrnorm

# Assuming `m2` is your model fitted with negative binomial family
# Step 1: Generate the lpmatrix for new or existing data
# Here, I assume you want to generate predictions for a range of `yse` values
newdata <- data.frame(yse = seq(min(dt5$yse), max(dt5$yse), length.out = 100))
lpmatrix <- predict(m2, newdata, type = "lpmatrix")

# Step 2: Obtain the covariance matrix of the coefficients
cov_mat <- vcov(m2)

# Step 3: Simulate new coefficients and predictions
n_sim <- 100000 # Number of simulations
simulated_coefs <- mvrnorm(n = n_sim, mu = coef(m2), Sigma = cov_mat)
simulated_preds <- lpmatrix %*% t(simulated_coefs) # Each column is one set of simulated predictions

# Step 4: Apply the inverse link function for the nb family (log link here)
# For the nb family, the expected counts are exp(predicted linear predictors)
simulated_counts <- exp(simulated_preds)

# Calculate summary statistics (e.g., mean, quantiles) from simulations
mean_counts <- apply(1/simulated_counts, 1, median)
lower_ci <- apply(1/simulated_counts, 1, quantile, probs = 0.025)
upper_ci <- apply(1/simulated_counts, 1, quantile, probs = 0.975)

# Combine with `newdata` for plotting or further analysis
results_df <- data.frame(yse = newdata$yse, mean_counts, lower_ci, upper_ci)

# Plotting
library(ggplot2)
ggplot(results_df, aes(x = yse)) +
  geom_line(aes(y = mean_counts), color = "blue") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2, fill = "blue")




