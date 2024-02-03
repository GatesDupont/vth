library(tidyverse)
library(vroom)
library(rcartocolor)
library(mgcv)
library(cowplot)
library(extrafont)
library(MASS)
select <- dplyr::select

# Number of simulations
n_sim <- 100000 


# ---- Organize data ----

df <- vroom("data/andre_hofiData_gcleaned_final_20240131.csv") %>%
  filter(status != "dead") %>%
  select(-1, -rpa_antibodies) %>% 
  select(-set, -aviary, -eyescore_left, -eyescore_right) %>%
  mutate(log_mg = log(mg_abundance)) %>%
  mutate(log_mg = if_else(mg_abundance == 0, 0, log_mg)) %>%
  mutate(year = str_sub(isolate, 3, -1L)) %>%
  select(-isolate, -mg_abundance) %>%
  mutate(year = as.numeric(as.character(year))) %>%
  mutate(first_year = if_else(strain == "east", 1993, 2005)) %>%
  mutate(yse = year - first_year) %>%
  select(-first_year, -year) %>% 
  na.omit() %>%
  mutate(bird_id = as.factor(bird_id)) %>%
  mutate(strain = as.factor(strain)) %>%
  mutate(sex = as.factor(sex))


# ---- BINOMIAL MODEL ----

# Did it reach eye score 5?
reached5 <- df %>%
  group_by(bird_id) %>%
  mutate(maxeye = max(eyescore_sum)) %>%
  ungroup() %>%
  select(strain, bird_id, sex, yse, maxeye) %>%
  distinct() %>%
  mutate(reached5 = if_else(maxeye >= 5, 1, 0))

# Model
m1 <- gam(reached5 ~ yse, data = reached5, family = "binomial", gamma = 1.4)

# Step 1: Assuming `newdata` is the data frame of new observations; 
# if you want to use the original data, replace `newdata` with your original dataset
newdata_m1 <- data.frame(yse = seq(min(df$yse), max(df$yse), length.out = 100))
lpmatrix_m1 <- predict(m1, newdata_m1, type = "lpmatrix")

# Step 2: Obtain the covariance matrix of the coefficients
cov_mat_m1 <- vcov(m1)

# Step 3: Simulate new coefficients and predictions
simulated_coefs_m1 <- MASS::mvrnorm(n = n_sim, mu = coef(m1), Sigma = cov_mat_m1)
simulated_preds_m1 <- lpmatrix_m1 %*% t(simulated_coefs_m1) # Each column is one set of simulated predictions

# Step 4: Apply the inverse link function for binomial family (logit link)
simulated_probs_m1 <- plogis(simulated_preds_m1)

# Calculate summary statistics (e.g., mean, quantiles) from simulations
mean_probs_m1 <- apply(simulated_probs_m1, 1, median)
lower_ci_m1 <- apply(simulated_probs_m1, 1, quantile, probs = 0.025)
upper_ci_m1 <- apply(simulated_probs_m1, 1, quantile, probs = 0.975)

# Combine with `newdata` for plotting or further analysis
results_df_m1 <- data.frame(yse = newdata_m1$yse, mean_probs_m1, lower_ci_m1, upper_ci_m1)

# Plotting (assuming `ggplot2` is loaded)
ggplot(results_df_m1, aes(x = yse)) +
  geom_line(aes(y = mean_probs_m1), color = "blue") +
  geom_ribbon(aes(ymin = lower_ci_m1, ymax = upper_ci_m1), alpha = 0.2) +
  lims(x = c(0,20), y = c(0,1)) +
  theme(aspect.ratio = 1)


# ---- NEGBIN MODEL ----

# If it reached 5, how many days?
dt5 <- df %>%
  group_by(bird_id) %>%
  mutate(lethal_eyescore = if_else(eyescore_sum >= 5, 1, 0)) %>%
  filter(lethal_eyescore == 1) %>%
  summarise(dt5 = min(DPI)) %>%
  ungroup() %>%
  left_join(df %>% select(bird_id, yse) %>% distinct())

# Model
m2 <- gam(dt5 ~ yse, data = dt5, family = "nb")

# Assuming `m2` is your model fitted with negative binomial family
# Step 1: Generate the lpmatrix for new or existing data
# Here, I assume you want to generate predictions for a range of `yse` values
newdata_m2 <- data.frame(yse = seq(min(df$yse), max(df$yse), length.out = 100))
lpmatrix_m2 <- predict(m2, newdata_m2, type = "lpmatrix")

# Step 2: Obtain the covariance matrix of the coefficients
cov_mat_m2 <- vcov(m2)

# Step 3: Simulate new coefficients and predictions
simulated_coefs_m2 <- mvrnorm(n = n_sim, mu = coef(m2), Sigma = cov_mat_m2)
simulated_preds_m2 <- lpmatrix_m2 %*% t(simulated_coefs_m2) # Each column is one set of simulated predictions

# Step 4: Apply the inverse link function for the nb family (log link here)
# For the nb family, the expected counts are exp(predicted linear predictors)
# take inverse for rate
simulated_counts_m2 <- 1/exp(simulated_preds_m2)

# Calculate summary statistics (e.g., mean, quantiles) from simulations
mean_counts_m2 <- apply(simulated_counts_m2, 1, median)
lower_ci_m2 <- apply(simulated_counts_m2, 1, quantile, probs = 0.025)
upper_ci_m2 <- apply(simulated_counts_m2, 1, quantile, probs = 0.975)

# Combine with `newdata` for plotting or further analysis
results_df_m2 <- data.frame(yse = newdata_m2$yse, mean_counts_m2, lower_ci_m2, upper_ci_m2)

# Plotting
ggplot(results_df_m2, aes(x = yse)) +
  geom_line(aes(y = mean_counts_m2), color = "blue") +
  geom_ribbon(aes(ymin = lower_ci_m2, ymax = upper_ci_m2), alpha = 0.2, fill = "blue")


# ---- COMBINE -----

# Simulated viurlence from posterior draws
simulated_virulence <- simulated_probs_m1 * simulated_counts_m2

# Calculate summary statistics (e.g., mean, quantiles) from simulations
mean_counts_vir <- apply(simulated_virulence, 1, median)
lower_ci_vir <- apply(simulated_virulence, 1, quantile, probs = 0.025)
upper_ci_vir <- apply(simulated_virulence, 1, quantile, probs = 0.975)

# Combine with `newdata` for plotting or further analysis
results_df_vir <- data.frame(yse = newdata_m1$yse, mean_counts_vir, lower_ci_vir, upper_ci_vir)

# Combined plot
ggplot(results_df_m2, aes(x = yse)) +
  geom_line(aes(y = mean_counts_vir), color = "blue", size = 1.3) +
  geom_ribbon(aes(ymin = lower_ci_vir, ymax = upper_ci_vir), alpha = 0.2, fill = "blue") +
  labs(x = "Years since disease emergence",
       title = "Virulence (?)",
       y = "Pr(reached lethal eyescore) × (1/Days_to_lethal_eyescore)") +
  theme_minimal(16) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        panel.border = element_rect(linewidth = 0.2, fill = NA),
        aspect.ratio = 1,
        text = element_text(family = "Lato"))
