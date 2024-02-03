library(tidyverse)
library(vroom)
library(rcartocolor)
library(mgcv)
library(cowplot)
library(extrafont)


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


# ---- Modeling whether or not they reached eyescore 5 ----

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

# Prediction data frame
pred_df1 <- data.frame(yse = seq(0,20,0.1))

# Predictions on link scale
pred_m1 <- predict(m1, pred_df1, type = "link", se = T)

pred_df1$mean <- plogis(pred_m1$fit)
pred_df1$lwr <- plogis(pred_m1$fit - 1.96 * pred_m1$se.fit)
pred_df1$upr <- plogis(pred_m1$fit + 1.96 * pred_m1$se.fit)

ggplot(pred_df1, aes(x = yse, y = mean)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4, fill = 8) +
  geom_line(size = 1)


# ---- Modeling length of infection ----

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

# Prediction data frame
pred_df2 <- data.frame(yse = seq(0,20,0.1))

# Predictions on link scale
pred_m2 <- predict(m2, pred_df2, type = "link", se = T)

# ---- Comibne ----

# Parameters

# - intercept
m1_int <- coef(m1)[1] # Coefficient estimate for yse
m1_int_se <- summary(m1)$se[1] # Standard error for yse coefficient

# - slope
m1_beta <- coef(m1)[2] # Coefficient estimate for yse
m1_beta_se <- summary(m1)$se[2] # Standard error for yse coefficient

# Setup loop
n_samples <- 10000 # Number of Monte Carlo samples
cov_yse <- seq(0,20,by=1)
results_list <- list()

# Loop
for(i in 1:n_samples){
  
  # Generate Monte Carlo samples from a normal distribution
  sample_int <- rnorm(n = 1, mean = m1_int, sd = m1_int_se)
  sample_beta <- rnorm(n = 1, mean = m1_beta, sd = m1_beta_se)
  
  # Predict line
  tmp_line <- plogis(sample_int + sample_beta * cov_yse)
  
  # Append results
  results_list[[i]] <- data.frame(
    predicted_prob = tmp_line,
    yse = cov_yse,
    iteration = i
  )
  
}

results <- do.call(rbind, results_list) %>%
  group_by(yse) %>%
  summarise(lwr = quantile(predicted_prob, 0.025),
            upr = quantile(predicted_prob, 0.975),
            mean = median(predicted_prob)) %>%
  ungroup()

ggplot() +
  geom_ribbon(data = results, aes(x = yse, y = mean, ymin = lwr, ymax = upr), fill = 2, alpha = 0.3) +
  geom_line(data = results, aes(x = yse, y = mean), color = 2) +
  geom_ribbon(data = pred_df1, aes(x = yse, y = mean, ymin = lwr, ymax = upr), color = 8, fill = NA, size = 1) +
  geom_line(data = pred_df1, aes(x = yse, y = mean), size = 1)




