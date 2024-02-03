library(tidyverse)
library(vroom)
library(rcartocolor)
library(cowplot)
library(scam)


# ---- Load data ----

df <- vroom("data/andre_hofiData_gcleaned_final_20240131.csv") %>%
  mutate(log_mg = log(mg_abundance)) %>%
  mutate(log_mg = if_else(mg_abundance == 0, 0, log_mg)) %>%
  mutate(year = str_sub(isolate, 3, -1L)) %>%
  mutate(year = as.numeric(as.character(year))) %>%
  mutate(first_year = if_else(strain == "east", 1993, 2005)) %>%
  mutate(years_since_emergence = year - first_year) %>%
  mutate(bird_id = as.factor(bird_id)) %>%
  mutate(strain = as.factor(strain))


# # ---- Fit moodel ----
# 
# m <- gam(log_mg ~ 0 + s(DPI) + s(years_since_emergence, k = 8) + s(strain, bs = "re") + s(bird_id, bs = "re"),
#          family = nb(),
#          data = df,
#          gamma = 1.4)
# 
# library(scam)
# 
# m <- scam(log_mg ~ 0 + s(DPI, bs = "mpi", k = 6) + s(years_since_emergence, k = 8) + s(strain, bs = "re") + s(bird_id, bs = "re"),
#           data = df,
#           gamma = 1.4)
# 
# 
# # ---- Predict ----
# 
# pred_df <- expand.grid(
#   DPI = seq(0, 71, by = 0.5),
#   years_since_emergence = seq(0, 20),
#   strain = sample(x = df$strain, size = 1),
#   bird_id = sample(x = df$bird_id, size = 1))
# 
# pred_df$log_mg <- predict(m, pred_df, type = "response", exclude = c("s(strain)", "s(bird_id)"))
# 
# 
# # ---- Figure ----
# 
# ggplot(pred_df, aes(x = DPI, y = log_mg, color = years_since_emergence, group = years_since_emergence)) +
#   # geom_rug(data = df, aes(y = log_mg, x = DPI), alpha = 0.5) +
#   geom_line(linewidth = 2) +
#   scale_color_gradientn(colors = c(2,6,4,3)) +
#   scale_y_continuous(breaks = scales::pretty_breaks()) +
  # labs(x = "Days post inoculation",
  #      y = "MG abundance (log)") +
  # guides(color = guide_legend(title.position = "top",
  #                             title = "Years since emergence",
  #                             title.hjust = 0.5)) +
  # theme_minimal(14) +
  # theme(panel.grid.minor = element_blank(),
  #       legend.position = "bottom",
  #       panel.border = element_rect(linewidth = 0.2, fill = NA),
  #       aspect.ratio = 1,
  #       text = element_text(family = "lato"))







df_sum <- df %>%
  filter(!is.na(mg_abundance)) %>%
  group_by(bird_id) %>%
  summarise(n = sum(mg_abundance > 0)) %>%
  ungroup() %>%
  filter(n > 2)

birds <- df_sum%>%
  pull(bird_id) %>%
  as.character()
  
n_birds <- length(birds)


bird_beta <- c()
bird_fit_i <- list()
for(i in 1:n_birds){
  
  tmp_df <- df %>%
    filter(bird_id == birds[i])
  
  # tmp_m <- scam(
  #   log_mg ~ 0 + DPI,
  #   data = tmp_df,
  #   gamma = 1.4)
  
  tmp_m <- gam(log_mg ~ 0 + DPI, data = tmp_df)
  
  tmp_pred_df <- data.frame(
    DPI = seq(0, 71, by = 0.01))
  
  bird_beta[i] <- coef(tmp_m)[1]
  
  tmp_pred_df$log_mg_pred <- predict(tmp_m, tmp_pred_df, type = "response")
  
  bird_fit_i[[i]] <- tmp_pred_df %>%
    mutate(bird_id = birds[i]) %>%
    mutate(years_since_emergence = unique(tmp_df$years_since_emergence))
  
}

bird_fits <- do.call(rbind, bird_fit_i)


ggplot(bird_fits, aes(x = DPI, y = log_mg_pred, color = years_since_emergence, group = bird_id)) +
  geom_line(linewidth = 1) +
  scale_color_gradientn(colors = c(2,6,4,3)) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  labs(x = "Days post inoculation", 
       y = "MG abundance (log)") + 
  guides(color = guide_legend(title.position = "top",
                              title = "Years since emergence",
                              title.hjust = 0.5)) +
  theme_minimal(14) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        panel.border = element_rect(linewidth = 0.2, fill = NA),
        aspect.ratio = 1,
        text = element_text(family = "lato"))



betas_df <- data.frame(
  bird = birds,
  beta = bird_beta) %>%
  left_join(bird_fits %>% 
              select(bird = bird_id, years_since_emergence) %>% 
              distinct(),
            by = "bird") %>%
  left_join(df_sum, by = c("bird" = "bird_id"))

beta_m <- glm(beta ~ years_since_emergence, data = betas_df, family = gaussian(link = "log"), weights = betas_df$n)

beta_pred_df <- data.frame(
  years_since_emergence = seq(0, max(betas_df$years_since_emergence), by = 0.1)
)

beta_preds <- predict(beta_m, beta_pred_df, type = "link", se = T)

beta_pred_df$mean <- exp(beta_preds$fit)
beta_pred_df$lwr <- exp(beta_preds$fit - 1.96 * beta_preds$se.fit)
beta_pred_df$upr <- exp(beta_preds$fit + 1.96 * beta_preds$se.fit)

ggplot(betas_df, aes(x = years_since_emergence, y = beta)) +
  geom_point(alpha = 0.4, aes(size = n)) +
  geom_ribbon(data = beta_pred_df, alpha = 0.2,
              aes(x = years_since_emergence, y = mean, ymin = lwr, ymax = upr)) +
  geom_line(data = beta_pred_df, aes(x = years_since_emergence, y = mean)) +
  labs(x = "Years since emergence", y = "Rate of increae of MG") +
  guides(size = guide_legend(title.position = "top",
                              title = "Number of MG observations\ngreater than 0",
                              title.hjust = 0.5)) +
  theme_minimal(14) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        panel.border = element_rect(linewidth = 0.2, fill = NA),
        aspect.ratio = 1,
        text = element_text(family = "lato"))




# Here's a more basic way of looking at it.
# 
# I fit a simple linear model to each bird individually: log_mg ~ days post infection (intercept fixed at (0,0)).
# 
# Then I just kept the slope of that line for each bird.
# 
# Then I plotted that here (y-axis) against how old the strain was relative to its emergence (x-axis).
# 
# Point size relates the number of experimental observations for that bird that had an MG abundance greater than 0.
























