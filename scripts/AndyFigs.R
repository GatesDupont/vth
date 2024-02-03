library(tidyverse)
library(cowplot)
library(vroom)
library(mgcv)
library(scam)

df <- vroom("data/andre_hofiData_gcleaned_final_20240131.csv") %>%
  mutate(log_mg = log(mg_abundance)) %>%
  mutate(log_mg = if_else(mg_abundance == 0, 0, log_mg))


# ---- Figure 1 ----

p1 <- ggplot(df, aes(x = log_mg, y = eyescore_sum, color = strain)) +
  geom_point(alpha = 0.1, size = 4) +
  geom_smooth(method = "gam", method.args = list(family = nb(), gamma = 1.4), se = F) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  scale_color_manual(values = c("east" = "blue", "west" = "red")) +
  labs(x = "MG abundance (log)", y = "Eyescore sum", color = NULL) +
  theme_minimal(14) +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 0.2, fill = NA),
        aspect.ratio = 1,
        text = element_text(family = "lato"))

ggplot(df, aes(x = log_mg, y = eyescore_sum, color = strain)) +
  # geom_point(alpha = 0.1, size = 4) +
  geom_smooth(method = "gam", method.args = list(family = nb(), gamma = 1.4), se = F) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  scale_color_manual(values = c("east" = "blue", "west" = "red")) +
  labs(x = "MG abundance (log)", y = "Eyescore sum", color = NULL) +
  theme_minimal(14) +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 0.2, fill = NA),
        aspect.ratio = 1,
        text = element_text(family = "lato"))


# ---- Figure 2 ----

tmp_pred <- list()
for(i in 1:6){
  
  tmp_df <- df %>%
    filter(eyescore_sum == i)
  
  # tmp_m <- lm(log_mg ~ 0 + DPI, data = tmp_df)
  tmp_m <- lm(log_mg ~ DPI, data = tmp_df)
  # tmp_m <- scam(log_mg ~ s(DPI, bs = "mpi", k = 4), data = tmp_df)
  
  tmp_pred[[i]] <- expand.grid(
    eyescore_sum = i, 
    DPI = seq(0, 71, by = 0.2))
  
  tmp_pred[[i]]$pred <- predict(
    object = tmp_m, 
    newdata = tmp_pred[[i]], 
    type = "response")
  
}
eyescore_pred_df <- do.call(rbind, tmp_pred)

tmp_pred <- list()
for(i in 1:2){
  
  strain <- c("east","west")
  
  if(i == 1){
    tmp_df <- df %>%
      filter(strain == "east")
  }
  
  if(i == 2){
    tmp_df <- df %>%
      filter(strain == "west")
  }

  
  tmp_m <- lm(log_mg ~ 0 + DPI, data = tmp_df)
  
  tmp_pred[[i]] <- expand.grid(
    strain = strain[i], 
    DPI = seq(0, 71, by = 0.2))
  
  tmp_pred[[i]]$pred <- predict(
    object = tmp_m, 
    newdata = tmp_pred[[i]], 
    type = "response")
  
}
strain_pred_df <- do.call(rbind, tmp_pred)

p2 <- ggplot(df, aes(x = DPI, y = log_mg)) +
  geom_point(alpha = 0.1, size = 4, na.rm = T) +
  geom_line(data = eyescore_pred_df, 
            aes(x = DPI, color = as.factor(eyescore_sum), group = eyescore_sum, y = pred), 
            linewidth = 1, linetype = "dashed") +
  geom_line(data = strain_pred_df %>% filter(strain == "east"), 
            aes(x = DPI, y = pred), 
            color = "blue",
            linewidth = 1) +
  geom_line(data = strain_pred_df %>% filter(strain == "west"), 
            aes(x = DPI, y = pred), 
            color = "red",
            linewidth = 1) +
  scale_color_viridis_d() +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  labs(x = "Days post inoculation", y = "MG abundance (log)", color = "Eyescore") +
  theme_minimal(14) +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 0.2, fill = NA),
        aspect.ratio = 1,
        text = element_text(family = "lato"));p2

plot_grid(p1, p2, nrow = 2, align = "hv")

