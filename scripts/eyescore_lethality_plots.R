library(tidyverse)
library(vroom)
library(rcartocolor)
library(mgcv)

df <- vroom("data/andre_hofiData_gcleaned_final_20240131.csv") %>%
  mutate(log_mg = log(mg_abundance)) %>%
  mutate(log_mg = if_else(mg_abundance == 0, 0, log_mg)) %>%
  mutate(year = str_sub(isolate, 3, -1L)) %>%
  mutate(year = as.numeric(as.character(year))) %>%
  mutate(first_year = if_else(strain == "east", 1993, 2005)) %>%
  mutate(years_since_emergence = year - first_year) %>%
  mutate(bird_id = as.factor(bird_id)) %>%
  mutate(strain = as.factor(strain)) %>%
  mutate(sex = as.factor(sex))

new_bird_id_order <- df %>% 
  group_by(bird_id) %>%
  mutate(maxeye = max(eyescore_sum)) %>%
  ungroup() %>%
  select(bird_id, years_since_emergence, maxeye) %>%
  distinct() %>%
  mutate(bird_id = as.character(bird_id)) %>%
  arrange(maxeye) %>%
  arrange(years_since_emergence) %>%
  pull(bird_id)

df2 <- df %>%
  mutate(bird_id = factor(bird_id, levels = new_bird_id_order))
  
ggplot(data = df2, aes(x = DPI, y = eyescore_sum, color = years_since_emergence)) +
  facet_wrap(~bird_id, nrow = 5) +
  geom_point(na.rm = T) +
  geom_line(aes(group = bird_id), na.rm = T, linewidth = 1) +
  scale_color_carto_c(palette = "Temps") +
  labs(x = "Days post inoculation",
       y = "Eyescore total") +
  guides(color = guide_legend(title.position = "top",
                              title = "Years since emergence",
                              title.hjust = 0.5)) +
  theme_minimal(14) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        panel.border = element_rect(linewidth = 0.2, fill = NA),
        aspect.ratio = 1,
        text = element_text(family = "lato"))



# m <- gam(eyescore_sum ~ 0 + DPI:years_since_emergence, data = df)

m <- gam(eyescore_sum ~ s(DPI, by = years_since_emergence, bs = "fs", m = 2) + s(bird_id, bs = "re"), 
         data = df, 
         gamma = 1.4, 
         family = gaussian(link = "log"))

pred_df <- expand.grid(
  bird_id = sample(unique(df$bird_id), 1),
  DPI = seq(0,71,0.1),
  years_since_emergence = seq(0,20)) %>%
  as_tibble()

pred_df$eyescore_sum_pred = predict(m, pred_df, type = "response", exclude = "s(bird_id)")


ggplot() +
  geom_line(data = pred_df, size = 2,
            aes(x = DPI, y = eyescore_sum_pred, color = years_since_emergence, group = years_since_emergence)) +
  geom_hline(yintercept = 5) +
  scale_color_carto_c(palette = "Temps") +
  labs(x = "Days post inoculation",
       y = "Eyescore total") +
  guides(color = guide_legend(title.position = "top",
                              title = "Years since emergence",
                              title.hjust = 0.5)) +
  theme_minimal(14) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        panel.border = element_rect(linewidth = 0.2, fill = NA),
        aspect.ratio = 1,
        text = element_text(family = "lato"))



ggplot(data = df, aes(x = log_mg, y = eyescore_sum)) +
  geom_point(na.rm = T, size = 4, alpha = 0.5, color = 8) +
  geom_smooth(method = "gam", method.args = list(gamma = 1.4)) +
  labs(x = "MG abundance (log)",
       y = "Eyescore total") +
  theme_minimal(14) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        panel.border = element_rect(linewidth = 0.2, fill = NA),
        aspect.ratio = 1,
        text = element_text(family = "lato"))


testm <- gam(eyescore_sum ~ s(log_mg) + s(bird_id, bs = "re"), data = df, family = nb(), gamma = 1.4)
# testm <- scam(eyescore_sum ~ s(log_mg, bs = "mpi") + s(bird_id, bs = "re"), data = df)
plot(testm, shade = T, shift = coef(testm)[1], trans = exp, select = 1, 
     xlab = "MG abundance (log)", ylab = "Predicted eyescore")
