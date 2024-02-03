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


m <- gam(reached5 ~ yse, data = reached5, family = "binomial")

pred_df <- data.frame(yse = seq(0,20,0.1))

pred_m <- predict(m, pred_df, type = "link", se = T)

pred_df$mean <- plogis(pred_m$fit)
pred_df$lwr <- plogis(pred_m$fit - 1.96 * pred_m$se.fit)
pred_df$upr <- plogis(pred_m$fit + 1.96 * pred_m$se.fit)

p1 <- ggplot(pred_df, aes(x = yse, y = mean)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4, fill = 8) +
  geom_line(size = 1) +
  geom_point(data = reached5, aes(x = yse, y = reached5), 
             size = 4, pch = 21, stroke = 2, alpha = 0.4, color = 2) +
  labs(x = "Years since disease emergence", caption = "p-val = 0.008 **",
       title = "Which birds reached lethal eyescore?",
       y = "Probability of reaching lethal eyescore") +
  theme_minimal(14) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        panel.border = element_rect(linewidth = 0.2, fill = NA),
        aspect.ratio = 1)



# ---- Modeling days to reach eyescore 5 ----

# If it reached 5, how many days?
dt5 <- df %>%
  group_by(bird_id) %>%
  mutate(lethal_eyescore = if_else(eyescore_sum >= 5, 1, 0)) %>%
  filter(lethal_eyescore == 1) %>%
  summarise(dt5 = min(DPI)) %>%
  ungroup() %>%
  left_join(df %>% select(bird_id, yse) %>% distinct())

p2 <- ggplot(dt5, aes(x = yse, y = dt5)) +
  geom_point(size = 4, pch = 21, stroke = 2, alpha = 0.4, color = 2) +
  geom_smooth(method = "gam", formula = y ~ x, method.args = list(family = nb()), size = 1, color = 1) +
  labs(x = "Years since disease emergence", y = "Days to lethal eyescore",
       title  = "For birds that did, how long did it take?") +
  theme_minimal(14) +
  theme(aspect.ratio = 1,
        # text = element_text(family = "Lato"),
        panel.border = element_rect(fill = NA, linewidth = 0.2),
        panel.grid.minor = element_blank()); p2

plot_grid(p1, p2, nrow = 1, align = "hv")

