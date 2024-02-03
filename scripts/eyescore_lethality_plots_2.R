library(tidyverse)
library(vroom)
library(rcartocolor)
library(mgcv)

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


ggplot(df, aes(x = DPI, y = eyescore_sum)) +
  facet_wrap(~bird_id) +
  geom_line(color = 2) +
  geom_point(color = 2) +
  geom_hline(yintercept = 5, linetype = "dashed")


# Did it reach eye score 5?
reached5 <- df %>%
  group_by(bird_id) %>%
  mutate(maxeye = max(eyescore_sum)) %>%
  ungroup() %>%
  select(strain, bird_id, sex, yse, maxeye) %>%
  distinct() %>%
  mutate(reached5 = if_else(maxeye >= 5, 1, 0))


m <- gam(reached5 ~ s(yse, k = 8), data = reached5, family = "binomial", gamma = 1.4)

plot(m, shift = coef(m)[1], trans = plogis, shade = T, ylim = c(0,1), rug = F,
     # residuals = T, cex = 2, rug = F, pch = 20,
     ylab = "Probability of reaching lethal eyescore", 
     xlab = "Years since emergence")

# If it reached 5, how many days?
dt5 <- df %>%
  group_by(bird_id) %>%
  mutate(maxeye = max(eyescore_sum)) %>%
  ungroup() %>%
  filter(maxeye >= 5) %>%
  group_by(bird_id) %>%
  summarise(dt5 = min(DPI)) %>%
  ungroup()





m <- gam(reached5 ~ yse, data = reached5, family = "binomial")

pred_df <- data.frame(yse = seq(0,20,0.1))

pred_m <- predict(m, pred_df, type = "link", se = T)

pred_df$mean <- plogis(pred_m$fit)
pred_df$lwr <- plogis(pred_m$fit - 1.96 * pred_m$se.fit)
pred_df$upr <- plogis(pred_m$fit + 1.96 * pred_m$se.fit)

ggplot(pred_df, aes(x = yse, y = mean)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4, fill = 8) +
  geom_line(size = 1) +
  geom_point(data = reached5, aes(x = yse, y = reached5), size = 3, pch = 21) +
  labs(x = "Years since disease emergence", caption = "p-val = 0.008 **",
       y = "Probability of reaching lethal eyescore") +
  theme_minimal(14) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        panel.border = element_rect(linewidth = 0.2, fill = NA),
        aspect.ratio = 1,
        text = element_text(family = "lato"))





