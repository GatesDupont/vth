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


ggplot(df, aes(x = DPI, y = log_mg, color = yse)) +
  facet_wrap(~bird_id) +
  geom_line() +
  geom_point()


days_to_max_mg <- df %>%
  group_by(bird_id) %>%
  mutate(max_log_mg = max(log_mg)) %>%
  filter(log_mg == max_log_mg) %>%
  summarise(days_to_max_mg = min(DPI),
            max_log_mg = unique(max_log_mg),
            yse = unique(yse)) %>%
  ungroup() %>%
  filter(max_log_mg != 0)


birds_by_daystomax <- days_to_max_mg %>%
  arrange((days_to_max_mg)) %>%
  pull(bird_id) %>%
  as.character()

df %>%
  mutate(bird_id = factor(bird_id, levels = birds_by_daystomax)) %>%
  na.omit() %>%
  ggplot(aes(x = DPI, y = log_mg, color = yse)) +
  facet_wrap(~bird_id, nrow = 5) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_carto_c(palette = "Temps") +
  labs(x = "Days post inoculation",
       y = "MG (log)") +
  guides(color = guide_legend(title.position = "top",
                              title = "Years since emergence",
                              title.hjust = 0.5)) +
  theme_minimal(14) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        panel.border = element_rect(linewidth = 0.2, fill = NA),
        aspect.ratio = 1,
        text = element_text(family = "lato"))



m1 <- gam(days_to_max_mg ~ yse, family = poisson(), data = days_to_max_mg)
m2 <- gam(days_to_max_mg ~ yse, family = gaussian(link = "log"), data = days_to_max_mg)
m3 <- gam(days_to_max_mg ~ yse, family = nb(), data = days_to_max_mg)
m4 <- gam(days_to_max_mg ~ yse, family = tw(), data = days_to_max_mg)

bbmle::AICctab(m1, m2, m3, m4)

summary(m1)$dev.exp
summary(m2)$dev.exp
summary(m3)$dev.exp
summary(m4)$dev.exp


ggplot(days_to_max_mg, aes(x = yse, y = days_to_max_mg)) +
  geom_point(size = 4, alpha = 0.5, pch = 21, stroke = 2) +
  geom_smooth(method = "gam", linewidth = 1, color = 2,
              formula = y ~ x, method.args = list(family = nb())) +
  labs(x = "Years since disease emergence",
       y = "Number of days to maximum MG load") +
  theme_minimal(14) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        panel.border = element_rect(linewidth = 0.2, fill = NA),
        aspect.ratio = 1,
        text = element_text(family = "lato"))




library(lmodel2)

m <- lmodel2(days_to_max_mg ~ yse, data = days_to_max_mg)

ggplot(days_to_max_mg, aes(x = yse, y = days_to_max_mg)) +
  geom_point(size = 4, alpha = 0.5, pch = 21, stroke = 2) +
  # mean
  geom_abline(intercept = m$regression.results[3,2], slope = m$regression.results[3,3]) +
  # lower
  geom_abline(intercept = m$confidence.intervals[3,2], slope = m$confidence.intervals[3,4]) +
  # upper
  geom_abline(intercept = m$confidence.intervals[3,3], slope = m$confidence.intervals[3,5]) +
  labs(x = "Years since disease emergence",
       y = "Number of days to maximum MG load") +
  theme_minimal(14) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        panel.border = element_rect(linewidth = 0.2, fill = NA),
        aspect.ratio = 1,
        text = element_text(family = "lato"))



ggplot(df, aes(x = log_mg, y = eyescore_sum, color = yse)) +
  facet_wrap(~bird_id, nrow = 5) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = F, color = 1) +
  scale_color_carto_c(palette = "Temps") +
  labs(y = "Eyescore",
       x = "MG (log)") +
  guides(color = guide_legend(title.position = "top",
                              title = "Years since emergence",
                              title.hjust = 0.5)) +
  theme_minimal(14) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        panel.border = element_rect(linewidth = 0.2, fill = NA),
        aspect.ratio = 1,
        text = element_text(family = "lato"))

