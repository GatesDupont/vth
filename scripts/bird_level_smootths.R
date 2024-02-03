
m <- gam(eyescore_sum ~ s(DPI, by = years_since_emergence, bs = "fs") + s(bird_id, bs = "re") + s(strain, bs = "re") + s(sex, bs = "re"), 
         data = df, 
         gamma = 1.4, 
         family = nb())



m <- gam(eyescore_sum ~ s(DPI, by = years_since_emergence, bs = "fs"),
         data = df, 
         gamma = 1.4, 
         family = poisson())

pred_df <- expand.grid(
  sex = sample(unique(df$sex), 1),
  strain = sample(unique(df$strain), 1),
  bird_id = sample(unique(df$bird_id), 1),
  DPI = seq(0,71,0.1),
  years_since_emergence = seq(0,20)) %>%
  as_tibble()

pred_df$eyescore_sum_pred = predict(m, pred_df, type = "response", exclude = c("s(bird_id)","s(strain)","s(sex)"))


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






pred_df <- expand.grid(
  bird_id = unique(df$bird_id),
  DPI = seq(0,71,0.1)) %>%
  as_tibble() %>%
  left_join(
    df2 %>% select(bird_id, strain, years_since_emergence, sex) %>% distinct()) %>%
  mutate(bird_id = factor(bird_id, levels = new_bird_id_order))

pred_df$eyescore_sum_pred = predict(m, pred_df, type = "response")


ggplot() +
  facet_wrap(~ bird_id, nrow = 5) +
  geom_line(data = pred_df, size = 1,
            aes(x = DPI, y = eyescore_sum_pred, color = years_since_emergence, group = years_since_emergence)) +
  geom_point(data = df2, aes(x = DPI, y = eyescore_sum, color = years_since_emergence), na.rm = T) +
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

