# Replication and duration
# Virulence and duration


# REPLICATION ~ DURATION
# Model
RepDurMod <- gam(replication ~ duration + s(isolate, bs = "re"), 
                 family = tw(),
                 gamma = 1.4,
                 data = df)

# Predictions
RepDurPreds <- get_predictions(RepDurMod)
RepDurPlot <- ggplot(RepDurPreds, aes(x = x, y = mean)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "gray80", alpha = 0.5) +
  geom_line() +
  theme_minimal(14) +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_continuous(breaks = pretty_breaks()) +
  labs(y = "Replication", x = "Duration") +
  theme(aspect.ratio = 1, panel.border = element_rect(linewidth = 0.5, fill = NA, color = 8), 
        panel.grid.minor = element_blank()); RepDurPlot




# VIRULENCE ~ DURATION
# Model
VirDurMod <- gam(virulence ~ duration + s(isolate, bs = "re"), 
                 family = tw(),
                 gamma = 1.4,
                 data = df)

# Predictions
VirDurPreds <- get_predictions(VirDurMod)
VirDurPlot <- ggplot(VirDurPreds, aes(x = x, y = mean)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "gray80", alpha = 0.5) +
  geom_line() +
  theme_minimal(14) +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_continuous(breaks = pretty_breaks()) +
  labs(y = "Virulence", x = "Duration") +
  theme(aspect.ratio = 1, panel.border = element_rect(linewidth = 0.5, fill = NA, color = 8), 
        panel.grid.minor = element_blank()); VirDurPlot

pfinal <- plot_grid(RepDurPlot, VirDurPlot, ncol = 2, align = 'hv'); pfinal

ggsave("VTH_duration_gates.pdf", plot = pfinal, path = "/Users/gatesdupont/Desktop", height = 4.5)


summary(RepDurMod)
summary(VirDurMod)
