# Replication and duration
# Virulence and duration


# REPLICATION ~ DURATION
# Model
DurRepMod <- gam(duration ~ replication + s(isolate, bs = "re"), 
                 family = nb(),
                 gamma = 1.4,
                 data = df)

# Predictions
DurRepPreds <- get_predictions(DurRepMod)
DurRepPlot <- ggplot(DurRepPreds, aes(x = x, y = mean)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "gray80", alpha = 0.5) +
  geom_line() +
  theme_minimal(14) +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_continuous(breaks = pretty_breaks()) +
  labs(y = "Duration", x = "Replication") +
  theme(aspect.ratio = 1, panel.border = element_rect(linewidth = 0.5, fill = NA, color = 8), 
        panel.grid.minor = element_blank()); DurRepPlot




# VIRULENCE ~ DURATION
# Model
DurVirMod <- gam(duration ~ virulence + s(isolate, bs = "re"), 
                 family = nb(),
                 gamma = 1.4,
                 data = df)

# Predictions
DurVirPreds <- get_predictions(DurVirMod)
DurVirPlot <- ggplot(DurVirPreds, aes(x = x, y = mean)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "gray80", alpha = 0.5) +
  geom_line() +
  theme_minimal(14) +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_continuous(breaks = pretty_breaks()) +
  labs(y = "Duration", x = "Virulence") +
  theme(aspect.ratio = 1, panel.border = element_rect(linewidth = 0.5, fill = NA, color = 8), 
        panel.grid.minor = element_blank()); DurVirPlot

pfinal <- plot_grid(DurRepPlot, DurVirPlot, ncol = 2, align = 'hv'); pfinal

ggsave("VTH_DurRep_DurVir_gates.pdf", plot = pfinal, path = "/Users/gatesdupont/Desktop", height = 4.5)


summary(DurRepMod)
summary(DurVirMod)
