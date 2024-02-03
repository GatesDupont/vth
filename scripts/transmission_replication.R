# TRANSMISSION ~ REPLICATION

# Model
TransRepMod <- gam(transmission ~ replication + s(isolate, bs = "re"), 
                   gamma = 1.4,
                   data = df)

# Predictions
TransRepPreds <- get_predictions(TransRepMod)
p2 <- ggplot(TransRepPreds, aes(x = x, y = mean)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "gray80", alpha = 0.5) +
  geom_line() +
  theme_minimal(14) +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_continuous(breaks = pretty_breaks()) +
  labs(y = "Transmission", x = "Replication") +
  theme(aspect.ratio = 1, panel.border = element_rect(linewidth = 0.5, fill = NA, color = 8), 
        panel.grid.minor = element_blank()); p2


ggplot(df, aes(x = replication, y = transmission)) +
  geom_point(aes(color = isolate), size = 4, pch = 1, stroke = 3) +
  geom_smooth(method = "gam",
              formula = y ~ s(x),
              method.args = list(
                gamma = 1.4,
                family = gaussian(link = "log"))) +
  theme_minimal()


newdat <- data.frame(replication = seq(0, 20, by = 0.1), isolate = "NA")

newdat$pred <- predict(TransRepMod, newdata = newdat, type = "response", exclude = "s(isolate)")

ggplot(newdat, aes(x = replication, y = pred)) +
  geom_line()

ggplot(df, aes(x = replication, y = transmission)) + 
  geom_point(aes(color = isolate), size = 4, pch = 1, stroke = 3) +
  geom_ribbon(data = TransRepPreds, aes(x = x, y = mean, ymin = lwr, ymax = upr), fill = "gray80", alpha = 0.5) +
  geom_line(data = TransRepPreds, aes(x = x, y = mean)) +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_continuous(breaks = pretty_breaks()) +
  labs(y = "Transmission", x = "Replication") +
  theme(aspect.ratio = 1, panel.border = element_rect(linewidth = 0.5, fill = NA, color = 8), 
        panel.grid.minor = element_blank())



ggplot(df, aes(x = replication, y = transmission)) + 
  geom_point(aes(color = isolate), size = 4, pch = 1, stroke = 3) +
  geom_line(data = newdat, aes(x = replication, y = pred)) +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_continuous(breaks = pretty_breaks()) +
  labs(y = "Transmission", x = "Replication") +
  theme(aspect.ratio = 1, panel.border = element_rect(linewidth = 0.5, fill = NA, color = 8), 
        panel.grid.minor = element_blank())
