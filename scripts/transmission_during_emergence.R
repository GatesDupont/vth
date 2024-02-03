library(tidyverse)
library(ggtext)

df <- read.csv("data/beta.csv") %>%
  mutate(geo = if_else(geo == "east", "East", "West")) %>%
  mutate(first_year = if_else(geo == "East", 1993, 2005)) %>%
  mutate(years_since_emergence = year - first_year)

m <- glm(mean ~ years_since_emergence, data = df, family = gaussian(link = "log"))

newdf <- data.frame(
  years_since_emergence = seq(0, max(df$years_since_emergence), by = 0.1)) %>%
  mutate(pred = predict(m, newdata = ., type = "response"))

# Transmission should be on the log axis
fig <- ggplot(df, aes(x = years_since_emergence, y = mean)) +
  
  geom_line(data = newdf, aes(x = years_since_emergence, y = pred)) +
  
  geom_pointrange(aes(ymin = lwr, ymax = upr, color = geo)) +
  geom_point(aes(color = geo), size = 3) +
  
  scale_color_manual(values = c("East" = "blue", "West" = "red")) +
  
  labs(x = "Years since emergence", 
       title = "Model: log(beta) = -3.5030 + years * 0.1068",
       y = "Transmission rate parameter beta (1/week)",
       caption = str_wrap("Figure: Estimates of transmission from Ellner based on aviary experiments. Note that intervals are 95% confidence intervals from Ellner's results. Although we show the uncertainty in the point estimates here, we ignore that uncertainty when fitting the exponential model. We use 1993 and 2005 as the year of emergenece in the east and west, respectively."),
       color = "Strain") +
  
  theme_minimal(10) +
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5, size = 11),
        plot.caption = element_text(hjust = 0),
        text = element_text(family = "Lato"),
        legend.title = element_text(size = 8),
        panel.grid.minor = element_blank()); fig

ggsave(plot = fig, path = "/Users/gatesdupont/Desktop", filename = "beta_change.pdf", 
       width = 6.5, height = 8, units = "in", device = cairo_pdf)




dflog <- df %>%
  mutate(lwr = log(lwr),
         mean = log(mean),
         upr = log(upr)) %>%
  mutate(sd = (upr - mean)/1.96)

simulated_data <- list()

for(i in 1:nrow(dflog)){
  
  # Simulate 1000 numbers from the normal distribution
  simulated_data[[i]] <- data.frame(
    years_since_emergence = dflog$years_since_emergence[i],
    geo = dflog$geo[i],
    sim_mean = exp(rnorm(n = 1000, mean = dflog$mean[i], sd = dflog$sd[i])))
  
}
    
simdat <- do.call(rbind, simulated_data) %>%
  as_tibble()

m2 <- glm(sim_mean ~ years_since_emergence, data = simdat, family = gaussian(link = "log"))

newdf2 <- data.frame(
  years_since_emergence = seq(0, max(df$years_since_emergence), by = 0.1)) %>%
  mutate(pred = predict(m2, newdata = ., type = "response"))


# Transmission should be on the log axis
fig <- ggplot(df, aes(x = years_since_emergence, y = mean)) +
  
  geom_line(data = newdf, aes(x = years_since_emergence, y = pred)) +
  geom_line(data = newdf2, aes(x = years_since_emergence, y = pred), color = 3, linewidth = 1) +
  
  geom_pointrange(aes(ymin = lwr, ymax = upr, color = geo)) +
  geom_point(aes(color = geo), size = 3) +
  
  scale_color_manual(values = c("East" = "blue", "West" = "red")) +
  
  labs(x = "Years since emergence", 
       title = "Model: log(beta) = -3.5030 + years * 0.1068",
       y = "Transmission rate parameter beta (1/week)",
       caption = str_wrap("Figure: Estimates of transmission from Ellner based on aviary experiments. Note that intervals are 95% confidence intervals from Ellner's results. Although we show the uncertainty in the point estimates here, we ignore that uncertainty when fitting the exponential model. We use 1993 and 2005 as the year of emergenece in the east and west, respectively."),
       color = "Strain") +
  
  theme_minimal(10) +
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5, size = 11),
        plot.caption = element_text(hjust = 0),
        text = element_text(family = "Lato"),
        legend.title = element_text(size = 8),
        panel.grid.minor = element_blank()); fig





       