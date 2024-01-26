library(tidyverse)
library(readxl)
library(mgcv)
library(gratia)
library(scales)
library(cowplot)

# Function to get predictions from model
get_predictions <- function(fitted.model){
  
  # Make prediction data frame
  x.vals <- fitted.model$model[,2]
  min.x <- min(x.vals)
  max.x <- max(x.vals)
  x.new <- seq(min.x, max.x, length = 100)
  pred.df <- data.frame(
    x = x.new,
    isolate = fitted.model$model$isolate[1]) %>%
    mutate(x2 = x)
  colnames(pred.df)[3] <- colnames(fitted.model$model)[2]
  
  # Predict to new data
  preds <- predict.gam(
    fitted.model, 
    newdata = pred.df, 
    exlucde = "s(isolate)",
    type = "link",
    se = T)
  
  # Record predictions
  pred.df$lwr <- fitted.model$family$linkinv(preds$fit - 1.96 * preds$se.fit)
  pred.df$upr <- fitted.model$family$linkinv(preds$fit + 1.96 * preds$se.fit)
  pred.df$mean <- fitted.model$family$linkinv(preds$fit)
  
  result <- pred.df
  
  return(result)
  
}


# --- Load data ----

df <- read_xlsx("data/Verified-data-for-Gates-25-Jan-2024.xlsx") %>%
  mutate(transmission = rate) %>%
  mutate(isolate = as.factor(isolate)); df


# ---- Models and figures ----

# FITNESS ~ REPLICATION
# Model
FitRepMod <- gam(fitness ~ s(replication) + s(isolate, bs = "re"), 
                 family = tw(),
                 gamma = 1.4,
                 data = df)

# Predictions
FitRepPreds <- get_predictions(FitRepMod)
p1 <- ggplot(FitRepPreds, aes(x = x, y = mean)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "gray80", alpha = 0.5) +
  geom_line() +
  theme_minimal(14) +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_continuous(breaks = pretty_breaks()) +
  labs(y = "Fitness", x = "Replication") +
  theme(aspect.ratio = 1, panel.border = element_rect(linewidth = 0.5, fill = NA, color = 8), 
        panel.grid.minor = element_blank()); p1


# FITNESS ~ TRANSMISSION
# Model
FitTransMod <- gam(fitness ~ transmission + s(isolate, bs = "re"), 
                   family = tw(),
                   gamma = 1.4,
                   data = df)

# Predictions
FitTransPreds <- get_predictions(FitTransMod)
p2 <- ggplot(FitTransPreds, aes(x = x, y = mean)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "gray80", alpha = 0.5) +
  geom_line() +
  theme_minimal(14) +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_continuous(breaks = pretty_breaks()) +
  labs(y = "Fitness", x = "Transmisson") +
  theme(aspect.ratio = 1, panel.border = element_rect(linewidth = 0.5, fill = NA, color = 8), 
        panel.grid.minor = element_blank()); p2


# FITNESS ~ VIRULENCE
# Model
FitVirMod <- gam(fitness ~ s(virulence) + s(isolate, bs = "re"), 
                 family = tw(),
                 gamma = 1.4,
                 data = df)

# Predictions
FitVirPreds <- get_predictions(FitVirMod)
p3 <- ggplot(FitVirPreds, aes(x = x, y = mean)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "gray80", alpha = 0.5) +
  geom_line() +
  theme_minimal(14) +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_continuous(breaks = pretty_breaks()) +
  labs(y = "Fitness", x = "Virulence") +
  theme(aspect.ratio = 1, panel.border = element_rect(linewidth = 0.5, fill = NA, color = 8), 
        panel.grid.minor = element_blank()); p3


# FITNESS ~ DURATION
# Model
FitDurMod <- gam(fitness ~ duration + s(isolate, bs = "re"), 
                 family = tw(),
                 gamma = 1.4,
                 data = df)

# Predictions
FitDurPreds <- get_predictions(FitDurMod)
p4 <- ggplot(FitDurPreds, aes(x = x, y = mean)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "gray80", alpha = 0.5) +
  geom_line() +
  theme_minimal(14) +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_continuous(breaks = pretty_breaks()) +
  labs(y = "Fitness", x = "Duration") +
  theme(aspect.ratio = 1, panel.border = element_rect(linewidth = 0.5, fill = NA, color = 8), 
        panel.grid.minor = element_blank()); p4


# ---- Output -----

pfinal <- plot_grid(p1,p2,p3,p4, ncol = 2, align = 'hv'); pfinal

ggsave("VTH_fitness_gates.pdf", plot = pfinal, path = "/Users/gatesdupont/Desktop")
