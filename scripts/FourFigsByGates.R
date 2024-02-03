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
    exclude = "s(isolate)",
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


# ---- Make models and figures ----

# VIRULENCE ~ REPLICATION

# Model
VirRepMod <- gam(virulence ~ s(replication) + s(isolate, bs = "re"), 
                 family = nb(),
                 gamma = 1.4,
                 data = df)

# Predictions
VirRepPreds <- get_predictions(VirRepMod)
p1 <- ggplot(VirRepPreds, aes(x = x, y = mean)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "gray80", alpha = 0.5) +
  geom_line() +
  theme_minimal(14) +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_continuous(breaks = pretty_breaks()) +
  labs(y = "Virulence", x = "Replication") +
  theme(aspect.ratio = 1, panel.border = element_rect(linewidth = 0.5, fill = NA, color = 8), 
        panel.grid.minor = element_blank()); p1


# TRANSMISSION ~ REPLICATION

# Model
TransRepMod <- gam(transmission ~ s(replication) + s(isolate, bs = "re"), 
                   gamma = 1.4,
                   data = df)

# Predictions
TransRepPreds <- get_predictions(TransRepMod)
p2 <- ggplot(TransRepPreds, aes(x = x, y = mean)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "gray80", alpha = 0.5) +
  geom_line() +
  theme_minimal(14) +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_continuous(breaks = pretty_breaks(), labels = scales::label_scientific()) +
  labs(y = "Transmission", x = "Replication") +
  theme(aspect.ratio = 1, panel.border = element_rect(linewidth = 0.5, fill = NA, color = 8), 
        panel.grid.minor = element_blank()); p2


# TRANSMISSION ~ VIRULENCE

# Model
TransVirMod <- gam(transmission ~ s(virulence) + s(isolate, bs = "re"), 
                   gamma = 1.4,
                   data = df)

# Predictions
TransVirPreds <- get_predictions(TransVirMod)
p3 <- ggplot(TransVirPreds, aes(x = x, y = mean)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "gray80", alpha = 0.5) +
  geom_line() +
  theme_minimal(14) +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_continuous(breaks = pretty_breaks(), labels = scales::label_scientific()) +
  labs(y = "Transmission", x = "Virulence") +
  theme(aspect.ratio = 1, panel.border = element_rect(linewidth = 0.5, fill = NA, color = 8), 
        panel.grid.minor = element_blank()); p3


# VIRULENCE ~ DURATION

# Model
VirDurMod <- gam(virulence ~ duration + s(isolate, bs = "re"), 
                 family = nb(),
                 gamma = 1.4,
                 data = df)

# Predictions
VirDurPreds <- get_predictions(VirDurMod)
p4 <- ggplot(VirDurPreds, aes(x = x, y = mean)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "gray80", alpha = 0.5) +
  geom_line() +
  theme_minimal(14) +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_continuous(breaks = pretty_breaks()) +
  labs(y = "Virulence", x = "Duration") +
  theme(aspect.ratio = 1, panel.border = element_rect(linewidth = 0.5, fill = NA, color = 8), 
        panel.grid.minor = element_blank()); p4


# ---- Make main figure ----

pfinal <- plot_grid(p1,p2,p3,p4, ncol = 2, align = 'hv'); pfinal

ggsave("VTH_gates.pdf", plot = pfinal, path = "/Users/gatesdupont/Desktop")







  

  



