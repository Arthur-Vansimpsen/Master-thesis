#function for normal distribution with left censoring

simulation_compare <- function(censoring_proportion) {
  
  # Data generation and substitution methods
  LOD = 1
  sigma = 1
  mu = LOD - qnorm(censoring_proportion)*sigma
  
  data = c(
    rnorm(200, mu, sigma)
  )
  
  y_obs <- ifelse(data < LOD, LOD, data)
  y_obs_sub = ifelse(data < LOD, LOD/sqrt(2) , data)
  censored <- data < LOD

  
  # True mean and SD
  true_mu = mu
  true_sigma = sigma
  
  # Normal MLE
  fit_norm = survreg(Surv(y_obs, !censored, type = "left") ~ 1, dist = "gaussian")
  
  # Normal ROS
  ros_fit = ros(y_obs, censored, forwardT = NULL, reverseT = NULL)

  # Kaplan–Meier
  km_fit <- survfit(Surv(data, !censored) ~ 1)
  
  t <- c(0, km_fit$time)
  S <- c(1, km_fit$surv)
  
  mu_km <- sum(diff(t) * S[-length(S)])
  
  # Mean and SD estimates of all methods
  mu_LOD = mean(y_obs)
  sigma_LOD = sd(y_obs)
  
  mu_LODSQRT2 = mean(y_obs_sub)
  sigma_LODSQRT2 = sd(y_obs_sub)
  
  mu_norm = fit_norm$coefficients[1]
  sigma_norm = fit_norm$scale
  
  mu_ros = mean(ros_fit$modeled)
  sigma_ros = sd(ros_fit$modeled)

  return(c(true_mu, true_sigma,
           mu_LOD, sigma_LOD,
           mu_LODSQRT2, sigma_LODSQRT2,
           mu_norm, sigma_norm,
           mu_ros, sigma_ros,
           mu_km
  ))
}

# Set censoring levels to evaluate
censoring_levels = c(0.00001, 0.25, 0.5, 0.75, 0.85, 0.95)

# Run 1000 simulations for each censoring level and methods
results_list = lapply(censoring_levels, function(cens_level) {
  
  res = replicate(1000, simulation_compare(censoring_proportion = cens_level))
  
  res_df = as.data.frame(t(res))
  colnames(res_df) = c(
    "true_mu", "true_sigma",
    "mu_LOD", "sigma_LOD",
    "mu_LODSQRT2", "sigma_LODSQRT2",
    "mu_norm", "sigma_norm",
    "mu_ros", "sigma_ros",
    "mu_km")
  
  res_df$censoring = cens_level
  
  return(res_df)
})

results_all = do.call(rbind, results_list)

# Bias
Bias <- aggregate(
  cbind(
    mu_LOD = mu_LOD - true_mu,
    mu_LODSQRT2 = mu_LODSQRT2 - true_mu,
    mu_norm = mu_norm - true_mu,
    mu_ros = mu_ros - true_mu,
    mu_km = mu_km - true_mu,
    
    sigma_LOD = sigma_LOD - true_sigma,
    sigma_LODSQRT2 = sigma_LODSQRT2 - true_sigma,
    sigma_norm = sigma_norm - true_sigma,
    sigma_ros = sigma_ros - true_sigma
  ) ~ censoring,
  data = results_all,
  FUN = mean
)
# RMSE
rmse = aggregate(
  cbind(mu_LOD = (mu_LOD - true_mu)^2,
        mu_LODSQRT2 = (mu_LODSQRT2 - true_mu)^2,
        mu_norm = (mu_norm - true_mu)^2,
        mu_ros = (mu_ros - true_mu)^2,
        mu_km = (mu_km - true_mu)^2,
        
        sigma_LOD = (sigma_LOD - true_sigma)^2,
        sigma_LODSQRT2 = (sigma_LODSQRT2 - true_sigma)^2,
        sigma_norm = (sigma_norm - true_sigma)^2,
        sigma_ros = (sigma_ros - true_sigma)^2
  ) ~ censoring, data = results_all, FUN = function(x) sqrt(mean(x)))

library(tidyr)
library(dplyr)

#Plots
plot_df <- Bias %>%
  pivot_longer(
    cols = -c(censoring),
    names_to = c("param", "method"),
    names_sep = "_",
    values_to = "bias"
  ) %>%
  mutate(
    method = recode(method,
                    "LOD" = "Substitution (LOD)",
                    "LODSQRT2" = "Substitution (LOD/√2)",
                    "norm" = "Normal MLE",
                    "ros" = "ROS",
                    "km" = "Kaplan–Meier"
    )
  )
ggplot(plot_df, aes(x = censoring, y = bias, color = method)) +
  
  geom_line(aes(linetype = param), linewidth = 1) +
  geom_point(aes(shape = param), size = 2) +
  
  scale_linetype_manual(values = c(mu = "solid", sigma = "dashed")) +
  scale_shape_manual(values = c(mu = 16, sigma = 17)) +  
  
  theme_minimal(base_size = 14) +
  
  theme(
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  ) +
  
  labs(
    title = "B: Line plot bias",
    x = "Censoring level",
    y = "bias",
    color = "Method",
    linetype = "Parameter",
    shape = "Parameter"
  )


#Plots
plot_df <- rmse %>%
  pivot_longer(
    cols = -c(censoring),
    names_to = c("param", "method"),
    names_sep = "_",
    values_to = "rmse"
  ) %>%
  mutate(
    method = recode(method,
                    "LOD" = "Substitution (LOD)",
                    "LODSQRT2" = "Substitution (LOD/√2)",
                    "norm" = "Normal MLE",
                    "ros" = "ROS",
                    "km" = "Kaplan–Meier"
    )
  )
ggplot(plot_df, aes(x = censoring, y = rmse, color = method)) +
  
  geom_line(aes(linetype = param), linewidth = 1) +
  geom_point(aes(shape = param), size = 2) +
  
  scale_linetype_manual(values = c(mu = "solid", sigma = "dashed")) +
  scale_shape_manual(values = c(mu = 16, sigma = 17)) +  
  
  theme_minimal(base_size = 14) +
  
  theme(
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  ) +
  
  labs(
    title = "B: Line plot RMSE",
    x = "Censoring level",
    y = "RMSE",
    color = "Method",
    linetype = "Parameter",
    shape = "Parameter"
  )
ggsave("bias_big.png", width = 12, height = 8)
ggsave("RMSE_big.png", width = 12, height = 8)

