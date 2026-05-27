library(NADA)
library(survival)
library(ggplot2)
library(sn)
library(qcr)
library(dplyr)

# Function to generate data with left censoring from a normal distribution
generate_data <- function(n, mu, sigma, LOD, delta = 0) {
  

  data <- rnorm(n, mu + delta, sigma)
  censored <- data < LOD
  y_obs <- ifelse(censored, LOD, data)
  
  list(y = y_obs, censored = censored)
}


#Estimation of mean and SD in Phase I using MLE or ROS.
estimate_phase1 <- function(y, censored, estimation = "MLE") {
  
  df <- data.frame(y = y, censored = censored)
  
  if (estimation == "MLE") {
    
    fit <- survreg(Surv(y, !censored, type = "left") ~ 1,
                   data = df,
                   dist = "gaussian")
    
    mu_hat <- fit$coefficients[1]
    sigma_hat <- fit$scale
    
  } else if (estimation == "ROS") {
    
    ros_fit <- ros(df$y, df$censored)
    
    vals <- ros_fit$modeled
    vals <- vals[is.finite(vals)]
    
    if (length(vals) < 5) {
      return(list(mu = NA, sigma = NA))
    }
    
    mu_hat <- mean(vals)
    sigma_hat <- sd(vals)
    
    if (is.na(sigma_hat) || sigma_hat <= 0) {
      return(list(mu = NA, sigma = NA))
    }
  }
  
  list(mu = mu_hat, sigma = sigma_hat)
}

# Function of different imputation methods for Phase II 
impute_data <- function(y, censored, method, mu, sigma, LOD) {
  
  if (method == "LOD") {
    y[censored] <- LOD
  }
  
  if (method == "LOD_SQRT2") {
    y[censored] <- LOD / sqrt(2)
  }
  
  if (method == "RANDOM") {
    
    z <- (LOD - mu) / sigma
    
    # draw from truncated normal below LOD
    u <- runif(sum(censored), min = 0, max = pnorm(z))
    
    y[censored] <- mu + sigma * qnorm(u)
  }
  
  return(y)
}

# Run CUSUM on imputed data and return the index of the first signal
run_cusum <- function(y, mu, sigma, k, h) {
  
  if (is.na(sigma) || sigma <= 0) return(NA)
  
  monitor <- data.frame(
    indicator = y,
    Index = seq_along(y)
  )
  
  qcs_monitor <- qcs.cusum(
    x = monitor,
    var.index = "indicator",
    sample.index = "Index",
    plot = FALSE,
    decision.interval = h,
    se.shift = k,
    center = mu,
    std.dev = sigma
  )
  
  signal_idx <- qcs_monitor$violations$upper[1]
  
  if (is.na(signal_idx)) {
    return(length(y))
  } else {
    return(signal_idx)
  }
}

#Generate data for Phase I, Phase II and the shift in phase II. Substitute censored observations when not running the LH CUSUM.
simulate_run <- function(n_phase1, n_phase2, n_phase3,
                         mu, sigma, LOD,
                         method, k, h,
                         delta = 0,
                         estimation) {
  
  # ---- Phase I ----
  d1 <- generate_data(n_phase1, mu, sigma, LOD)
  est <- estimate_phase1(d1$y, d1$censored, estimation)
  
  if (is.na(est$mu) || is.na(est$sigma) || est$sigma <= 0) {
    return(NA)
  }
  
  # ---- Phase II ----
  d2 <- generate_data(n_phase2, mu, sigma, LOD)
  
  y2 <- impute_data(d2$y, d2$censored,
                    method,
                    est$mu, est$sigma,
                    LOD)
  
  #Introduce shift
  delta_abs <- delta * sigma
  
  d3 <- generate_data(n_phase3, mu, sigma, LOD, delta_abs)
  
  y3 = impute_data(d3$y, d3$censored,
                   method,
                   est$mu, est$sigma,
                   LOD)
  y_all <- c(d2$y, d3$y)
  censored_all <- c(d2$censored, d3$censored)
  
  if (method == "LIK") {
    
    return(run_lik_cusum(
      y_all,
      censored_all,
      est$mu,
      est$sigma,
      k,
      h
    ))
    
  } else {
    
    y_imp <- impute_data(y_all, censored_all,
                         method,
                         est$mu, est$sigma,
                         LOD)
    
    return(run_cusum(y_imp, est$mu, est$sigma, k, h))
  }
}

#Likelihood scores per observations for LH CUSUM and return the index of the first signal
lik_cusum <- function(y, censored, mu0, sigma, delta, h) {
  
  mu1 <- mu0 + delta
  
  S <- 0
  S_path <- numeric(length(y))  # store cumulative sum
  
  for (t in seq_along(y)) {
    
    if (!censored[t]) {
      
      l0 <- dnorm(y[t], mu0, sigma, log = TRUE)
      l1 <- dnorm(y[t], mu1, sigma, log = TRUE)
      
    } else {
      
      l0 <- pnorm(y[t], mu0, sigma, log.p = TRUE)
      l1 <- pnorm(y[t], mu1, sigma, log.p = TRUE)
    }
    
    S <- max(0, S + (l1 - l0))
    S_path[t] <- S
    
    if (S > h) {
      return(list(
        signal = t,
        S = S_path
      ))
    }
  }
  
  return(list(
    signal = length(y),
    S = S_path
  ))
}

# LH CUSUM function
run_lik_cusum <- function(y, censored, mu, sigma, k, h) {
  delta <- 2 * k/2 * sigma
  lik_cusum(y, censored, mu, sigma, delta, h)$signal
}


#Computation of ARL0 and ARL1 scores
compute_arl <- function(nsim, ...) {
  
  runs <- numeric(nsim)
  
  for (i in 1:nsim) {
    val <- simulate_run(...)
    if (is.na(val)) next
    runs[i] <- val
  }
  
  q <- quantile(runs, probs = c(0.1, 0.25, 0.75, 0.9), na.rm = TRUE)
  
  c(
    ARL = mean(runs),
    ARL_med = median(runs),
    SDRL = sd(runs),
    q10 = q[[1]],
    q25 = q[[2]],
    q75 = q[[3]],
    q90 = q[[4]]
    
  )
}

# Simulation parameters. Use very small value when simulating the 0% censoring scenario to avoid issues with the qnorm function.ns
methods <- c("LOD", "LOD_SQRT2", "RANDOM", "LIK")
censoring_levels <- c(0.0001, 0.25, 0.5, 0.85, 0.95)
deltas <- c(0, 0.5, 1)  # 0 = ARL0
estimation <- c("MLE")

results <- list()

#Run simulation
for (cens in censoring_levels) {
  
  LOD <- 1
  sigma <- 1
  mu <- LOD - qnorm(cens) * sigma
  
  for (method in methods) {
    
    for (delta in deltas) {
      
      arl <- compute_arl(
        nsim = 200,
        n_phase1 = 200,
        n_phase2 = 25,
        n_phase3 = 3000,
        mu = mu,
        sigma = sigma,
        LOD = LOD,
        method = method,
        k = 1,
        h = 5,
        delta = delta,
        estimation = estimation
      )
      
      results[[length(results) + 1]] <- data.frame(
        method = method,
        censoring = cens,
        k = 1,
        h = 5,
        delta = delta,
        ARL = as.numeric(arl["ARL"]),
        ARL_med = as.numeric(arl["ARL_med"]),
        SDRL = as.numeric(arl["SDRL"]),
        q10 = as.numeric(arl["q10"]),
        q25 = as.numeric(arl["q25"]),
        q75 = as.numeric(arl["q75"]),
        q90 = as.numeric(arl["q90"])
      )
    }
  }
}

final_results <- do.call(rbind, results)
MLE_normal3000  = final_results


#Plots
plot_data <- MLE_normal3000 %>%
  dplyr::select(method, censoring, delta, ARL, q25, q75) %>%
  filter(
    censoring %in% c(0.0001, 0.5, 0.85),
    method %in% c("LOD", "LOD_SQRT2", "RANDOM", "LIK"),
    delta %in% c(0, 0.5, 1)
  ) %>%
  mutate(
    censoring = ifelse(censoring == 0.0001, 0, censoring),
    censoring = factor(censoring, labels = c("0% censored", "50% censored", "85% censored")),
    method = recode(method,
                    "LOD" = "Substitution (LOD)",
                    "LOD_SQRT2" = "Substitution (LOD/√2)",
                    "LIK" = "Likelihood-based CUSUM",
                    "RANDOM" = "Stochastic Random Imputation"
    ),
                    # 🔧 correction step
    across(
      c(ARL, q25, q75),
      ~ ifelse(delta %in% c(0.5, 1), .x - 25, .x)
    ),
    across(
      c(ARL, q25, q75),
      ~ ifelse(.x < 0, NA, .x)
    ))



ggplot(plot_data, aes(x = delta, y = ARL, color = method)) +
  
  geom_line(linewidth = 1.5) +
  geom_point(size = 2) +
  
  geom_errorbar(
    aes(ymin = q25, ymax = q75),
    width = 0.2,
    alpha = 0.7,
    linewidth = 1.0
  ) +
  
  facet_wrap(~ censoring) +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  
  theme_minimal(base_size = 14) +
  
  labs(
    x = "Delta",
    y = "ARL",
    color = "Method",
    title = "ARL across different levels of censoring"
  )


