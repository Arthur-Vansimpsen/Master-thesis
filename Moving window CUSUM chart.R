#------------------------------------------------------------#
# Moving CUSUM Window for category FC026 (Minced meat plain) #
#------------------------------------------------------------#

" Long section of comments

- Important tests in this category are the following:

 Counts:
  +   Aerobe mesophile Koloniezahl_KbE/g
  +   Enterobacteriaceae_KbE/g
  +   Listeria spp._Nachweis
  +   Escherichia coli_KbE/g
  +   Koagulase-positive Staphylokokken_KbE/g
  +   Milchsäurebakterien_KbE/g
  +   Hefen_KbE/g
  +   Listeria spp._KbE/g
  +   Listeria monocytogenes_KbE/g
  +   Präsumtive Pseudomonas_KbE/g
  
Attribute (presence/absence):  
  +   Thermophile Campylobacter_Nachweis
  +   Listeria monocytogenes_Nachweis
  +   Salmonella_Nachweis

- Function analyze_critical_day performs the following steps:
  1. Filter data for the specified customer, product, and indicator.
  2. Identify the first and last observation of the target day.
  3. Define baseline (observations 250 till 50 before target day) and monitoring periods (observations 50 till last observation target date).
  4. Perform Phase I outlier removal on the baseline data using individual Shewhart chart with a 3-sigma limit.
  5. Compute final control limits from the cleaned baseline.
  6. Apply  individual Shewhart chart and CUSUM to the monitoring data with the final limits and specified se_shift.
  7. Identify critical points in the monitoring data where the indicator exceeds the critical value (M).
  8. Plot the CUSUM chart and highlight critical points in green and points from the target day in red.
  Comment: change plot = TRUE to false if you don't want to see the plots.
"

analyze_critical_day = function(customer_nr,
                                 product_name,
                                 indicator,
                                 critical_value,
                                 se_shift,
                                 day,
                                 month,
                                 year,
                                 baseline_n_start = 250,
                                 baseline_n_end   = 50) {
  
  target_date = as.Date(sprintf("%04d-%02d-%02d", year, month, day))
  
  # Filter customer/product
  product_data = data_unique_or_na %>%
    filter(
      CustomerID == customer_nr,
      Productname_QI == product_name
    ) %>%
    arrange(ReceivedDate)
  
  # Filter indicator and remove NA
  product_data_indicator = product_data %>%
    dplyr::select(SampleID, all_of(indicator), ReceivedDate) %>%
    filter(!is.na(.data[[indicator]])) %>%
    arrange(ReceivedDate)
  
  # Compute critical flag
  product_data_indicator = product_data_indicator %>%
    mutate(
      indicator_critical = if_else(.data[[indicator]] >= critical_value, 1, 0)
    )
  
  colnames(product_data_indicator)[2] = "indicator"
  product_data_indicator$indicator = as.numeric(product_data_indicator$indicator)
  product_data_indicator$Index = seq_len(nrow(product_data_indicator))
  
  # Identify first and last observation of target day 
  day_idx = which(product_data_indicator$ReceivedDate == target_date)
  
  if (length(day_idx) == 0) {
    return(NULL)
  }
  
  first_day_idx = min(day_idx)
  last_day_idx  = max(day_idx)

  # Define baseline and monitoring indices
  baseline_start_idx = max(1, first_day_idx - baseline_n_start)
  baseline_end_idx   = max(1, first_day_idx - baseline_n_end)
  
  monitor_start_idx  = max(1, first_day_idx - baseline_n_end + 1)
  monitor_end_idx    = last_day_idx
  
  baseline_data = product_data_indicator[baseline_start_idx:baseline_end_idx, ]
  monitor_data  = product_data_indicator[monitor_start_idx:monitor_end_idx, ]
  
  if (nrow(baseline_data) < 10 | nrow(monitor_data) < 5) {
    return(NULL)
  }
  
  baseline_data$Index = seq_len(nrow(baseline_data))
  monitor_data$Index  = seq_len(nrow(monitor_data))
  
  # Phase I outlier removal on baseline
  data_clean = baseline_data
  iteration = 1
  removed_rows = list() 
  
  repeat {
    data_clean$Index = seq_len(nrow(data_clean))
    
    qcs_baseline = qcs.one(
      x = data_clean,
      var.index = "indicator",
      sample.index = "Index",
      conf.nsigma = 3,
      plot = FALSE
    )
    
    beyond_idx = qcs_baseline$violations$beyond.limits
    
    if (is.null(beyond_idx) || length(beyond_idx) == 0) break
    
    removed_rows[[iteration]] = beyond_idx 
    data_clean = data_clean[-beyond_idx, , drop = FALSE] 
    iteration = iteration + 1
  }
  
  data_clean$Index = seq_len(nrow(data_clean)) 
  
  final_qcs = qcs.one(
    x = data_clean,
    var.index = "indicator",
    sample.index = "Index",
    conf.nsigma = 3,
    plot = FALSE,
    data.name = paste0(
      "Baseline individual Shewhart - Customer ", customer_nr,
      ", Product ", product_name,
      "\nTarget date: ", target_date
    )
  )
  
  monitor_qcs = qcs.one(
    x = monitor_data,
    var.index = "indicator",
    sample.index = "Index",
    conf.nsigma = 3,
    limits = final_qcs$limits,
    center = final_qcs$center,
    plot = FALSE,
    data.name = paste0(
      "Monitoring individual Shewhart - Customer ", customer_nr,
      ", Product ", product_name,
      "\nTarget date: ", target_date
    )
  )
  
  # Apply CUSUM to monitoring
  qcs_monitor = qcs.cusum(
    x = monitor_data,
    var.index = "indicator",
    sample.index = "Index",
    plot = FALSE,
    decision.interval = 3,
    se.shift = se_shift,
    center = final_qcs$center,
    std.dev = final_qcs$std.dev,
    data.name = paste0(
      "CUSUM - Customer ", customer_nr,
      ", Product ", product_name,
      "\nTarget date: ", target_date
    )
  )
  
  # Identify critical points in monitoring
  idx_list = list(
    critical = which(monitor_data$indicator_critical == 1)
  )
  
  cusum_vals = qcs_monitor$pos
  h = qcs_monitor$decision.interval
  
  # Count above/below limits
  green_idx = which(
    monitor_data$indicator_critical == 1 &
      monitor_data$ReceivedDate == target_date
  )
  above_limit = if (length(green_idx) > 0) sum(cusum_vals[green_idx] > h) else 0 
  below_limit = if (length(green_idx) > 0) sum(cusum_vals[green_idx] < -h) else 0 
  
  # Plot critical points on CUSUM in green
  if (length(green_idx) > 0) {
    points(
      monitor_data$Index[green_idx],
      cusum_vals[green_idx],
      col = "green",
      pch = 19,
      cex = 1.2
    )
  }
  
  # COLOR POINTS OF LAST DAY IN RED
  last_day_idx_in_monitor = which(monitor_data$ReceivedDate == target_date)
  if (length(last_day_idx_in_monitor) > 0) {
    points(
      monitor_data$Index[last_day_idx_in_monitor],
      cusum_vals[last_day_idx_in_monitor],
      col = "red",
      pch = 19,
      cex = 0.7
    )
  }
  
  legend (
    "topleft",
    legend = c("Critical points (green)", "Points from target day (red)"),
    col = c("green", "red"),
    pch = 19,
    cex = 0.8
  )
  
  if (length(green_idx) > 0) {
    
    first_green = min(green_idx)
    
    if (first_green > 1) {
      
      before_vals = cusum_vals[1:(first_green - 1)]
      
      above_h = before_vals > h
      
      run_length = 0
      current_run = 0
      
      for (val in above_h) {
        if (val) {
          current_run = current_run + 1
          run_length = current_run
        } else {
          current_run = 0
        }
      }
      
    } else {
      run_length = 0
    }
    
  } else {
    run_length = NA
  }
  
  return(list(
    removed_indices = removed_rows,
    above_limit = above_limit,
    below_limit = below_limit, 
    n_green = length(green_idx),
    number_of_points_above_CUSUM_limit = run_length
  ))
}

#EXAMPLE#
analyze_critical_day(
  customer_nr = 77,
  product_name = "RINDERHACKFLEISCH",
  indicator = "Aerobe mesophile Koloniezahl_KbE/g",
  critical_value = 6.7,
  se_shift = 0.5,
  day = 6,
  month = 11,
  year = 2020
)

filter_day_data = function(data, customer_nr, product_name) {
  
  data %>%
    filter(CustomerID == customer_nr,
           Productname_QI == product_name) %>%
    distinct(ReceivedDate) %>%
    arrange(ReceivedDate) %>%
    mutate(
      year  = as.integer(format(ReceivedDate, "%Y")),
      month = as.integer(format(ReceivedDate, "%m")),
      day   = as.integer(format(ReceivedDate, "%d"))
    )
}

run_backtest_daily = function(customer_nr,
                               product_name,
                               indicator,
                               critical_value,
                               se_shift) {
  
  grid = filter_day_data(data_wide,
                       customer_nr,
                       product_name)
  
  results = list()
  
  for (i in seq_len(nrow(grid))) {
    
    y = grid$year[i]
    m = grid$month[i]
    d = grid$day[i]
    
    cat("Running:", y, "-", m, "-", d, "\n")
    
    res = tryCatch({
      
      analyze_critical_day(
        customer_nr = customer_nr,
        product_name = product_name,
        indicator = indicator,
        critical_value = critical_value,
        se_shift = se_shift,
        day = d,
        month = m,
        year = y
      )
      
    }, error = function(e) NULL)
    
    if (!is.null(res) && is.list(res)) {
      
      n_green = res$n_green
      n_above = res$above_limit
      
      pct = if (n_green > 0) {
        (n_above / n_green) * 100
      } else {
        NA
      }
      
      pct_label = if (!is.na(pct)) {
        paste0(round(pct, 1), "%")
      } else {
        NA
      }
      
      results[[i]] = data.frame(
        customer = customer_nr,
        date = as.Date(sprintf("%04d-%02d-%02d", y, m, d)),
        n_green = n_green,
        n_green_above_limit = n_above,
        pct_green_above = pct_label,
        number_of_points_above_CUSUM_limit = res$number_of_points_above_CUSUM_limit
      )
    }
  }
  
  dplyr::bind_rows(results)
}

#EXAMPLE#
performance_table = run_backtest_daily(
  customer_nr = 77,
  product_name = "RINDERHACKFLEISCH",
  indicator = "Aerobe mesophile Koloniezahl_KbE/g",
  critical_value = 6.7,
  se_shift = 0.5
)

#print rows for which n_green > 0
performance_table[performance_table$n_green > 0, ]
