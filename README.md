# Master-thesis

Moving CUSUM Window for category FC026 (Minced meat plain) 

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

- function run_backtest_daily runs the analyze_critical_day function for each day with available data for the specified customer and product, and compiles a performance table summarizing the results.
  output: Gives for every day the number of points above the traditional M limit (green colored), if they are above the CUSUM limit and how many censequent previous points are above the CUSUM limit (early warning)
