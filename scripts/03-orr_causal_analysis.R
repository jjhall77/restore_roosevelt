#==============================================================================
# Operation Restore Roosevelt - Causal Analysis (Final Version)
# 
# Primary Method: gsynth (generalized synthetic control with factor model)
# Robustness: Synthetic DiD, Augmented SCM
# 
# Unit of analysis: Census Tracts (donors) vs Zone (rosie_buffer aggregate)
# Note: Zone is ~5-10x larger than individual tracts, so methods that
#       require level-matching (pure SC) will fail. gsynth handles this
#       via interactive fixed effects.
#==============================================================================

library(here)
library(tidyverse)
library(lubridate)
library(sf)
library(gsynth)
library(synthdid)
library(augsynth)
library(panelView)
library(future)
library(furrr)
library(patchwork)

# Set up parallel processing
n_workers <- min(parallel::detectCores() - 1, 8)
plan(multisession, workers = n_workers)
cat(sprintf("Using %d parallel workers\n", n_workers))

# Create output directory
dir.create(here("output", "causal_tract"), showWarnings = FALSE, recursive = TRUE)

#==============================================================================
# 1. PARAMETERS AND DATA LOADING
#==============================================================================

cat("\n=== SETUP ===\n")

# Load EDA objects
if (file.exists(here("output", "eda_objects.rds"))) {
  eda_objects <- readRDS(here("output", "eda_objects.rds"))
  list2env(eda_objects, envir = .GlobalEnv)
} else {
  stop("Run EDA scripts first to generate eda_objects.rds")
}

# Key dates
op_start <- as.Date("2024-10-15")
op_end <- op_start + 89
treatment_month <- floor_date(op_start, "month")  # 2024-10-01

# Analysis window
data_start <- as.Date("2022-01-01")
data_end <- as.Date("2025-09-01")  # Last complete month

all_months <- seq(data_start, data_end, by = "month")
n_months <- length(all_months)
n_pre <- sum(all_months < treatment_month)
n_post <- n_months - n_pre

cat(sprintf("Analysis window: %s to %s\n", data_start, data_end))
cat(sprintf("Treatment month: %s\n", treatment_month))
cat(sprintf("Periods: %d total (%d pre, %d post)\n", n_months, n_pre, n_post))

#==============================================================================
# 2. IDENTIFY ZONE AND DONOR TRACTS
#==============================================================================

cat("\n=== IDENTIFYING TREATMENT AND DONOR UNITS ===\n")

# Tracts that intersect the zone = excluded from donors
zone_tracts <- nyct %>%
  st_filter(rosie_buffer, .predicate = st_intersects) %>%
  pull(boro_ct2020)

cat(sprintf("Tracts intersecting zone (excluded): %d\n", length(zone_tracts)))

# Donor pool: all tracts NOT intersecting zone
donor_tracts <- nyct %>%
  filter(!boro_ct2020 %in% zone_tracts) %>%
  pull(boro_ct2020)

cat(sprintf("Potential donor tracts: %d\n", length(donor_tracts)))

#==============================================================================
# 3. BUILD PANEL DATA
#==============================================================================

cat("\n=== BUILDING PANEL DATA ===\n")

# Monthly street violent crime by tract (for donors)
street_violent_tract <- street_violent_crime %>%
  filter(date >= data_start) %>%
  st_join(nyct %>% select(boro_ct2020)) %>%
  st_drop_geometry() %>%
  filter(!is.na(boro_ct2020)) %>%
  mutate(month = floor_date(date, "month")) %>%
  filter(month <= data_end) %>%
  count(boro_ct2020, month, name = "crime_count")

# Zone aggregate: crimes within rosie_buffer (not tract-based)
zone_monthly <- street_violent_crime %>%
  filter(date >= data_start) %>%
  st_filter(rosie_buffer, .predicate = st_intersects) %>%
  st_drop_geometry() %>%
  mutate(month = floor_date(date, "month")) %>%
  filter(month <= data_end) %>%
  count(month, name = "crime_count") %>%
  right_join(tibble(month = all_months), by = "month") %>%
  replace_na(list(crime_count = 0)) %>%
  mutate(boro_ct2020 = "ZONE") %>%
  arrange(month)

cat(sprintf("\nZone monthly crime: mean=%.1f, sd=%.1f, range=%d-%d\n",
            mean(zone_monthly$crime_count),
            sd(zone_monthly$crime_count),
            min(zone_monthly$crime_count),
            max(zone_monthly$crime_count)))

# Complete donor panel
donor_panel <- expand_grid(
  boro_ct2020 = donor_tracts,
  month = all_months
) %>%
  left_join(street_violent_tract, by = c("boro_ct2020", "month")) %>%
  replace_na(list(crime_count = 0))

# Donor statistics
donor_stats <- donor_panel %>%
  filter(month < treatment_month) %>%
  group_by(boro_ct2020) %>%
  summarize(
    mean_crime = mean(crime_count),
    total_crime = sum(crime_count),
    nonzero_months = sum(crime_count > 0),
    .groups = "drop"
  )

cat(sprintf("Donor tract monthly crime: mean=%.1f, median=%.1f, max=%.1f\n",
            mean(donor_stats$mean_crime),
            median(donor_stats$mean_crime),
            max(donor_stats$mean_crime)))

#==============================================================================
# 4. FILTER DONOR POOL
#==============================================================================

cat("\n=== FILTERING DONOR POOL ===\n")

# Remove tracts with very little crime (too sparse for estimation)
# Keep tracts with at least some activity
min_total_crime <- 10
min_nonzero_months <- n_pre * 0.3  # At least 30% non-zero months

good_donors <- donor_stats %>%
  filter(
    total_crime >= min_total_crime,
    nonzero_months >= min_nonzero_months
  ) %>%
  pull(boro_ct2020)

cat(sprintf("Donors after filtering: %d (from %d)\n", 
            length(good_donors), length(donor_tracts)))
cat(sprintf("  Min total crime: %d\n", min_total_crime))
cat(sprintf("  Min non-zero months: %.0f (%.0f%%)\n", 
            min_nonzero_months, min_nonzero_months/n_pre*100))

#==============================================================================
# 5. BUILD FINAL PANEL (gsynth format)
#==============================================================================

cat("\n=== BUILDING FINAL PANEL ===\n")

# Combine zone and donors
panel_final <- donor_panel %>%
  filter(boro_ct2020 %in% good_donors) %>%
  mutate(tag = "donor") %>%
  bind_rows(
    zone_monthly %>% mutate(tag = "treated")
  ) %>%
  mutate(
    treat = as.integer(tag == "treated" & month >= treatment_month),
    log_crime = log(crime_count + 0.5),
    post = month >= treatment_month
  ) %>%
  arrange(boro_ct2020, month)

# Verify balanced panel
n_units <- n_distinct(panel_final$boro_ct2020)
n_periods <- n_distinct(panel_final$month)

cat(sprintf("Final panel: %d units × %d periods = %d rows\n",
            n_units, n_periods, nrow(panel_final)))

# Check balance
balance_check <- panel_final %>% count(boro_ct2020) %>% pull(n) %>% unique()
if (length(balance_check) == 1) {
  cat("✓ Panel is balanced\n")
} else {
  warning("Panel is UNBALANCED")
}

#==============================================================================
# 6. VISUALIZE DATA
#==============================================================================

cat("\n=== VISUALIZING DATA ===\n")

# Zone time series
p_zone <- ggplot(zone_monthly, aes(x = month, y = crime_count)) +
  geom_line(linewidth = 1, color = "darkred") +
  geom_point(color = "darkred", size = 1.5) +
  geom_vline(xintercept = treatment_month, linetype = "dashed") +
  annotate("rect", xmin = treatment_month, xmax = max(all_months),
           ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "red") +
  labs(
    title = "Zone Monthly Street Violence",
    subtitle = "Crimes within Roosevelt Ave corridor (rosie_buffer)",
    x = "Month", y = "Crime Count"
  ) +
  theme_minimal()

ggsave(here("output", "causal_tract", "00_zone_time_series.png"), p_zone,
       width = 10, height = 6, dpi = 300)

# Zone vs donor distribution
donor_monthly_avg <- panel_final %>%
  filter(tag == "donor") %>%
  group_by(month) %>%
  summarize(
    mean = mean(crime_count),
    median = median(crime_count),
    p25 = quantile(crime_count, 0.25),
    p75 = quantile(crime_count, 0.75),
    p90 = quantile(crime_count, 0.90),
    .groups = "drop"
  )

p_comparison <- ggplot() +
  geom_ribbon(data = donor_monthly_avg,
              aes(x = month, ymin = p25, ymax = p75),
              alpha = 0.2, fill = "steelblue") +
  geom_line(data = donor_monthly_avg,
            aes(x = month, y = median, color = "Donor Median"),
            linewidth = 0.8) +
  geom_line(data = donor_monthly_avg,
            aes(x = month, y = p90, color = "Donor 90th %ile"),
            linewidth = 0.6, linetype = "dashed") +
  geom_line(data = zone_monthly,
            aes(x = month, y = crime_count, color = "Zone"),
            linewidth = 1.2) +
  geom_vline(xintercept = treatment_month, linetype = "dashed") +
  scale_color_manual(values = c("Zone" = "darkred", 
                                "Donor Median" = "steelblue",
                                "Donor 90th %ile" = "gray50")) +
  labs(
    title = "Zone vs Donor Tracts",
    subtitle = "Zone is aggregated corridor; donors are individual tracts (note scale difference)",
    x = "Month", y = "Crime Count", color = ""
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(here("output", "causal_tract", "01_zone_vs_donors.png"), p_comparison,
       width = 10, height = 6, dpi = 300)

# Treatment pattern visualization
png(here("output", "causal_tract", "02_treatment_pattern.png"),
    width = 12, height = 8, units = "in", res = 300)
panelview(crime_count ~ treat, data = panel_final, 
          index = c("boro_ct2020", "month"),
          type = "treat", main = "Treatment Pattern",
          by.timing = TRUE)
dev.off()

#==============================================================================
# 7. PRIMARY METHOD: GSYNTH
#==============================================================================

cat("\n=== GSYNTH (PRIMARY METHOD) ===\n")
cat("Using generalized synthetic control with interactive fixed effects.\n")
cat("This method handles level differences between zone and individual tracts.\n\n")

results <- list()

# --- Raw counts ---
cat("Fitting gsynth (raw counts)...\n")

tryCatch({
  gs_raw <- gsynth(
    crime_count ~ treat,
    data      = panel_final,
    index     = c("boro_ct2020", "month"),
    force     = "two-way",      # Unit + time fixed effects
    CV        = TRUE,           # Cross-validate number of factors
    r         = c(0, 5),        # Search 0-5 factors
    se        = TRUE,
    inference = "parametric",
    nboots    = 500,
    parallel  = TRUE,
    seed      = 42
  )
  
  results$gs_raw <- gs_raw
  
  cat(sprintf("\nGsynth Results (raw):\n"))
  cat(sprintf("  Factors selected: %d\n", gs_raw$r.cv))
  cat(sprintf("  Average ATT: %.2f (SE: %.2f)\n", gs_raw$att.avg, sqrt(gs_raw$att.avg.var)))
  
  # Counterfactual plot
  png(here("output", "causal_tract", "03_gsynth_raw_counterfactual.png"),
      width = 10, height = 6, units = "in", res = 300)
  plot(gs_raw, type = "counterfactual", raw = "none",
       main = "Gsynth: Treated vs Synthetic Control (Raw Counts)")
  dev.off()
  
  # Gap plot
  png(here("output", "causal_tract", "04_gsynth_raw_gap.png"),
      width = 10, height = 6, units = "in", res = 300)
  plot(gs_raw, type = "gap", main = "Gsynth: Treatment Effect Over Time")
  dev.off()
  
  # CT plot (average treatment effect with CI)
  png(here("output", "causal_tract", "05_gsynth_raw_att.png"),
      width = 10, height = 6, units = "in", res = 300)
  plot(gs_raw, type = "ct", main = "Gsynth: Average Treatment Effect")
  dev.off()
  
}, error = function(e) {
  cat(sprintf("Gsynth (raw) failed: %s\n", e$message))
})

# --- Log counts ---
cat("\nFitting gsynth (log counts)...\n")

tryCatch({
  gs_log <- gsynth(
    log_crime ~ treat,
    data      = panel_final,
    index     = c("boro_ct2020", "month"),
    force     = "two-way",
    CV        = TRUE,
    r         = c(0, 5),
    se        = TRUE,
    inference = "parametric",
    nboots    = 500,
    parallel  = TRUE,
    seed      = 42
  )
  
  results$gs_log <- gs_log
  
  cat(sprintf("\nGsynth Results (log):\n"))
  cat(sprintf("  Factors selected: %d\n", gs_log$r.cv))
  cat(sprintf("  Average ATT: %.3f (SE: %.3f)\n", gs_log$att.avg, sqrt(gs_log$att.avg.var)))
  cat(sprintf("  Approx %% change: %.1f%%\n", (exp(gs_log$att.avg) - 1) * 100))
  
  png(here("output", "causal_tract", "06_gsynth_log_counterfactual.png"),
      width = 10, height = 6, units = "in", res = 300)
  plot(gs_log, type = "counterfactual", raw = "none",
       main = "Gsynth: Treated vs Synthetic Control (Log Counts)")
  dev.off()
  
  png(here("output", "causal_tract", "07_gsynth_log_gap.png"),
      width = 10, height = 6, units = "in", res = 300)
  plot(gs_log, type = "gap", main = "Gsynth: Treatment Effect Over Time (Log)")
  dev.off()
  
}, error = function(e) {
  cat(sprintf("Gsynth (log) failed: %s\n", e$message))
})

#==============================================================================
# 8. ROBUSTNESS: SYNTHETIC DIFF-IN-DIFF
#==============================================================================

cat("\n=== SYNTHETIC DIFF-IN-DIFF (ROBUSTNESS) ===\n")
cat("Note: SDID uses differencing which mitigates level mismatch.\n")
cat("The plot may look poor but the estimate can still be valid.\n\n")

# Prepare matrices for synthdid
# Select subset of donors for computational feasibility
n_donors_sdid <- min(200, length(good_donors))

# Select by pre-treatment correlation with zone
zone_pre <- zone_monthly %>%
  filter(month < treatment_month) %>%
  arrange(month) %>%
  pull(crime_count)

donor_cors <- panel_final %>%
  filter(tag == "donor", month < treatment_month) %>%
  group_by(boro_ct2020) %>%
  arrange(month) %>%
  summarize(
    series = list(crime_count),
    .groups = "drop"
  ) %>%
  mutate(
    cor_zone = map_dbl(series, ~cor(.x, zone_pre, use = "complete.obs"))
  ) %>%
  filter(!is.na(cor_zone)) %>%
  slice_max(cor_zone, n = n_donors_sdid)

sdid_donors <- donor_cors$boro_ct2020

cat(sprintf("Using %d donors for SDID (top correlated)\n", length(sdid_donors)))

# Build matrix
panel_sdid <- panel_final %>%
  filter(boro_ct2020 == "ZONE" | boro_ct2020 %in% sdid_donors)

panel_wide_raw <- panel_sdid %>%
  select(boro_ct2020, month, crime_count) %>%
  pivot_wider(names_from = month, values_from = crime_count, values_fill = 0) %>%
  column_to_rownames("boro_ct2020")

panel_wide_log <- panel_sdid %>%
  select(boro_ct2020, month, log_crime) %>%
  pivot_wider(names_from = month, values_from = log_crime, values_fill = log(0.5)) %>%
  column_to_rownames("boro_ct2020")

# Reorder: controls first, treated last
prepare_matrix <- function(wide_df) {
  mat <- as.matrix(wide_df)
  zone_row <- mat["ZONE", , drop = FALSE]
  control_rows <- mat[rownames(mat) != "ZONE", , drop = FALSE]
  rbind(control_rows, zone_row)
}

Y_raw <- prepare_matrix(panel_wide_raw)
Y_log <- prepare_matrix(panel_wide_log)

N0 <- nrow(Y_raw) - 1
T0 <- sum(as.Date(colnames(Y_raw)) < treatment_month)

cat(sprintf("Matrix: %d controls × %d periods (T0=%d)\n", N0, ncol(Y_raw), T0))

# --- SDID Raw ---
cat("\nFitting SDID (raw counts)...\n")

tryCatch({
  sdid_raw <- synthdid_estimate(Y_raw, N0, T0)
  se_raw <- sqrt(vcov(sdid_raw, method = "placebo"))
  
  results$sdid_raw <- sdid_raw
  results$sdid_raw_se <- se_raw
  
  cat(sprintf("SDID estimate (raw): %.2f (SE: %.2f)\n", sdid_raw, se_raw))
  cat(sprintf("95%% CI: [%.2f, %.2f]\n", 
              sdid_raw - 1.96*se_raw, sdid_raw + 1.96*se_raw))
  
  png(here("output", "causal_tract", "08_sdid_raw.png"),
      width = 10, height = 6, units = "in", res = 300)
  synthdid_plot(sdid_raw, se.method = "placebo")
  dev.off()
  
}, error = function(e) {
  cat(sprintf("SDID (raw) failed: %s\n", e$message))
})

# --- SDID Log ---
cat("\nFitting SDID (log counts)...\n")

tryCatch({
  sdid_log <- synthdid_estimate(Y_log, N0, T0)
  se_log <- sqrt(vcov(sdid_log, method = "placebo"))
  
  results$sdid_log <- sdid_log
  results$sdid_log_se <- se_log
  
  cat(sprintf("SDID estimate (log): %.3f (SE: %.3f)\n", sdid_log, se_log))
  cat(sprintf("Approx %% change: %.1f%%\n", (exp(sdid_log) - 1) * 100))
  
  png(here("output", "causal_tract", "09_sdid_log.png"),
      width = 10, height = 6, units = "in", res = 300)
  synthdid_plot(sdid_log, se.method = "placebo")
  dev.off()
  
}, error = function(e) {
  cat(sprintf("SDID (log) failed: %s\n", e$message))
})

# --- Standard SC and DiD for comparison ---
cat("\nFitting SC and DiD for comparison...\n")

tryCatch({
  sc_raw <- sc_estimate(Y_raw, N0, T0)
  did_raw <- did_estimate(Y_raw, N0, T0)
  
  results$sc_raw <- sc_raw
  results$did_raw <- did_raw
  
  cat(sprintf("SC estimate (raw): %.2f\n", sc_raw))
  cat(sprintf("DiD estimate (raw): %.2f\n", did_raw))
  
}, error = function(e) {
  cat(sprintf("SC/DiD failed: %s\n", e$message))
})

#==============================================================================
# 9. ROBUSTNESS: AUGMENTED SYNTHETIC CONTROL
#==============================================================================

cat("\n=== AUGMENTED SYNTHETIC CONTROL (ROBUSTNESS) ===\n")
cat("Ridge augmentation corrects for imperfect pre-treatment fit.\n\n")

# Use same donor subset as SDID
augsynth_data <- panel_sdid %>%
  mutate(
    unit = boro_ct2020,
    time = as.integer(factor(month))
  )

# --- ASCM Raw ---
cat("Fitting ASCM (raw counts)...\n")

tryCatch({
  ascm_raw <- augsynth(
    crime_count ~ treat,
    unit = unit,
    time = time,
    data = augsynth_data,
    progfunc = "Ridge",
    scm = TRUE,
    CV = FALSE,
    lambda = 1e-3
  )
  
  results$ascm_raw <- ascm_raw
  ascm_raw_summ <- summary(ascm_raw)
  
  cat(sprintf("ASCM estimate (raw): %.2f (SE: %.2f)\n",
              ascm_raw_summ$average_att$Estimate,
              ascm_raw_summ$average_att$Std.Error))
  
  p_ascm_raw <- plot(ascm_raw) +
    labs(title = "Augmented SCM (Raw Counts)") +
    theme_minimal()
  
  ggsave(here("output", "causal_tract", "10_ascm_raw.png"), p_ascm_raw,
         width = 10, height = 6, dpi = 300)
  
}, error = function(e) {
  cat(sprintf("ASCM (raw) failed: %s\n", e$message))
})

# --- ASCM Log ---
cat("\nFitting ASCM (log counts)...\n")

tryCatch({
  ascm_log <- augsynth(
    log_crime ~ treat,
    unit = unit,
    time = time,
    data = augsynth_data,
    progfunc = "Ridge",
    scm = TRUE,
    CV = FALSE,
    lambda = 1e-3
  )
  
  results$ascm_log <- ascm_log
  ascm_log_summ <- summary(ascm_log)
  
  cat(sprintf("ASCM estimate (log): %.3f (SE: %.3f)\n",
              ascm_log_summ$average_att$Estimate,
              ascm_log_summ$average_att$Std.Error))
  
  p_ascm_log <- plot(ascm_log) +
    labs(title = "Augmented SCM (Log Counts)") +
    theme_minimal()
  
  ggsave(here("output", "causal_tract", "11_ascm_log.png"), p_ascm_log,
         width = 10, height = 6, dpi = 300)
  
}, error = function(e) {
  cat(sprintf("ASCM (log) failed: %s\n", e$message))
})

#==============================================================================
# 10. EXTRACT AND COMPARE ESTIMATES
#==============================================================================

cat("\n=== RESULTS SUMMARY ===\n")

# Build results table
results_table <- tibble(
  
  Method = character(),
  Scale = character(),
  Estimate = numeric(),
  SE = numeric(),
  CI_lower = numeric(),
  CI_upper = numeric(),
  Pct_Change = numeric()
)

# Gsynth raw
if (!is.null(results$gs_raw)) {
  se <- sqrt(results$gs_raw$att.avg.var)
  est <- results$gs_raw$att.avg
  results_table <- results_table %>%
    add_row(Method = "Gsynth", Scale = "Raw",
            Estimate = est, SE = se,
            CI_lower = est - 1.96*se, CI_upper = est + 1.96*se,
            Pct_Change = NA)
}

# Gsynth log
if (!is.null(results$gs_log)) {
  se <- sqrt(results$gs_log$att.avg.var)
  est <- results$gs_log$att.avg
  results_table <- results_table %>%
    add_row(Method = "Gsynth", Scale = "Log",
            Estimate = est, SE = se,
            CI_lower = est - 1.96*se, CI_upper = est + 1.96*se,
            Pct_Change = (exp(est) - 1) * 100)
}

# SDID raw
if (!is.null(results$sdid_raw)) {
  est <- as.numeric(results$sdid_raw)
  se <- results$sdid_raw_se
  results_table <- results_table %>%
    add_row(Method = "Synthetic DiD", Scale = "Raw",
            Estimate = est, SE = se,
            CI_lower = est - 1.96*se, CI_upper = est + 1.96*se,
            Pct_Change = NA)
}

# SDID log
if (!is.null(results$sdid_log)) {
  est <- as.numeric(results$sdid_log)
  se <- results$sdid_log_se
  results_table <- results_table %>%
    add_row(Method = "Synthetic DiD", Scale = "Log",
            Estimate = est, SE = se,
            CI_lower = est - 1.96*se, CI_upper = est + 1.96*se,
            Pct_Change = (exp(est) - 1) * 100)
}

# ASCM raw
if (!is.null(results$ascm_raw)) {
  summ <- summary(results$ascm_raw)
  est <- summ$average_att$Estimate
  se <- summ$average_att$Std.Error
  results_table <- results_table %>%
    add_row(Method = "Augmented SCM", Scale = "Raw",
            Estimate = est, SE = se,
            CI_lower = est - 1.96*se, CI_upper = est + 1.96*se,
            Pct_Change = NA)
}

# ASCM log
if (!is.null(results$ascm_log)) {
  summ <- summary(results$ascm_log)
  est <- summ$average_att$Estimate
  se <- summ$average_att$Std.Error
  results_table <- results_table %>%
    add_row(Method = "Augmented SCM", Scale = "Log",
            Estimate = est, SE = se,
            CI_lower = est - 1.96*se, CI_upper = est + 1.96*se,
            Pct_Change = (exp(est) - 1) * 100)
}

# SC and DiD (no SE)
if (!is.null(results$sc_raw)) {
  results_table <- results_table %>%
    add_row(Method = "Synthetic Control", Scale = "Raw",
            Estimate = as.numeric(results$sc_raw), SE = NA,
            CI_lower = NA, CI_upper = NA, Pct_Change = NA)
}

if (!is.null(results$did_raw)) {
  results_table <- results_table %>%
    add_row(Method = "Diff-in-Diff", Scale = "Raw",
            Estimate = as.numeric(results$did_raw), SE = NA,
            CI_lower = NA, CI_upper = NA, Pct_Change = NA)
}

# Print results
cat("\n")
print(results_table, n = 20)

# Save results
write_csv(results_table, here("output", "causal_tract", "results_summary.csv"))

#==============================================================================
# 11. FOREST PLOT
#==============================================================================

cat("\n=== GENERATING PLOTS ===\n")

# Forest plot - Raw scale
p_forest_raw <- results_table %>%
  filter(Scale == "Raw", !is.na(SE)) %>%
  ggplot(aes(x = Estimate, y = reorder(Method, Estimate))) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title = "Treatment Effect Estimates (Raw Scale)",
    subtitle = "Point estimates with 95% CI",
    x = "Estimated Effect (change in monthly crime count)", 
    y = ""
  ) +
  theme_minimal()

ggsave(here("output", "causal_tract", "12_forest_plot_raw.png"), p_forest_raw,
       width = 10, height = 5, dpi = 300)

# Forest plot - Log scale
p_forest_log <- results_table %>%
  filter(Scale == "Log", !is.na(SE)) %>%
  ggplot(aes(x = Estimate, y = reorder(Method, Estimate))) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title = "Treatment Effect Estimates (Log Scale)",
    subtitle = sprintf("Point estimates with 95%% CI. Negative = crime reduction."),
    x = "Estimated Effect (log scale)", 
    y = ""
  ) +
  theme_minimal()

ggsave(here("output", "causal_tract", "13_forest_plot_log.png"), p_forest_log,
       width = 10, height = 5, dpi = 300)

#==============================================================================
# 12. PRE-TREATMENT FIT DIAGNOSTICS (GSYNTH)
#==============================================================================

cat("\n=== PRE-TREATMENT FIT DIAGNOSTICS ===\n")

if (!is.null(results$gs_raw)) {
  
  # Extract treated and counterfactual series
  Y_treated <- results$gs_raw$Y.tr
  Y_counterfactual <- results$gs_raw$Y.ct
  
  # Time index
  time_index <- as.Date(rownames(results$gs_raw$Y.tr))
  
  fit_df <- tibble(
    month = time_index,
    treated = as.vector(Y_treated),
    counterfactual = as.vector(Y_counterfactual),
    gap = treated - counterfactual,
    period = if_else(month < treatment_month, "Pre-treatment", "Post-treatment")
  )
  
  # Pre-treatment fit statistics
  pre_fit <- fit_df %>% filter(period == "Pre-treatment")
  
  rmse_pre <- sqrt(mean(pre_fit$gap^2))
  mae_pre <- mean(abs(pre_fit$gap))
  cor_pre <- cor(pre_fit$treated, pre_fit$counterfactual)
  
  cat(sprintf("\nGsynth Pre-treatment Fit:\n"))
  cat(sprintf("  RMSE: %.2f\n", rmse_pre))
  cat(sprintf("  MAE: %.2f\n", mae_pre))
  cat(sprintf("  Correlation: %.3f\n", cor_pre))
  
  # Fit plot
  p_fit <- ggplot(fit_df, aes(x = month)) +
    geom_line(aes(y = treated, color = "Treated (Observed)"), linewidth = 1) +
    geom_line(aes(y = counterfactual, color = "Synthetic Control"), 
              linewidth = 1, linetype = "dashed") +
    geom_vline(xintercept = treatment_month, linetype = "dashed", color = "black") +
    annotate("rect", xmin = min(fit_df$month), xmax = treatment_month,
             ymin = -Inf, ymax = Inf, alpha = 0.05, fill = "blue") +
    annotate("text", x = treatment_month - months(6), y = max(fit_df$treated) * 0.95,
             label = sprintf("Pre-treatment\nRMSE=%.1f, r=%.2f", rmse_pre, cor_pre),
             hjust = 1, size = 3) +
    scale_color_manual(values = c("Treated (Observed)" = "darkred", 
                                  "Synthetic Control" = "steelblue")) +
    labs(
      title = "Gsynth: Pre-treatment Fit Diagnostic",
      subtitle = "Blue shading = pre-treatment period used for matching",
      x = "Month", y = "Crime Count", color = ""
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  ggsave(here("output", "causal_tract", "14_pretreatment_fit.png"), p_fit,
         width = 10, height = 6, dpi = 300)
  
  # Gap plot with CI
  p_gap_detail <- ggplot(fit_df, aes(x = month, y = gap)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_line(linewidth = 0.8) +
    geom_point(aes(color = period), size = 2) +
    geom_vline(xintercept = treatment_month, linetype = "dashed") +
    scale_color_manual(values = c("Pre-treatment" = "gray50", "Post-treatment" = "darkred")) +
    labs(
      title = "Treatment Effect Over Time (Actual - Counterfactual)",
      subtitle = "Negative values = crime reduction relative to synthetic",
      x = "Month", y = "Gap (Treated - Synthetic)", color = ""
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  ggsave(here("output", "causal_tract", "15_gap_over_time.png"), p_gap_detail,
         width = 10, height = 6, dpi = 300)
  
  # Save fit data
  write_csv(fit_df, here("output", "causal_tract", "gsynth_fit_data.csv"))
}

#==============================================================================
# 13. INTERPRETATION GUIDANCE
#==============================================================================

cat("\n")
cat("================================================================\n")
cat("INTERPRETATION GUIDANCE\n")
cat("================================================================\n")

cat("
PRIMARY RESULTS (Gsynth):
- Gsynth uses interactive fixed effects which absorb level differences
  between the aggregated zone and individual tract donors.
- The counterfactual plot should show good pre-treatment fit.
- Negative ATT = crime reduction in zone relative to counterfactual.

ROBUSTNESS CHECKS:
- SDID: Uses differencing which mitigates level mismatch. The visual
  plot may look poor but the estimate can still be valid.
- ASCM: Ridge augmentation corrects for imperfect SCM fit.
- If all methods agree on direction, that's reassuring.
- If they disagree, investigate why and report with caveats.

CAVEATS:
- Pre-intervention spike may bias results (regression to mean).
- Zone is geographically unique (el train + commercial corridor).
- No perfect donor exists; we're constructing a weighted counterfactual.

")

#==============================================================================
# 14. SAVE ALL RESULTS
#==============================================================================

saveRDS(results, here("output", "causal_tract", "all_results.rds"))
saveRDS(panel_final, here("output", "causal_tract", "panel_final.rds"))

# Clean up parallel
plan(sequential)

cat("\n================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("================================================================\n")
cat(sprintf("\nOutputs saved to: %s\n", here("output", "causal_tract")))

list.files(here("output", "causal_tract")) %>%
  paste(" -", .) %>%
  cat(sep = "\n")

cat("\n")