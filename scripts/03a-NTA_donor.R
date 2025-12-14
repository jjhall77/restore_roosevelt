#03a- causal analysis for NTA's
#==============================================================================
# Operation Restore Roosevelt - Causal Analysis (NTA Version)
# 
# Primary Method: gsynth (generalized synthetic control with factor model)
# Robustness: Synthetic DiD, Augmented SCM
# 
# Unit of analysis: NTAs (donors) vs Zone (rosie_buffer aggregate)
# NTAs are larger than tracts, providing better scale match with zone.
#==============================================================================

library(here)
library(tidyverse)
library(lubridate)
library(sf)
library(gsynth)
library(synthdid)
library(augsynth)
library(future)
library(furrr)
library(patchwork)

# Set up parallel processing
n_workers <- min(parallel::detectCores() - 1, 8)
plan(multisession, workers = n_workers)
cat(sprintf("Using %d parallel workers\n", n_workers))

# Create output directory
dir.create(here("output", "causal_nta"), showWarnings = FALSE, recursive = TRUE)

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

# Load NTA shapefile if not already loaded
# Adjust path as needed - check what you have available
if (!exists("nynta")) {
  # Try common locations
  if (file.exists(here("data", "nynta.shp"))) {
    nynta <- st_read(here("data", "nynta.shp")) %>% st_transform(2263)
  } else if (file.exists(here("data", "nynta2020.shp"))) {
    nynta <- st_read(here("data", "nynta2020.shp")) %>% st_transform(2263)
  } else {
    stop("NTA shapefile not found. Please load nynta manually.")
  }
}

# Check NTA ID column name - adapt as needed
nta_id_col <- intersect(c("ntacode", "nta2020", "ntaname", "NTACode", "NTA2020"), names(nynta))[1]
if (is.na(nta_id_col)) {
  
  cat("Available NTA columns:", paste(names(nynta), collapse = ", "), "\n")
  stop("Cannot find NTA ID column. Please set nta_id_col manually.")
}
cat(sprintf("Using NTA ID column: %s\n", nta_id_col))

# Rename for consistency
nynta <- nynta %>% rename(nta_id = all_of(nta_id_col))

# Key dates
op_start <- as.Date("2024-10-15")
op_end <- op_start + 89
treatment_month <- floor_date(op_start, "month")  # 2024-10-01

# Analysis window
data_start <- as.Date("2022-01-01")
data_end <- as.Date("2025-09-30")  # Last complete month

all_months <- seq(data_start, data_end, by = "month")
n_months <- length(all_months)
n_pre <- sum(all_months < treatment_month)
n_post <- n_months - n_pre

cat(sprintf("Analysis window: %s to %s\n", data_start, data_end))
cat(sprintf("Treatment month: %s\n", treatment_month))
cat(sprintf("Periods: %d total (%d pre, %d post)\n", n_months, n_pre, n_post))
cat(sprintf("Total NTAs in NYC: %d\n", nrow(nynta)))

#==============================================================================
# 2. IDENTIFY ZONE AND DONOR NTAs
#==============================================================================

cat("\n=== IDENTIFYING TREATMENT AND DONOR UNITS ===\n")

# NTAs that intersect the zone = excluded from donors
zone_ntas <- nynta %>%
  st_filter(rosie_buffer, .predicate = st_intersects) %>%
  pull(nta_id)

cat(sprintf("NTAs intersecting zone (excluded): %d\n", length(zone_ntas)))
cat("  ", paste(zone_ntas, collapse = ", "), "\n")

# Donor pool: all NTAs NOT intersecting zone
donor_ntas <- nynta %>%
  filter(!nta_id %in% zone_ntas) %>%
  pull(nta_id)

cat(sprintf("Potential donor NTAs: %d\n", length(donor_ntas)))

#==============================================================================
# 3. BUILD PANEL DATA
#==============================================================================

cat("\n=== BUILDING PANEL DATA ===\n")

# Monthly street violent crime by NTA (for donors)
street_violent_nta <- street_violent_crime %>%
  filter(date >= data_start) %>%
  st_join(nynta %>% select(nta_id)) %>%
  st_drop_geometry() %>%
  filter(!is.na(nta_id)) %>%
  mutate(month = floor_date(date, "month")) %>%
  filter(month <= data_end) %>%
  count(nta_id, month, name = "crime_count")

# Zone aggregate: crimes within rosie_buffer (not NTA-based)
zone_monthly <- street_violent_crime %>%
  filter(date >= data_start) %>%
  st_filter(rosie_buffer, .predicate = st_intersects) %>%
  st_drop_geometry() %>%
  mutate(month = floor_date(date, "month")) %>%
  filter(month <= data_end) %>%
  count(month, name = "crime_count") %>%
  right_join(tibble(month = all_months), by = "month") %>%
  replace_na(list(crime_count = 0)) %>%
  mutate(nta_id = "ZONE") %>%
  arrange(month)

cat(sprintf("\nZone monthly crime: mean=%.1f, sd=%.1f, range=%d-%d\n",
            mean(zone_monthly$crime_count),
            sd(zone_monthly$crime_count),
            min(zone_monthly$crime_count),
            max(zone_monthly$crime_count)))

# Complete donor panel
donor_panel <- expand_grid(
  nta_id = donor_ntas,
  month = all_months
) %>%
  left_join(street_violent_nta, by = c("nta_id", "month")) %>%
  replace_na(list(crime_count = 0))

# Donor statistics
donor_stats <- donor_panel %>%
  filter(month < treatment_month) %>%
  group_by(nta_id) %>%
  summarize(
    mean_crime = mean(crime_count),
    total_crime = sum(crime_count),
    nonzero_months = sum(crime_count > 0),
    .groups = "drop"
  )

cat(sprintf("Donor NTA monthly crime: mean=%.1f, median=%.1f, max=%.1f\n",
            mean(donor_stats$mean_crime),
            median(donor_stats$mean_crime),
            max(donor_stats$mean_crime)))

# Compare scale
cat(sprintf("\nScale comparison:\n"))
cat(sprintf("  Zone mean: %.1f crimes/month\n", mean(zone_monthly$crime_count)))
cat(sprintf("  Donor NTA mean: %.1f crimes/month\n", mean(donor_stats$mean_crime)))
cat(sprintf("  Ratio: %.1fx\n", mean(zone_monthly$crime_count) / mean(donor_stats$mean_crime)))

#==============================================================================
# 4. FILTER DONOR POOL
#==============================================================================

cat("\n=== FILTERING DONOR POOL ===\n")

# For NTAs, we can be less aggressive with filtering since they're larger
min_total_crime <- 20
min_nonzero_months <- n_pre * 0.5  # At least 50% non-zero months

good_donors <- donor_stats %>%
  filter(
    total_crime >= min_total_crime,
    nonzero_months >= min_nonzero_months
  ) %>%
  pull(nta_id)

cat(sprintf("Donors after filtering: %d (from %d)\n", 
            length(good_donors), length(donor_ntas)))
cat(sprintf("  Min total crime: %d\n", min_total_crime))
cat(sprintf("  Min non-zero months: %.0f (%.0f%%)\n", 
            min_nonzero_months, min_nonzero_months/n_pre*100))

#==============================================================================
# 5. BUILD FINAL PANEL (gsynth format)
#==============================================================================

cat("\n=== BUILDING FINAL PANEL ===\n")

# Combine zone and donors
panel_final <- donor_panel %>%
  filter(nta_id %in% good_donors) %>%
  mutate(tag = "donor") %>%
  bind_rows(
    zone_monthly %>% mutate(tag = "treated")
  ) %>%
  mutate(
    treat = as.integer(tag == "treated" & month >= treatment_month),
    log_crime = log(crime_count + 0.5),
    post = month >= treatment_month
  ) %>%
  arrange(nta_id, month)

# Verify balanced panel
n_units <- n_distinct(panel_final$nta_id)
n_periods <- n_distinct(panel_final$month)

cat(sprintf("Final panel: %d units × %d periods = %d rows\n",
            n_units, n_periods, nrow(panel_final)))

# Check balance
balance_check <- panel_final %>% count(nta_id) %>% pull(n) %>% unique()
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

ggsave(here("output", "causal_nta", "00_zone_time_series.png"), p_zone,
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
    title = "Zone vs Donor NTAs",
    subtitle = "NTAs provide better scale match than census tracts",
    x = "Month", y = "Crime Count", color = ""
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

p_comparison

ggsave(here("output", "causal_nta", "01_zone_vs_donors.png"), p_comparison,
       width = 10, height = 6, dpi = 300)

# Treatment pattern
p_treatment_pattern <- panel_final %>%
  mutate(
    status = case_when(
      tag == "treated" & post ~ "Treated (Post)",
      tag == "treated" & !post ~ "Treated (Pre)",
      TRUE ~ "Donor"
    )
  ) %>%
  count(month, status) %>%
  ggplot(aes(x = month, y = n, fill = status)) +
  geom_col() +
  geom_vline(xintercept = treatment_month, linetype = "dashed") +
  scale_fill_manual(values = c("Donor" = "gray70", 
                               "Treated (Pre)" = "steelblue",
                               "Treated (Post)" = "darkred")) +
  labs(
    title = "Treatment Pattern (NTA Level)",
    subtitle = "1 treated unit (Zone), rest are donor NTAs",
    x = "Month", y = "Number of Units", fill = ""
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
p_treatment_pattern

ggsave(here("output", "causal_nta", "02_treatment_pattern.png"), 
       p_treatment_pattern, width = 10, height = 6, dpi = 300)

#==============================================================================
# 7. PRIMARY METHOD: GSYNTH
#==============================================================================

cat("\n=== GSYNTH (PRIMARY METHOD) ===\n")
cat("Using generalized synthetic control with interactive fixed effects.\n")
cat("NTAs provide better scale match than census tracts.\n\n")

results <- list()

# Helper function to safely extract SE from gsynth
get_gsynth_se <- function(gs_obj) {
  if (is.null(gs_obj)) return(NA_real_)
  if (!is.null(gs_obj$est.avg) && is.matrix(gs_obj$est.avg)) {
    if ("S.E." %in% colnames(gs_obj$est.avg)) {
      return(gs_obj$est.avg[1, "S.E."])
    }
  }
  if (!is.null(gs_obj$att.avg.se) && is.numeric(gs_obj$att.avg.se) && length(gs_obj$att.avg.se) == 1) {
    return(gs_obj$att.avg.se)
  }
  if (!is.null(gs_obj$att.avg.var) && is.numeric(gs_obj$att.avg.var) && length(gs_obj$att.avg.var) == 1) {
    return(sqrt(gs_obj$att.avg.var))
  }
  return(NA_real_)
}

# --- Raw counts ---
cat("Fitting gsynth (raw counts)...\n")

tryCatch({
  gs_raw <- gsynth(
    crime_count ~ treat,
    data      = panel_final,
    index     = c("nta_id", "month"),
    force     = "two-way",
    CV        = TRUE,
    r         = c(0, 5),
    se        = TRUE,
    inference = "parametric",
    nboots    = 500,
    parallel  = TRUE,
    seed      = 42
  )
  
  results$gs_raw <- gs_raw
  
  cat(sprintf("\nGsynth Results (raw):\n"))
  cat(sprintf("  Factors selected: %d\n", gs_raw$r.cv))
  cat(sprintf("  Average ATT: %.2f\n", gs_raw$att.avg))
  
  se_raw <- get_gsynth_se(gs_raw)
  if (!is.na(se_raw)) {
    cat(sprintf("  SE: %.2f\n", se_raw))
    cat(sprintf("  95%% CI: [%.2f, %.2f]\n", gs_raw$att.avg - 1.96*se_raw, gs_raw$att.avg + 1.96*se_raw))
  }
  results$gs_raw_se <- se_raw
  
  # Plots
  png(here("output", "causal_nta", "03_gsynth_raw_counterfactual.png"),
      width = 10, height = 6, units = "in", res = 300)
  plot(gs_raw, type = "counterfactual", raw = "none",
       main = "Gsynth: Treated vs Synthetic Control (Raw Counts)")
  dev.off()
  
  png(here("output", "causal_nta", "04_gsynth_raw_gap.png"),
      width = 10, height = 6, units = "in", res = 300)
  plot(gs_raw, type = "gap", main = "Gsynth: Treatment Effect Over Time")
  dev.off()
  
  png(here("output", "causal_nta", "05_gsynth_raw_att.png"),
      width = 10, height = 6, units = "in", res = 300)
  plot(gs_raw, type = "ct", main = "Gsynth: Average Treatment Effect")
  dev.off()
  
  cat("  Plots saved.\n")
  
}, error = function(e) {
  cat(sprintf("Gsynth (raw) failed: %s\n", e$message))
})

# --- Log counts ---
cat("\nFitting gsynth (log counts)...\n")

tryCatch({
  gs_log <- gsynth(
    log_crime ~ treat,
    data      = panel_final,
    index     = c("nta_id", "month"),
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
  cat(sprintf("  Average ATT: %.3f\n", gs_log$att.avg))
  
  se_log <- get_gsynth_se(gs_log)
  if (!is.na(se_log)) {
    cat(sprintf("  SE: %.3f\n", se_log))
  }
  cat(sprintf("  Approx %% change: %.1f%%\n", (exp(gs_log$att.avg) - 1) * 100))
  results$gs_log_se <- se_log
  
  png(here("output", "causal_nta", "06_gsynth_log_counterfactual.png"),
      width = 10, height = 6, units = "in", res = 300)
  plot(gs_log, type = "counterfactual", raw = "none",
       main = "Gsynth: Treated vs Synthetic Control (Log Counts)")
  dev.off()
  
  png(here("output", "causal_nta", "07_gsynth_log_gap.png"),
      width = 10, height = 6, units = "in", res = 300)
  plot(gs_log, type = "gap", main = "Gsynth: Treatment Effect Over Time (Log)")
  dev.off()
  
  cat("  Plots saved.\n")
  
}, error = function(e) {
  cat(sprintf("Gsynth (log) failed: %s\n", e$message))
})

#==============================================================================
# 8. ROBUSTNESS: SYNTHETIC DIFF-IN-DIFF
#==============================================================================

cat("\n=== SYNTHETIC DIFF-IN-DIFF (ROBUSTNESS) ===\n")

# With NTAs we have fewer units, so use all of them
n_donors_sdid <- length(good_donors)

# Select by pre-treatment correlation with zone
zone_pre <- zone_monthly %>%
  filter(month < treatment_month) %>%
  arrange(month) %>%
  pull(crime_count)

donor_cors <- panel_final %>%
  filter(tag == "donor", month < treatment_month) %>%
  group_by(nta_id) %>%
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

sdid_donors <- donor_cors$nta_id

cat(sprintf("Using %d donors for SDID\n", length(sdid_donors)))

# Build matrix
panel_sdid <- panel_final %>%
  filter(nta_id == "ZONE" | nta_id %in% sdid_donors)

panel_wide_raw <- panel_sdid %>%
  select(nta_id, month, crime_count) %>%
  pivot_wider(names_from = month, values_from = crime_count, values_fill = 0) %>%
  column_to_rownames("nta_id")

panel_wide_log <- panel_sdid %>%
  select(nta_id, month, log_crime) %>%
  pivot_wider(names_from = month, values_from = log_crime, values_fill = log(0.5)) %>%
  column_to_rownames("nta_id")

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
  
  png(here("output", "causal_nta", "08_sdid_raw.png"),
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
  
  png(here("output", "causal_nta", "09_sdid_log.png"),
      width = 10, height = 6, units = "in", res = 300)
  synthdid_plot(sdid_log, se.method = "placebo")
  dev.off()
  
}, error = function(e) {
  cat(sprintf("SDID (log) failed: %s\n", e$message))
})

# --- Standard SC and DiD ---
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

augsynth_data <- panel_sdid %>%
  mutate(
    unit = nta_id,
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
    labs(title = "Augmented SCM (Raw Counts) - NTA") +
    theme_minimal()
  
  ggsave(here("output", "causal_nta", "10_ascm_raw.png"), p_ascm_raw,
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
    labs(title = "Augmented SCM (Log Counts) - NTA") +
    theme_minimal()
  
  ggsave(here("output", "causal_nta", "11_ascm_log.png"), p_ascm_log,
         width = 10, height = 6, dpi = 300)
  
}, error = function(e) {
  cat(sprintf("ASCM (log) failed: %s\n", e$message))
})

#==============================================================================
# 10. EXTRACT AND COMPARE ESTIMATES
#==============================================================================

cat("\n=== RESULTS SUMMARY ===\n")

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
  est <- results$gs_raw$att.avg
  se <- get_gsynth_se(results$gs_raw)
  results_table <- results_table %>%
    add_row(Method = "Gsynth", Scale = "Raw",
            Estimate = est, SE = se,
            CI_lower = if (!is.na(se)) est - 1.96*se else NA_real_,
            CI_upper = if (!is.na(se)) est + 1.96*se else NA_real_,
            Pct_Change = NA_real_)
}

# Gsynth log
if (!is.null(results$gs_log)) {
  est <- results$gs_log$att.avg
  se <- get_gsynth_se(results$gs_log)
  results_table <- results_table %>%
    add_row(Method = "Gsynth", Scale = "Log",
            Estimate = est, SE = se,
            CI_lower = if (!is.na(se)) est - 1.96*se else NA_real_,
            CI_upper = if (!is.na(se)) est + 1.96*se else NA_real_,
            Pct_Change = (exp(est) - 1) * 100)
}

# SDID raw
if (!is.null(results$sdid_raw)) {
  est <- as.numeric(results$sdid_raw)
  se <- if (!is.null(results$sdid_raw_se)) as.numeric(results$sdid_raw_se) else NA_real_
  results_table <- results_table %>%
    add_row(Method = "Synthetic DiD", Scale = "Raw",
            Estimate = est, SE = se,
            CI_lower = if (!is.na(se)) est - 1.96*se else NA_real_,
            CI_upper = if (!is.na(se)) est + 1.96*se else NA_real_,
            Pct_Change = NA_real_)
}

# SDID log
if (!is.null(results$sdid_log)) {
  est <- as.numeric(results$sdid_log)
  se <- if (!is.null(results$sdid_log_se)) as.numeric(results$sdid_log_se) else NA_real_
  results_table <- results_table %>%
    add_row(Method = "Synthetic DiD", Scale = "Log",
            Estimate = est, SE = se,
            CI_lower = if (!is.na(se)) est - 1.96*se else NA_real_,
            CI_upper = if (!is.na(se)) est + 1.96*se else NA_real_,
            Pct_Change = (exp(est) - 1) * 100)
}

# ASCM raw
if (!is.null(results$ascm_raw)) {
  summ <- summary(results$ascm_raw)
  est <- summ$average_att$Estimate
  se <- summ$average_att$Std.Error
  se <- if (!is.null(se) && is.numeric(se)) as.numeric(se) else NA_real_
  results_table <- results_table %>%
    add_row(Method = "Augmented SCM", Scale = "Raw",
            Estimate = est, SE = se,
            CI_lower = if (!is.na(se)) est - 1.96*se else NA_real_,
            CI_upper = if (!is.na(se)) est + 1.96*se else NA_real_,
            Pct_Change = NA_real_)
}

# ASCM log
if (!is.null(results$ascm_log)) {
  summ <- summary(results$ascm_log)
  est <- summ$average_att$Estimate
  se <- summ$average_att$Std.Error
  se <- if (!is.null(se) && is.numeric(se)) as.numeric(se) else NA_real_
  results_table <- results_table %>%
    add_row(Method = "Augmented SCM", Scale = "Log",
            Estimate = est, SE = se,
            CI_lower = if (!is.na(se)) est - 1.96*se else NA_real_,
            CI_upper = if (!is.na(se)) est + 1.96*se else NA_real_,
            Pct_Change = (exp(est) - 1) * 100)
}

# SC and DiD
if (!is.null(results$sc_raw)) {
  results_table <- results_table %>%
    add_row(Method = "Synthetic Control*", Scale = "Raw",
            Estimate = as.numeric(results$sc_raw), SE = NA_real_,
            CI_lower = NA_real_, CI_upper = NA_real_, Pct_Change = NA_real_)
}

if (!is.null(results$did_raw)) {
  results_table <- results_table %>%
    add_row(Method = "Diff-in-Diff", Scale = "Raw",
            Estimate = as.numeric(results$did_raw), SE = NA_real_,
            CI_lower = NA_real_, CI_upper = NA_real_, Pct_Change = NA_real_)
}

cat("\n")
print(results_table, n = 20)

cat("\n* Synthetic Control may be unreliable if level mismatch persists.\n")

write_csv(results_table, here("output", "causal_nta", "results_summary.csv"))

#==============================================================================
# 11. FOREST PLOT
#==============================================================================

cat("\n=== GENERATING PLOTS ===\n")

p_forest_raw <- results_table %>%
  filter(Scale == "Raw", !is.na(SE)) %>%
  ggplot(aes(x = Estimate, y = reorder(Method, Estimate))) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title = "Treatment Effect Estimates (Raw Scale) - NTA Level",
    subtitle = "Point estimates with 95% CI",
    x = "Estimated Effect (change in monthly crime count)", y = ""
  ) +
  theme_minimal()

p_forest_raw

ggsave(here("output", "causal_nta", "12_forest_plot_raw.png"), p_forest_raw,
       width = 10, height = 5, dpi = 300)

p_forest_log <- results_table %>%
  filter(Scale == "Log", !is.na(SE)) %>%
  ggplot(aes(x = Estimate, y = reorder(Method, Estimate))) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title = "Treatment Effect Estimates (Log Scale) - NTA Level",
    subtitle = "Point estimates with 95% CI. Negative = crime reduction.",
    x = "Estimated Effect (log scale)", y = ""
  ) +
  theme_minimal()

p_forest_log

ggsave(here("output", "causal_nta", "13_forest_plot_log.png"), p_forest_log,
       width = 10, height = 5, dpi = 300)

#==============================================================================
# 12. PRE-TREATMENT FIT DIAGNOSTICS
#==============================================================================

cat("\n=== PRE-TREATMENT FIT DIAGNOSTICS ===\n")

if (!is.null(results$gs_raw)) {
  
  gs <- results$gs_raw
  att <- gs$att
  rel_time <- gs$time
  
  unique_months <- sort(unique(panel_final$month))
  treatment_index <- which(unique_months == treatment_month)
  
  att_df <- tibble(
    relative_period = rel_time,
    ATT = att,
    month = unique_months[treatment_index + rel_time]
  ) %>%
    mutate(period = if_else(relative_period < 0, "Pre-treatment", "Post-treatment"))
  
  if (!is.null(gs$est.att)) {
    att_df <- att_df %>%
      mutate(
        SE = gs$est.att[, "S.E."],
        CI_lower = gs$est.att[, "CI.lower"],
        CI_upper = gs$est.att[, "CI.upper"],
        p_value = gs$est.att[, "p.value"]
      )
  }
  
  pre_att <- att_df %>% filter(period == "Pre-treatment")
  post_att <- att_df %>% filter(period == "Post-treatment")
  
  cat(sprintf("\nPre-treatment periods: %d\n", nrow(pre_att)))
  cat(sprintf("Pre-treatment ATT mean: %.2f (should be ~0 if good fit)\n", mean(pre_att$ATT)))
  cat(sprintf("Pre-treatment ATT SD: %.2f\n", sd(pre_att$ATT)))
  cat(sprintf("Pre-treatment ATT range: [%.2f, %.2f]\n", min(pre_att$ATT), max(pre_att$ATT)))
  
  cat(sprintf("\nPost-treatment periods: %d\n", nrow(post_att)))
  cat(sprintf("Post-treatment ATT mean: %.2f\n", mean(post_att$ATT)))
  
  p_att_time <- ggplot(att_df, aes(x = month, y = ATT)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    {if ("CI_lower" %in% names(att_df)) 
      geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2, fill = "steelblue")} +
    geom_line(linewidth = 0.8) +
    geom_point(aes(color = period), size = 2) +
    geom_vline(xintercept = treatment_month, linetype = "dashed") +
    annotate("rect", xmin = op_start, xmax = op_end, ymin = -Inf, ymax = Inf,
             alpha = 0.1, fill = "red") +
    scale_color_manual(values = c("Pre-treatment" = "gray50", "Post-treatment" = "darkred")) +
    labs(
      title = "Gsynth: Treatment Effect by Period (NTA Level)",
      subtitle = "Pre-treatment ATT should be ~0 if model fits well",
      x = "Month", y = "ATT (Treated - Counterfactual)", color = ""
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  p_att_time
  
  ggsave(here("output", "causal_nta", "14_att_by_period.png"), p_att_time,
         width = 10, height = 6, dpi = 300)
  
  write_csv(att_df, here("output", "causal_nta", "gsynth_att_by_period.csv"))
  
  pre_rmse <- sqrt(mean(pre_att$ATT^2))
  pre_rmse
  cat(sprintf("\nPre-treatment RMSE: %.2f\n", pre_rmse))
}

#==============================================================================
# 13. SAVE RESULTS
#==============================================================================

saveRDS(results, here("output", "causal_nta", "all_results.rds"))
saveRDS(panel_final, here("output", "causal_nta", "panel_final.rds"))

plan(sequential)

cat("\n================================================================\n")
cat("NTA-LEVEL ANALYSIS COMPLETE\n")
cat("================================================================\n")
cat(sprintf("\nOutputs saved to: %s\n", here("output", "causal_nta")))

list.files(here("output", "causal_nta")) %>%
  paste(" -", .) %>%
  cat(sep = "\n")

cat("\n")

