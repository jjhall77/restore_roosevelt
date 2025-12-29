#04-robustness
#==============================================================================
# Operation Restore Roosevelt - ROBUSTNESS CHECKS
# 
# Run after main analysis (03b-causal_final.R)
# Assumes all objects from main analysis are still in environment
#==============================================================================

library(here)
library(tidyverse)
library(lubridate)
library(sf)
library(gsynth)
library(synthdid)
library(future)
library(furrr)
library(patchwork)

# Output directory
dir.create(here("output", "causal_robustness"), showWarnings = FALSE, recursive = TRUE)

cat("\n")
cat("================================================================\n")
cat("ROBUSTNESS CHECKS FOR OPERATION RESTORE ROOSEVELT\n")
cat("================================================================\n")

# Store all robustness results
robustness_results <- list()

#==============================================================================
# 1. PRE-SPIKE MATCHING
#==============================================================================
# Re-estimate matching only on Jan 2022 - May 2024 (before summer spike)

cat("\n")
cat("================================================================\n")
cat("1. PRE-SPIKE MATCHING\n")
cat("================================================================\n")
cat("Matching on Jan 2022 - May 2024 only (avoids summer 2024 spike)\n\n")

pre_spike_end <- as.Date("2024-06-01")

# --- SDID with pre-spike weights ---
# Create matrix with only pre-spike periods for weight estimation
pre_spike_months <- all_months[all_months < pre_spike_end]
n_pre_spike <- length(pre_spike_months)

cat(sprintf("Pre-spike periods: %d (Jan 2022 - May 2024)\n", n_pre_spike))

# Build pre-spike panel for weight estimation
pre_spike_panel <- donor_panel %>%
  
  filter(nta_id %in% good_donors, month < pre_spike_end) %>%
  bind_rows(zone_panel %>% filter(month < pre_spike_end))

pre_spike_wide <- pre_spike_panel %>%
  select(nta_id, month, crime_count) %>%
  pivot_wider(names_from = month, values_from = crime_count, values_fill = 0) %>%
  column_to_rownames("nta_id")

mat_pre <- as.matrix(pre_spike_wide)
zone_row_pre <- mat_pre["ZONE", , drop = FALSE]
control_rows_pre <- mat_pre[rownames(mat_pre) != "ZONE", , drop = FALSE]
Y_pre_spike <- rbind(control_rows_pre, zone_row_pre)

N0_pre <- nrow(Y_pre_spike) - 1
T0_pre <- ncol(Y_pre_spike)  # All periods are pre-treatment in this subset

# Estimate weights on pre-spike data
cat("Estimating SDID weights on pre-spike data...\n")

tryCatch({
  # Get weights from pre-spike period
  sdid_pre_weights <- synthdid_estimate(Y_pre_spike, N0_pre, T0_pre - 1)
  omega_pre <- attr(sdid_pre_weights, "weights")$omega
  
  # Now apply these weights to full data
  # Manually compute ATT using pre-spike weights
  full_wide <- main_panel %>%
    select(nta_id, month, crime_count) %>%
    pivot_wider(names_from = month, values_from = crime_count, values_fill = 0) %>%
    column_to_rownames("nta_id")
  
  Y_full <- as.matrix(full_wide)
  zone_full <- Y_full["ZONE", ]
  controls_full <- Y_full[rownames(Y_full) != "ZONE", ]
  
  # Synthetic control using pre-spike weights
  synth_pre <- as.vector(omega_pre %*% controls_full)
  
  # ATT
  post_periods <- which(as.Date(colnames(Y_full)) >= treatment_month)
  att_pre_spike <- mean(zone_full[post_periods] - synth_pre[post_periods])
  
  robustness_results$pre_spike_sdid <- list(
    estimate = att_pre_spike,
    method = "SDID (pre-spike weights)",
    success = TRUE
  )
  
  cat(sprintf("  SDID (pre-spike weights) ATT: %.2f\n", att_pre_spike))
  
}, error = function(e) {
  cat(sprintf("  Pre-spike SDID failed: %s\n", e$message))
  robustness_results$pre_spike_sdid <- list(success = FALSE, error = e$message)
})

# --- Gsynth with pre-spike matching ---
cat("\nFitting gsynth with pre-spike emphasis...\n")

# Weight pre-spike periods more heavily by duplicating them
# Alternative: subset to pre-spike + post-treatment only
pre_spike_post_panel <- main_panel %>%
  filter(month < pre_spike_end | month >= treatment_month)

tryCatch({
  gs_pre_spike <- gsynth(
    crime_count ~ treat,
    data = pre_spike_post_panel %>%
      mutate(treat = as.integer(nta_id == "ZONE" & month >= treatment_month)),
    index = c("nta_id", "month"),
    force = "two-way",
    CV = TRUE,
    r = c(0, 5),
    se = TRUE,
    inference = "parametric",
    nboots = 500,
    parallel = TRUE,
    seed = 42
  )
  
  se <- if (!is.null(gs_pre_spike$est.avg) && is.matrix(gs_pre_spike$est.avg)) {
    gs_pre_spike$est.avg[1, "S.E."]
  } else NA_real_
  
  robustness_results$pre_spike_gsynth <- list(
    estimate = gs_pre_spike$att.avg,
    se = se,
    ci_lower = gs_pre_spike$att.avg - 1.96 * se,
    ci_upper = gs_pre_spike$att.avg + 1.96 * se,
    p_value = if (!is.na(se)) 2 * pnorm(-abs(gs_pre_spike$att.avg / se)) else NA,
    factors = gs_pre_spike$r.cv,
    method = "Gsynth (pre-spike only)",
    success = TRUE
  )
  
  cat(sprintf("  Gsynth (pre-spike) ATT: %.2f (SE: %.2f)\n", 
              gs_pre_spike$att.avg, se))
  
}, error = function(e) {
  cat(sprintf("  Pre-spike gsynth failed: %s\n", e$message))
  robustness_results$pre_spike_gsynth <- list(success = FALSE, error = e$message)
})

#==============================================================================
# 2. PLACEBO TIMING TEST
#==============================================================================
# Assign fake treatment at Oct 2023 (1 year earlier)
# Should find NO effect if model is well-specified

cat("\n")
cat("================================================================\n")
cat("2. PLACEBO TIMING TEST\n")
cat("================================================================\n")
cat("Fake treatment at Oct 2023 (should find NO effect)\n\n")

placebo_treatment <- as.Date("2023-10-01")

# Only use data up to actual treatment (avoid contamination)
placebo_panel <- main_panel %>%
  filter(month < treatment_month) %>%
  mutate(
    treat_placebo = as.integer(nta_id == "ZONE" & month >= placebo_treatment)
  )

placebo_n_pre <- sum(all_months < placebo_treatment & all_months < treatment_month)
placebo_n_post <- sum(all_months >= placebo_treatment & all_months < treatment_month)

cat(sprintf("Placebo treatment: %s\n", placebo_treatment))
cat(sprintf("Placebo periods: %d pre, %d post (before actual treatment)\n", 
            placebo_n_pre, placebo_n_post))

# --- Gsynth placebo ---
cat("\nFitting gsynth placebo...\n")

tryCatch({
  gs_placebo <- gsynth(
    crime_count ~ treat_placebo,
    data = placebo_panel,
    index = c("nta_id", "month"),
    force = "two-way",
    CV = TRUE,
    r = c(0, 5),
    se = TRUE,
    inference = "parametric",
    nboots = 500,
    parallel = TRUE,
    seed = 42
  )
  
  se <- if (!is.null(gs_placebo$est.avg) && is.matrix(gs_placebo$est.avg)) {
    gs_placebo$est.avg[1, "S.E."]
  } else NA_real_
  
  p_val <- if (!is.na(se)) 2 * pnorm(-abs(gs_placebo$att.avg / se)) else NA
  
  robustness_results$placebo_gsynth <- list(
    estimate = gs_placebo$att.avg,
    se = se,
    ci_lower = gs_placebo$att.avg - 1.96 * se,
    ci_upper = gs_placebo$att.avg + 1.96 * se,
    p_value = p_val,
    method = "Gsynth (placebo Oct 2023)",
    success = TRUE
  )
  
  cat(sprintf("  Gsynth placebo ATT: %.2f (SE: %.2f, p=%.4f)\n", 
              gs_placebo$att.avg, se, p_val))
  
  if (p_val > 0.1) {
    cat("  ✓ PASS: No significant placebo effect (as expected)\n")
  } else {
    cat("  ⚠ WARNING: Significant placebo effect detected - model may be misspecified\n")
  }
  
  # Save placebo plot
  p_placebo <- plot(gs_placebo, type = "gap", main = "Placebo Test: Fake Treatment Oct 2023")
  print(p_placebo)
  ggsave(here("output", "causal_robustness", "placebo_gsynth_gap.png"), p_placebo,
         width = 10, height = 6, dpi = 300)
  
}, error = function(e) {
  cat(sprintf("  Gsynth placebo failed: %s\n", e$message))
  robustness_results$placebo_gsynth <- list(success = FALSE, error = e$message)
})

# --- SDID placebo ---
cat("\nFitting SDID placebo...\n")

tryCatch({
  placebo_wide <- placebo_panel %>%
    select(nta_id, month, crime_count) %>%
    pivot_wider(names_from = month, values_from = crime_count, values_fill = 0) %>%
    column_to_rownames("nta_id")
  
  mat_placebo <- as.matrix(placebo_wide)
  zone_row_placebo <- mat_placebo["ZONE", , drop = FALSE]
  control_rows_placebo <- mat_placebo[rownames(mat_placebo) != "ZONE", , drop = FALSE]
  Y_placebo <- rbind(control_rows_placebo, zone_row_placebo)
  
  N0_placebo <- nrow(Y_placebo) - 1
  T0_placebo <- sum(as.Date(colnames(Y_placebo)) < placebo_treatment)
  
  sdid_placebo <- synthdid_estimate(Y_placebo, N0_placebo, T0_placebo)
  se_placebo <- sqrt(vcov(sdid_placebo, method = "placebo"))
  p_val_placebo <- 2 * pnorm(-abs(as.numeric(sdid_placebo) / se_placebo))
  
  robustness_results$placebo_sdid <- list(
    estimate = as.numeric(sdid_placebo),
    se = se_placebo,
    ci_lower = as.numeric(sdid_placebo) - 1.96 * se_placebo,
    ci_upper = as.numeric(sdid_placebo) + 1.96 * se_placebo,
    p_value = p_val_placebo,
    method = "SDID (placebo Oct 2023)",
    success = TRUE
  )
  
  cat(sprintf("  SDID placebo ATT: %.2f (SE: %.2f, p=%.4f)\n", 
              as.numeric(sdid_placebo), se_placebo, p_val_placebo))
  
  if (p_val_placebo > 0.1) {
    cat("  ✓ PASS: No significant placebo effect\n")
  } else {
    cat("  ⚠ WARNING: Significant placebo effect detected\n")
  }
  
  # Save placebo plot
  p_sdid_placebo <- synthdid_plot(sdid_placebo, se.method = "placebo")
  print(p_sdid_placebo)
  ggsave(here("output", "causal_robustness", "placebo_sdid.png"), p_sdid_placebo,
         width = 10, height = 6, dpi = 300)
  
}, error = function(e) {
  cat(sprintf("  SDID placebo failed: %s\n", e$message))
  robustness_results$placebo_sdid <- list(success = FALSE, error = e$message)
})

#==============================================================================
# 3. LEAVE-ONE-OUT DONORS
#==============================================================================
# Remove top 5 most-weighted donors, re-estimate

cat("\n")
cat("================================================================\n")
cat("3. LEAVE-ONE-OUT DONORS\n")
cat("================================================================\n")
cat("Dropping top 5 weighted donors to check sensitivity\n\n")

# Get top 5 donors by SDID weight
if (exists("sdid_weights")) {
  top5_donors <- sdid_weights %>%
    filter(sdid_weight > 0) %>%
    slice_max(sdid_weight, n = 5) %>%
    pull(nta_id)
  
  cat("Top 5 donors being dropped:\n")
  print(sdid_weights %>% filter(nta_id %in% top5_donors))
  
  # Reduced donor pool
  reduced_donors <- setdiff(good_donors, top5_donors)
  cat(sprintf("\nDonors remaining: %d (from %d)\n", length(reduced_donors), length(good_donors)))
  
  # --- Gsynth LOO ---
  cat("\nFitting gsynth with reduced donors...\n")
  
  loo_panel <- donor_panel %>%
    filter(nta_id %in% reduced_donors) %>%
    bind_rows(zone_panel)
  
  tryCatch({
    gs_loo <- gsynth(
      crime_count ~ treat,
      data = loo_panel %>%
        mutate(treat = as.integer(nta_id == "ZONE" & month >= treatment_month)),
      index = c("nta_id", "month"),
      force = "two-way",
      CV = TRUE,
      r = c(0, 5),
      se = TRUE,
      inference = "parametric",
      nboots = 500,
      parallel = TRUE,
      seed = 42
    )
    
    se <- if (!is.null(gs_loo$est.avg) && is.matrix(gs_loo$est.avg)) {
      gs_loo$est.avg[1, "S.E."]
    } else NA_real_
    
    robustness_results$loo_gsynth <- list(
      estimate = gs_loo$att.avg,
      se = se,
      ci_lower = gs_loo$att.avg - 1.96 * se,
      ci_upper = gs_loo$att.avg + 1.96 * se,
      p_value = if (!is.na(se)) 2 * pnorm(-abs(gs_loo$att.avg / se)) else NA,
      method = "Gsynth (LOO top 5)",
      success = TRUE
    )
    
    cat(sprintf("  Gsynth LOO ATT: %.2f (SE: %.2f)\n", gs_loo$att.avg, se))
    cat(sprintf("  Original estimate: %.2f\n", gs_zone$estimate))
    cat(sprintf("  Difference: %.2f (%.1f%%)\n", 
                gs_loo$att.avg - gs_zone$estimate,
                (gs_loo$att.avg - gs_zone$estimate) / abs(gs_zone$estimate) * 100))
    
  }, error = function(e) {
    cat(sprintf("  Gsynth LOO failed: %s\n", e$message))
    robustness_results$loo_gsynth <- list(success = FALSE, error = e$message)
  })
  
  # --- SDID LOO ---
  cat("\nFitting SDID with reduced donors...\n")
  
  tryCatch({
    loo_matrix <- prepare_sdid_matrix(zone_panel, donor_panel, reduced_donors)
    sdid_loo <- synthdid_estimate(loo_matrix$Y, loo_matrix$N0, loo_matrix$T0)
    se_loo <- sqrt(vcov(sdid_loo, method = "placebo"))
    
    robustness_results$loo_sdid <- list(
      estimate = as.numeric(sdid_loo),
      se = se_loo,
      ci_lower = as.numeric(sdid_loo) - 1.96 * se_loo,
      ci_upper = as.numeric(sdid_loo) + 1.96 * se_loo,
      p_value = 2 * pnorm(-abs(as.numeric(sdid_loo) / se_loo)),
      method = "SDID (LOO top 5)",
      success = TRUE
    )
    
    cat(sprintf("  SDID LOO ATT: %.2f (SE: %.2f)\n", as.numeric(sdid_loo), se_loo))
    cat(sprintf("  Original estimate: %.2f\n", sdid_zone$estimate))
    cat(sprintf("  Difference: %.2f (%.1f%%)\n", 
                as.numeric(sdid_loo) - sdid_zone$estimate,
                (as.numeric(sdid_loo) - sdid_zone$estimate) / abs(sdid_zone$estimate) * 100))
    
  }, error = function(e) {
    cat(sprintf("  SDID LOO failed: %s\n", e$message))
    robustness_results$loo_sdid <- list(success = FALSE, error = e$message)
  })
  
} else {
  cat("SDID weights not available - skipping LOO analysis\n")
}

#==============================================================================
# 4. ALTERNATIVE DONOR POOLS
#==============================================================================

cat("\n")
cat("================================================================\n")
cat("4. ALTERNATIVE DONOR POOLS\n")
cat("================================================================\n")

# --- 4a. Queens-only donors ---
cat("\n--- 4a. Queens-Only Donors ---\n")

queens_donors <- nynta %>%
  filter(boro_name == "Queens") %>%
  pull(nta_id) %>%
  intersect(good_donors)

cat(sprintf("Queens donors: %d\n", length(queens_donors)))

if (length(queens_donors) >= 10) {
  tryCatch({
    queens_panel <- donor_panel %>%
      filter(nta_id %in% queens_donors) %>%
      bind_rows(zone_panel)
    
    gs_queens <- gsynth(
      crime_count ~ treat,
      data = queens_panel %>%
        mutate(treat = as.integer(nta_id == "ZONE" & month >= treatment_month)),
      index = c("nta_id", "month"),
      force = "two-way",
      CV = TRUE,
      r = c(0, 5),
      se = TRUE,
      inference = "parametric",
      nboots = 500,
      parallel = TRUE,
      seed = 42
    )
    
    se <- if (!is.null(gs_queens$est.avg) && is.matrix(gs_queens$est.avg)) {
      gs_queens$est.avg[1, "S.E."]
    } else NA_real_
    
    robustness_results$queens_gsynth <- list(
      estimate = gs_queens$att.avg,
      se = se,
      ci_lower = gs_queens$att.avg - 1.96 * se,
      ci_upper = gs_queens$att.avg + 1.96 * se,
      p_value = if (!is.na(se)) 2 * pnorm(-abs(gs_queens$att.avg / se)) else NA,
      method = "Gsynth (Queens only)",
      success = TRUE
    )
    
    cat(sprintf("  Gsynth (Queens) ATT: %.2f (SE: %.2f)\n", gs_queens$att.avg, se))
    
  }, error = function(e) {
    cat(sprintf("  Queens gsynth failed: %s\n", e$message))
    robustness_results$queens_gsynth <- list(success = FALSE, error = e$message)
  })
  
  # SDID Queens
  tryCatch({
    queens_matrix <- prepare_sdid_matrix(zone_panel, donor_panel, queens_donors)
    sdid_queens <- synthdid_estimate(queens_matrix$Y, queens_matrix$N0, queens_matrix$T0)
    se_queens <- sqrt(vcov(sdid_queens, method = "placebo"))
    
    robustness_results$queens_sdid <- list(
      estimate = as.numeric(sdid_queens),
      se = se_queens,
      ci_lower = as.numeric(sdid_queens) - 1.96 * se_queens,
      ci_upper = as.numeric(sdid_queens) + 1.96 * se_queens,
      p_value = 2 * pnorm(-abs(as.numeric(sdid_queens) / se_queens)),
      method = "SDID (Queens only)",
      success = TRUE
    )
    
    cat(sprintf("  SDID (Queens) ATT: %.2f (SE: %.2f)\n", as.numeric(sdid_queens), se_queens))
    
  }, error = function(e) {
    cat(sprintf("  Queens SDID failed: %s\n", e$message))
    robustness_results$queens_sdid <- list(success = FALSE, error = e$message)
  })
  
} else {
  cat("  Not enough Queens donors for analysis\n")
}

# --- 4b. Non-Queens donors ---
cat("\n--- 4b. Non-Queens Donors ---\n")

non_queens_donors <- nynta %>%
  filter(boro_name != "Queens") %>%
  pull(nta_id) %>%
  intersect(good_donors)

cat(sprintf("Non-Queens donors: %d\n", length(non_queens_donors)))

if (length(non_queens_donors) >= 10) {
  tryCatch({
    non_queens_panel <- donor_panel %>%
      filter(nta_id %in% non_queens_donors) %>%
      bind_rows(zone_panel)
    
    gs_non_queens <- gsynth(
      crime_count ~ treat,
      data = non_queens_panel %>%
        mutate(treat = as.integer(nta_id == "ZONE" & month >= treatment_month)),
      index = c("nta_id", "month"),
      force = "two-way",
      CV = TRUE,
      r = c(0, 5),
      se = TRUE,
      inference = "parametric",
      nboots = 500,
      parallel = TRUE,
      seed = 42
    )
    
    se <- if (!is.null(gs_non_queens$est.avg) && is.matrix(gs_non_queens$est.avg)) {
      gs_non_queens$est.avg[1, "S.E."]
    } else NA_real_
    
    robustness_results$non_queens_gsynth <- list(
      estimate = gs_non_queens$att.avg,
      se = se,
      ci_lower = gs_non_queens$att.avg - 1.96 * se,
      ci_upper = gs_non_queens$att.avg + 1.96 * se,
      p_value = if (!is.na(se)) 2 * pnorm(-abs(gs_non_queens$att.avg / se)) else NA,
      method = "Gsynth (Non-Queens)",
      success = TRUE
    )
    
    cat(sprintf("  Gsynth (Non-Queens) ATT: %.2f (SE: %.2f)\n", gs_non_queens$att.avg, se))
    
  }, error = function(e) {
    cat(sprintf("  Non-Queens gsynth failed: %s\n", e$message))
    robustness_results$non_queens_gsynth <- list(success = FALSE, error = e$message)
  })
}

# --- 4c. High-crime NTAs only ---
cat("\n--- 4c. High-Crime NTAs Only ---\n")

# NTAs with mean crime >= zone's 25th percentile
zone_p25 <- quantile(zone_panel$crime_count, 0.25)

high_crime_donors <- donor_stats %>%
  filter(nta_id %in% good_donors) %>%
  filter(mean_crime >= zone_p25 * 0.5) %>%  # At least half of zone's 25th %ile
  pull(nta_id)

cat(sprintf("High-crime donors (mean >= %.1f): %d\n", zone_p25 * 0.5, length(high_crime_donors)))

if (length(high_crime_donors) >= 10) {
  tryCatch({
    high_crime_panel <- donor_panel %>%
      filter(nta_id %in% high_crime_donors) %>%
      bind_rows(zone_panel)
    
    gs_high_crime <- gsynth(
      crime_count ~ treat,
      data = high_crime_panel %>%
        mutate(treat = as.integer(nta_id == "ZONE" & month >= treatment_month)),
      index = c("nta_id", "month"),
      force = "two-way",
      CV = TRUE,
      r = c(0, 5),
      se = TRUE,
      inference = "parametric",
      nboots = 500,
      parallel = TRUE,
      seed = 42
    )
    
    se <- if (!is.null(gs_high_crime$est.avg) && is.matrix(gs_high_crime$est.avg)) {
      gs_high_crime$est.avg[1, "S.E."]
    } else NA_real_
    
    robustness_results$high_crime_gsynth <- list(
      estimate = gs_high_crime$att.avg,
      se = se,
      ci_lower = gs_high_crime$att.avg - 1.96 * se,
      ci_upper = gs_high_crime$att.avg + 1.96 * se,
      p_value = if (!is.na(se)) 2 * pnorm(-abs(gs_high_crime$att.avg / se)) else NA,
      method = "Gsynth (High-crime NTAs)",
      success = TRUE
    )
    
    cat(sprintf("  Gsynth (High-crime) ATT: %.2f (SE: %.2f)\n", gs_high_crime$att.avg, se))
    
  }, error = function(e) {
    cat(sprintf("  High-crime gsynth failed: %s\n", e$message))
    robustness_results$high_crime_gsynth <- list(success = FALSE, error = e$message)
  })
}

#==============================================================================
# 5. DIFFERENT OUTCOME MEASURES
#==============================================================================

cat("\n")
cat("================================================================\n")
cat("5. DIFFERENT OUTCOME MEASURES\n")
cat("================================================================\n")

# --- 5a. Log transformation ---
cat("\n--- 5a. Log Transformation ---\n")

log_panel <- main_panel %>%
  mutate(
    log_crime = log(crime_count + 0.5),
    treat = as.integer(nta_id == "ZONE" & month >= treatment_month)
  )

tryCatch({
  gs_log <- gsynth(
    log_crime ~ treat,
    data = log_panel,
    index = c("nta_id", "month"),
    force = "two-way",
    CV = TRUE,
    r = c(0, 5),
    se = TRUE,
    inference = "parametric",
    nboots = 500,
    parallel = TRUE,
    seed = 42
  )
  
  se <- if (!is.null(gs_log$est.avg) && is.matrix(gs_log$est.avg)) {
    gs_log$est.avg[1, "S.E."]
  } else NA_real_
  
  pct_change <- (exp(gs_log$att.avg) - 1) * 100
  
  robustness_results$log_gsynth <- list(
    estimate = gs_log$att.avg,
    se = se,
    ci_lower = gs_log$att.avg - 1.96 * se,
    ci_upper = gs_log$att.avg + 1.96 * se,
    p_value = if (!is.na(se)) 2 * pnorm(-abs(gs_log$att.avg / se)) else NA,
    pct_change = pct_change,
    method = "Gsynth (Log)",
    success = TRUE
  )
  
  cat(sprintf("  Gsynth (Log) ATT: %.3f (SE: %.3f)\n", gs_log$att.avg, se))
  cat(sprintf("  Implied %% change: %.1f%%\n", pct_change))
  
}, error = function(e) {
  cat(sprintf("  Log gsynth failed: %s\n", e$message))
  robustness_results$log_gsynth <- list(success = FALSE, error = e$message)
})

# --- 5b. Specific crime types ---
cat("\n--- 5b. Crime Type: Robbery Only ---\n")

# Build robbery-only panel
robbery_zone <- street_violent_crime %>%
  filter(date >= data_start) %>%
  filter(str_detect(tolower(ofns_desc), "robbery")) %>%
  st_filter(rosie_buffer, .predicate = st_intersects) %>%
  st_drop_geometry() %>%
  mutate(month = floor_date(date, "month")) %>%
  filter(month <= data_end) %>%
  count(month, name = "crime_count") %>%
  right_join(tibble(month = all_months), by = "month") %>%
  replace_na(list(crime_count = 0)) %>%
  mutate(nta_id = "ZONE") %>%
  arrange(month)

robbery_by_nta <- crime_with_nta %>%
  st_drop_geometry() %>%
  filter(str_detect(tolower(ofns_desc), "robbery")) %>%
  filter(!cmplnt_num %in% excluded_crimes) %>%
  count(nta_id, month, name = "crime_count")

robbery_donor <- expand_grid(nta_id = good_donors, month = all_months) %>%
  left_join(robbery_by_nta, by = c("nta_id", "month")) %>%
  replace_na(list(crime_count = 0))

robbery_panel <- robbery_donor %>%
  bind_rows(robbery_zone) %>%
  mutate(treat = as.integer(nta_id == "ZONE" & month >= treatment_month))

cat(sprintf("Zone robbery: mean=%.1f/month\n", mean(robbery_zone$crime_count)))

tryCatch({
  gs_robbery <- gsynth(
    crime_count ~ treat,
    data = robbery_panel,
    index = c("nta_id", "month"),
    force = "two-way",
    CV = TRUE,
    r = c(0, 5),
    se = TRUE,
    inference = "parametric",
    nboots = 500,
    parallel = TRUE,
    seed = 42
  )
  
  se <- if (!is.null(gs_robbery$est.avg) && is.matrix(gs_robbery$est.avg)) {
    gs_robbery$est.avg[1, "S.E."]
  } else NA_real_
  
  robustness_results$robbery_gsynth <- list(
    estimate = gs_robbery$att.avg,
    se = se,
    ci_lower = gs_robbery$att.avg - 1.96 * se,
    ci_upper = gs_robbery$att.avg + 1.96 * se,
    p_value = if (!is.na(se)) 2 * pnorm(-abs(gs_robbery$att.avg / se)) else NA,
    method = "Gsynth (Robbery)",
    success = TRUE
  )
  
  cat(sprintf("  Gsynth (Robbery) ATT: %.2f (SE: %.2f)\n", gs_robbery$att.avg, se))
  
}, error = function(e) {
  cat(sprintf("  Robbery gsynth failed: %s\n", e$message))
  robustness_results$robbery_gsynth <- list(success = FALSE, error = e$message)
})

# --- 5c. Assault only ---
cat("\n--- 5c. Crime Type: Assault Only ---\n")

assault_zone <- street_violent_crime %>%
  filter(date >= data_start) %>%
  filter(str_detect(tolower(ofns_desc), "assault")) %>%
  st_filter(rosie_buffer, .predicate = st_intersects) %>%
  st_drop_geometry() %>%
  mutate(month = floor_date(date, "month")) %>%
  filter(month <= data_end) %>%
  count(month, name = "crime_count") %>%
  right_join(tibble(month = all_months), by = "month") %>%
  replace_na(list(crime_count = 0)) %>%
  mutate(nta_id = "ZONE") %>%
  arrange(month)

assault_by_nta <- crime_with_nta %>%
  st_drop_geometry() %>%
  filter(str_detect(tolower(ofns_desc), "assault")) %>%
  filter(!cmplnt_num %in% excluded_crimes) %>%
  count(nta_id, month, name = "crime_count")

assault_donor <- expand_grid(nta_id = good_donors, month = all_months) %>%
  left_join(assault_by_nta, by = c("nta_id", "month")) %>%
  replace_na(list(crime_count = 0))

assault_panel <- assault_donor %>%
  bind_rows(assault_zone) %>%
  mutate(treat = as.integer(nta_id == "ZONE" & month >= treatment_month))

cat(sprintf("Zone assault: mean=%.1f/month\n", mean(assault_zone$crime_count)))

tryCatch({
  gs_assault <- gsynth(
    crime_count ~ treat,
    data = assault_panel,
    index = c("nta_id", "month"),
    force = "two-way",
    CV = TRUE,
    r = c(0, 5),
    se = TRUE,
    inference = "parametric",
    nboots = 500,
    parallel = TRUE,
    seed = 42
  )
  
  se <- if (!is.null(gs_assault$est.avg) && is.matrix(gs_assault$est.avg)) {
    gs_assault$est.avg[1, "S.E."]
  } else NA_real_
  
  robustness_results$assault_gsynth <- list(
    estimate = gs_assault$att.avg,
    se = se,
    ci_lower = gs_assault$att.avg - 1.96 * se,
    ci_upper = gs_assault$att.avg + 1.96 * se,
    p_value = if (!is.na(se)) 2 * pnorm(-abs(gs_assault$att.avg / se)) else NA,
    method = "Gsynth (Assault)",
    success = TRUE
  )
  
  cat(sprintf("  Gsynth (Assault) ATT: %.2f (SE: %.2f)\n", gs_assault$att.avg, se))
  
}, error = function(e) {
  cat(sprintf("  Assault gsynth failed: %s\n", e$message))
  robustness_results$assault_gsynth <- list(success = FALSE, error = e$message)
})

#==============================================================================
# 6. INTERRUPTED TIME SERIES (NO DONORS)
#==============================================================================

cat("\n")
cat("================================================================\n")
cat("6. INTERRUPTED TIME SERIES\n")
cat("================================================================\n")
cat("Zone-only analysis with segmented regression\n\n")

# Prepare ITS data
its_data <- zone_panel %>%
  mutate(
    time = row_number(),
    post = as.integer(month >= treatment_month),
    time_since_treat = pmax(0, time - n_pre)
  )

# --- Simple ITS model ---
cat("--- 6a. Simple Level Change ---\n")

its_simple <- lm(crime_count ~ time + post, data = its_data)
summary_its <- summary(its_simple)

cat("\nITS Model: crime ~ time + post\n")
cat(sprintf("  Pre-trend: %.2f crimes/month\n", coef(its_simple)["time"]))
cat(sprintf("  Level change at treatment: %.2f (SE: %.2f)\n", 
            coef(its_simple)["post"], 
            summary_its$coefficients["post", "Std. Error"]))
cat(sprintf("  p-value: %.4f\n", summary_its$coefficients["post", "Pr(>|t|)"]))

robustness_results$its_simple <- list(
  estimate = coef(its_simple)["post"],
  se = summary_its$coefficients["post", "Std. Error"],
  p_value = summary_its$coefficients["post", "Pr(>|t|)"],
  method = "ITS (level change)",
  success = TRUE
)

# --- ITS with slope change ---
cat("\n--- 6b. Level + Slope Change ---\n")

its_slope <- lm(crime_count ~ time + post + time_since_treat, data = its_data)
summary_slope <- summary(its_slope)

cat("\nITS Model: crime ~ time + post + time_since_treat\n")
cat(sprintf("  Pre-trend: %.2f crimes/month\n", coef(its_slope)["time"]))
cat(sprintf("  Level change: %.2f (SE: %.2f, p=%.4f)\n", 
            coef(its_slope)["post"],
            summary_slope$coefficients["post", "Std. Error"],
            summary_slope$coefficients["post", "Pr(>|t|)"]))
cat(sprintf("  Slope change: %.2f (SE: %.2f, p=%.4f)\n", 
            coef(its_slope)["time_since_treat"],
            summary_slope$coefficients["time_since_treat", "Std. Error"],
            summary_slope$coefficients["time_since_treat", "Pr(>|t|)"]))

robustness_results$its_slope <- list(
  level_change = coef(its_slope)["post"],
  slope_change = coef(its_slope)["time_since_treat"],
  level_se = summary_slope$coefficients["post", "Std. Error"],
  slope_se = summary_slope$coefficients["time_since_treat", "Std. Error"],
  level_p = summary_slope$coefficients["post", "Pr(>|t|)"],
  slope_p = summary_slope$coefficients["time_since_treat", "Pr(>|t|)"],
  method = "ITS (level + slope)",
  success = TRUE
)

# --- ITS Plot ---
its_data <- its_data %>%
  mutate(
    fitted_simple = predict(its_simple),
    fitted_slope = predict(its_slope),
    counterfactual = coef(its_slope)["(Intercept)"] + coef(its_slope)["time"] * time
  )

p_its <- ggplot(its_data, aes(x = month)) +
  geom_point(aes(y = crime_count), color = "gray40") +
  geom_line(aes(y = crime_count), color = "gray40", alpha = 0.5) +
  geom_line(aes(y = fitted_slope, color = "Fitted (level + slope)"), linewidth = 1) +
  geom_line(aes(y = counterfactual, color = "Counterfactual"), 
            linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = treatment_month, linetype = "dashed", color = "black") +
  annotate("rect", xmin = op_start, xmax = op_end, ymin = -Inf, ymax = Inf,
           alpha = 0.15, fill = "red") +
  scale_color_manual(values = c("Fitted (level + slope)" = "steelblue", 
                                "Counterfactual" = "darkred")) +
  labs(
    title = "Interrupted Time Series: Zone Crime",
    subtitle = sprintf("Level change: %.1f, Slope change: %.2f", 
                       coef(its_slope)["post"], coef(its_slope)["time_since_treat"]),
    x = "Month", y = "Crime Count", color = ""
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p_its)
ggsave(here("output", "causal_robustness", "its_analysis.png"), p_its,
       width = 10, height = 6, dpi = 300)

#==============================================================================
# 7. SUMMARY TABLE
#==============================================================================

cat("\n")
cat("================================================================\n")
cat("ROBUSTNESS SUMMARY\n")
cat("================================================================\n")

# Build summary table
summary_rows <- list()

# Main estimates (from original analysis)
summary_rows[[1]] <- tibble(
  Check = "Main Analysis",
  Method = "Gsynth",
  Estimate = gs_zone$estimate,
  SE = gs_zone$se,
  P_Value = gs_zone$p_value
)

summary_rows[[2]] <- tibble(
  Check = "Main Analysis",
  Method = "SDID",
  Estimate = sdid_zone$estimate,
  SE = sdid_zone$se,
  P_Value = sdid_zone$p_value
)

# Add robustness results
for (name in names(robustness_results)) {
  r <- robustness_results[[name]]
  if (r$success && !is.null(r$estimate)) {
    summary_rows[[length(summary_rows) + 1]] <- tibble(
      Check = r$method,
      Method = name,
      Estimate = r$estimate,
      SE = if (!is.null(r$se)) r$se else NA_real_,
      P_Value = if (!is.null(r$p_value)) r$p_value else NA_real_
    )
  }
}

robustness_table <- bind_rows(summary_rows) %>%
  mutate(
    Significant = case_when(
      P_Value < 0.01 ~ "***",
      P_Value < 0.05 ~ "**",
      P_Value < 0.1 ~ "*",
      TRUE ~ ""
    )
  )

print(robustness_table, n = 30)

write_csv(robustness_table, here("output", "causal_robustness", "robustness_summary.csv"))

#==============================================================================
# 8. FOREST PLOT OF ALL ESTIMATES
#==============================================================================

cat("\n=== GENERATING ROBUSTNESS FOREST PLOT ===\n")

# Filter to comparable estimates (raw scale, not log)
forest_data <- robustness_table %>%
  filter(!is.na(SE), !str_detect(Check, "Log|ITS")) %>%
  mutate(
    CI_Lower = Estimate - 1.96 * SE,
    CI_Upper = Estimate + 1.96 * SE,
    Category = case_when(
      str_detect(Check, "Main") ~ "Main",
      str_detect(Check, "placebo") ~ "Placebo",
      str_detect(Check, "LOO|Queens|Non-Queens|High-crime") ~ "Donor Pool",
      str_detect(Check, "pre-spike") ~ "Pre-spike",
      TRUE ~ "Other"
    )
  )

p_robustness_forest <- ggplot(forest_data, aes(x = Estimate, y = reorder(Check, Estimate))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = gs_zone$estimate, linetype = "dotted", color = "darkred", alpha = 0.5) +
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0.2) +
  geom_point(aes(color = Category), size = 3) +
  scale_color_manual(values = c(
    "Main" = "darkred",
    "Placebo" = "purple",
    "Donor Pool" = "steelblue",
    "Pre-spike" = "forestgreen",
    "Other" = "gray50"
  )) +
  labs(
    title = "Robustness Check: Treatment Effect Estimates",
    subtitle = "Dotted line = main gsynth estimate. All should be negative if effect is robust.",
    x = "Estimated Effect (change in monthly crimes)", y = ""
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p_robustness_forest)
ggsave(here("output", "causal_robustness", "robustness_forest.png"), p_robustness_forest,
       width = 12, height = 8, dpi = 300)

#==============================================================================
# 9. SAVE ALL RESULTS
#==============================================================================

saveRDS(robustness_results, here("output", "causal_robustness", "all_robustness_results.rds"))

cat("\n================================================================\n")
cat("ROBUSTNESS CHECKS COMPLETE\n")
cat("================================================================\n")

cat("\nKey findings:\n")

# Placebo check
if (!is.null(robustness_results$placebo_gsynth$p_value)) {
  if (robustness_results$placebo_gsynth$p_value > 0.1) {
    cat("✓ Placebo test PASSED (no false positive)\n")
  } else {
    cat("⚠ Placebo test FAILED (spurious effect detected)\n")
  }
}

# Sensitivity to donors
if (!is.null(robustness_results$loo_gsynth$estimate)) {
  pct_diff <- abs(robustness_results$loo_gsynth$estimate - gs_zone$estimate) / abs(gs_zone$estimate) * 100
  if (pct_diff < 25) {
    cat(sprintf("✓ Leave-one-out: Robust (%.1f%% change)\n", pct_diff))
  } else {
    cat(sprintf("⚠ Leave-one-out: Sensitive to top donors (%.1f%% change)\n", pct_diff))
  }
}

# Cross-method consistency
main_estimates <- c(gs_zone$estimate, sdid_zone$estimate)
if (all(main_estimates < 0)) {
  cat("✓ All main methods agree on direction (negative effect)\n")
}

cat("\nOutputs saved to:", here("output", "causal_robustness"), "\n")

list.files(here("output", "causal_robustness")) %>%
  paste(" -", .) %>%
  cat(sep = "\n")