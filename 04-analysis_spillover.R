#==============================================================================
# Operation Restore Roosevelt - Causal Analysis (Final Design)
# 
# Primary Method: Synthetic Difference-in-Differences (SDID)
# Backup Method: Augmented Synthetic Control (ASCM)
# 
# Spillover Analysis:
#   1. Thematic: High-prostitution BGs in Queens/Brooklyn (held out)
#   2. Spatial: Buffer zone around ORR corridor (held out)
#
# Robustness Checks:
#   - Alternative donor pools
#   - Pre-spike matching window
#   - Placebo timing tests
#   - Leave-one-out donors
#==============================================================================

library(here)
library(tidyverse)
library(lubridate)
library(sf)
library(synthdid)
library(augsynth)
library(future)
library(furrr)
library(patchwork)

# Parallel setup
n_workers <- min(parallel::detectCores() - 1, 8)
plan(multisession, workers = n_workers)
cat(sprintf("Using %d parallel workers\n", n_workers))

# Output directory
dir.create(here("output", "causal_final"), showWarnings = FALSE, recursive = TRUE)

#==============================================================================
# 1. PARAMETERS AND DATA LOADING
#==============================================================================

cat("\n=== SETUP ===\n")

# Load EDA objects
eda_objects <- readRDS(here("output", "eda_objects.rds"))
list2env(eda_objects, envir = .GlobalEnv)

# Key dates
op_start <- as.Date("2024-10-15")
op_end <- op_start + 89
treatment_month <- floor_date(op_start, "month")  # 2024-10-01

# Analysis window
data_start <- as.Date("2022-01-01")
data_end <- as.Date("2025-09-01")

# Pre-spike cutoff (for robustness - before the summer 2024 spike)
pre_spike_end <- as.Date("2024-06-01")

all_months <- seq(data_start, data_end, by = "month")
n_months <- length(all_months)
n_pre <- sum(all_months < treatment_month)
n_post <- n_months - n_pre

cat(sprintf("Analysis window: %s to %s\n", data_start, data_end))
cat(sprintf("Treatment month: %s\n", treatment_month))
cat(sprintf("Periods: %d total (%d pre, %d post)\n", n_months, n_pre, n_post))

#==============================================================================
# 2. DEFINE GEOGRAPHIC UNITS
#==============================================================================

cat("\n=== DEFINING GEOGRAPHIC UNITS ===\n")

# --- 2a. Treatment Zone ---
# Crimes within rosie_buffer
zone_monthly <- street_violent_crime %>%
  filter(date >= data_start) %>%
  st_filter(rosie_buffer, .predicate = st_intersects) %>%
  st_drop_geometry() %>%
  mutate(month = floor_date(date, "month")) %>%
  filter(month <= data_end) %>%
  count(month, name = "crime_count") %>%
  right_join(tibble(month = all_months), by = "month") %>%
  replace_na(list(crime_count = 0)) %>%
  arrange(month)

cat(sprintf("Zone monthly crime: mean=%.1f, range=%d-%d\n",
            mean(zone_monthly$crime_count),
            min(zone_monthly$crime_count),
            max(zone_monthly$crime_count)))

# --- 2b. Spatial Spillover Buffer ---
# Create buffer around the zone (e.g., 0.5 mile / ~800m beyond rosie_buffer)
spillover_distance <- 2640  # feet (0.5 mile in EPSG 2263)

spatial_spillover_buffer <- rosie_buffer %>%
  st_buffer(spillover_distance) %>%
  st_difference(rosie_buffer)  # Donut: spillover zone minus treatment zone

# BGs in spatial spillover zone
spatial_spillover_bgs <- nyc_bgs %>%
  st_filter(spatial_spillover_buffer, .predicate = st_intersects) %>%
  pull(geoid)

cat(sprintf("Spatial spillover BGs: %d\n", length(spatial_spillover_bgs)))

# --- 2c. Thematic Spillover: High-Prostitution BGs ---
# Identify prostitution hot spots in Queens and Brooklyn (outside zone)
# Using arrest data for prostitution-related offenses

# Prostitution offense codes (adjust based on your data)
prostitution_codes <- c(567, 568, 569, 570, 571, 572)  # Adjust as needed

# If you have arrests_sf with offense codes:
if (exists("arrests_sf")) {
  prostitution_arrests <- arrests_sf %>%
    filter(date >= data_start, date < treatment_month) %>%
    filter(pd_cd %in% prostitution_codes | 
             str_detect(tolower(ofns_desc), "prost|promot")) %>%
    st_join(nyc_bgs %>% select(geoid, boro_name)) %>%
    filter(boro_name %in% c("Queens", "Brooklyn")) %>%
    st_drop_geometry()
  
  # Count by BG
  prostitution_by_bg <- prostitution_arrests %>%
    count(geoid, name = "prostitution_arrests") %>%
    arrange(desc(prostitution_arrests))
  
  # Top prostitution BGs (e.g., top 50 or those with 5+ arrests)
  thematic_spillover_bgs <- prostitution_by_bg %>%
    filter(prostitution_arrests >= 5) %>%  # Threshold
    pull(geoid)
  
  # Remove any that overlap with zone or spatial spillover
  thematic_spillover_bgs <- setdiff(thematic_spillover_bgs, 
                                    c(spatial_spillover_bgs, 
                                      nyc_bgs %>% st_filter(rosie_buffer) %>% pull(geoid)))
  
  cat(sprintf("Thematic spillover BGs (prostitution hotspots): %d\n", 
              length(thematic_spillover_bgs)))
  
} else {
  # Fallback: use complaint data if arrests not available
  # Look for prostitution-related complaints
  prostitution_complaints <- complaints_sf %>%
    filter(date >= data_start, date < treatment_month) %>%
    filter(str_detect(tolower(ofns_desc), "prost|promot")) %>%
    st_join(nyc_bgs %>% select(geoid, boro_name)) %>%
    filter(boro_name %in% c("Queens", "Brooklyn")) %>%
    st_drop_geometry()
  
  prostitution_by_bg <- prostitution_complaints %>%
    count(geoid, name = "prostitution_count") %>%
    arrange(desc(prostitution_count))
  
  thematic_spillover_bgs <- prostitution_by_bg %>%
    slice_max(prostitution_count, n = 50) %>%
    pull(geoid)
  
  thematic_spillover_bgs <- setdiff(thematic_spillover_bgs,
                                    c(spatial_spillover_bgs,
                                      nyc_bgs %>% st_filter(rosie_buffer) %>% pull(geoid)))
  
  cat(sprintf("Thematic spillover BGs (prostitution hotspots): %d\n",
              length(thematic_spillover_bgs)))
}



# --- 2c. Thematic Spillover: High-Prostitution BGs ---
cat("\n=== IDENTIFYING PROSTITUTION HOT SPOTS ===\n")

# Prostitution arrests (PL 230 is the penal code for prostitution offenses)
# Filter to Queens (Q) and Brooklyn (K)
prostitution_arrests <- arrests_sf %>%
  filter(date >= data_start, date < treatment_month) %>%
  filter(arrest_boro %in% c("Q", "K")) %>%
  filter(str_detect(law_code, "PL 230"))

cat(sprintf("Prostitution arrests in Queens/Brooklyn (2022 - pre-treatment): %d\n", 
            nrow(prostitution_arrests)))

# Spatial join to BGs
prostitution_arrests_bg <- prostitution_arrests %>%
  st_join(nyc_bgs %>% select(geoid)) %>%
  st_drop_geometry() %>%
  filter(!is.na(geoid))

# Count by BG
prostitution_by_bg <- prostitution_arrests_bg %>%
  count(geoid, name = "prostitution_arrests") %>%
  arrange(desc(prostitution_arrests))

cat("\nTop 10 BGs by prostitution arrests:\n")
print(head(prostitution_by_bg, 10))

# Top prostitution BGs (those with 5+ arrests)
thematic_spillover_bgs <- prostitution_by_bg %>%
  filter(prostitution_arrests >= 5) %>%
  pull(geoid)

# If too few, take top 50 instead
if (length(thematic_spillover_bgs) < 20) {
  thematic_spillover_bgs <- prostitution_by_bg %>%
    slice_max(prostitution_arrests, n = 50, with_ties = FALSE) %>%
    pull(geoid)
}

# Remove any that overlap with zone or spatial spillover
thematic_spillover_bgs <- setdiff(thematic_spillover_bgs, 
                                  c(spatial_spillover_bgs, zone_bgs))

cat(sprintf("Thematic spillover BGs (after exclusions): %d\n", 
            length(thematic_spillover_bgs)))

# --- 2d. Zone BGs (excluded from donors) ---
zone_bgs <- nyc_bgs %>%
  st_filter(rosie_buffer, .predicate = st_intersects) %>%
  pull(geoid)

cat(sprintf("Zone BGs (excluded): %d\n", length(zone_bgs)))

# --- 2e. Donor Pool ---
# All NYC BGs EXCEPT: zone, spatial spillover, thematic spillover
all_excluded <- unique(c(zone_bgs, spatial_spillover_bgs, thematic_spillover_bgs))

donor_bgs <- nyc_bgs %>%
  filter(!geoid %in% all_excluded) %>%
  pull(geoid)

cat(sprintf("Donor BGs: %d\n", length(donor_bgs)))

#==============================================================================
# 3. BUILD PANEL DATA FOR ALL UNITS
#==============================================================================

cat("\n=== BUILDING PANEL DATA ===\n")

# Crime counts by BG
crime_by_bg <- street_violent_crime %>%
  filter(date >= data_start) %>%
  st_join(nyc_bgs %>% select(geoid)) %>%
  st_drop_geometry() %>%
  filter(!is.na(geoid)) %>%
  mutate(month = floor_date(date, "month")) %>%
  filter(month <= data_end) %>%
  count(geoid, month, name = "crime_count")

# Function to build panel for a set of BGs (aggregated)
build_aggregate_panel <- function(bg_ids, unit_name) {
  if (length(bg_ids) == 0) return(NULL)
  
  crime_by_bg %>%
    filter(geoid %in% bg_ids) %>%
    group_by(month) %>%
    summarize(crime_count = sum(crime_count), .groups = "drop") %>%
    right_join(tibble(month = all_months), by = "month") %>%
    replace_na(list(crime_count = 0)) %>%
    mutate(unit_id = unit_name) %>%
    arrange(month)
}

# Zone panel (already have this)
zone_panel <- zone_monthly %>%
  mutate(unit_id = "ZONE")

# Spatial spillover panel
spatial_spillover_panel <- build_aggregate_panel(spatial_spillover_bgs, "SPATIAL_SPILLOVER")

# Thematic spillover panel
thematic_spillover_panel <- build_aggregate_panel(thematic_spillover_bgs, "THEMATIC_SPILLOVER")

cat(sprintf("\nSpatial spillover monthly crime: mean=%.1f\n",
            mean(spatial_spillover_panel$crime_count)))
cat(sprintf("Thematic spillover monthly crime: mean=%.1f\n",
            mean(thematic_spillover_panel$crime_count)))

# --- Donor panel (individual BGs) ---
# Filter to BGs with sufficient activity
donor_stats <- crime_by_bg %>%
  filter(geoid %in% donor_bgs, month < treatment_month) %>%
  group_by(geoid) %>%
  summarize(
    mean_crime = mean(crime_count),
    total_crime = sum(crime_count),
    nonzero_months = sum(crime_count > 0),
    .groups = "drop"
  )

# Keep donors with at least some activity
min_total_crime <- 5
min_nonzero_pct <- 0.3

good_donors <- donor_stats %>%
  filter(
    total_crime >= min_total_crime,
    nonzero_months >= n_pre * min_nonzero_pct
  ) %>%
  pull(geoid)

cat(sprintf("\nDonors after filtering: %d (from %d)\n", 
            length(good_donors), length(donor_bgs)))

# Complete donor panel
donor_panel <- expand_grid(
  geoid = good_donors,
  month = all_months
) %>%
  left_join(crime_by_bg, by = c("geoid", "month")) %>%
  replace_na(list(crime_count = 0)) %>%
  rename(unit_id = geoid)

#==============================================================================
# 4. SELECT TOP CORRELATED DONORS
#==============================================================================

cat("\n=== SELECTING DONORS BY CORRELATION ===\n")

# Pre-treatment zone series
zone_pre <- zone_panel %>%
  filter(month < treatment_month) %>%
  arrange(month) %>%
  pull(crime_count)

# Compute correlation with zone for each donor
donor_cors <- donor_panel %>%
  filter(month < treatment_month) %>%
  group_by(unit_id) %>%
  arrange(month) %>%
  summarize(
    series = list(crime_count),
    .groups = "drop"
  ) %>%
  mutate(
    cor_zone = map_dbl(series, ~cor(.x, zone_pre, use = "complete.obs"))
  ) %>%
  filter(!is.na(cor_zone))

# Select top N donors
n_donors_select <- 150

top_donors <- donor_cors %>%
  slice_max(cor_zone, n = n_donors_select) %>%
  pull(unit_id)

cat(sprintf("Selected %d donors (correlation range: %.3f to %.3f)\n",
            length(top_donors),
            min(donor_cors$cor_zone[donor_cors$unit_id %in% top_donors]),
            max(donor_cors$cor_zone[donor_cors$unit_id %in% top_donors])))

#==============================================================================
# 5. HELPER FUNCTIONS
#==============================================================================

# Prepare matrix for synthdid
prepare_sdid_matrix <- function(treated_panel, donor_panel, 
                                donor_ids, outcome_col = "crime_count") {
  
  # Combine
  panel_combined <- donor_panel %>%
    filter(unit_id %in% donor_ids) %>%
    bind_rows(treated_panel)
  
  # Wide format
  panel_wide <- panel_combined %>%
    select(unit_id, month, !!sym(outcome_col)) %>%
    pivot_wider(names_from = month, values_from = !!sym(outcome_col), values_fill = 0) %>%
    column_to_rownames("unit_id")
  
  # Reorder: controls first, treated last
  mat <- as.matrix(panel_wide)
  treated_name <- unique(treated_panel$unit_id)
  treated_row <- mat[treated_name, , drop = FALSE]
  control_rows <- mat[rownames(mat) != treated_name, , drop = FALSE]
  
  Y <- rbind(control_rows, treated_row)
  N0 <- nrow(Y) - 1
  T0 <- sum(as.Date(colnames(Y)) < treatment_month)
  
  list(Y = Y, N0 = N0, T0 = T0)
}

# Run SDID with SE
run_sdid <- function(Y, N0, T0, name = "") {
  result <- list(name = name)
  
  tryCatch({
    est <- synthdid_estimate(Y, N0, T0)
    se <- sqrt(vcov(est, method = "placebo"))
    
    result$estimate <- as.numeric(est)
    result$se <- se
    result$ci_lower <- result$estimate - 1.96 * se
    result$ci_upper <- result$estimate + 1.96 * se
    result$p_value <- 2 * pnorm(-abs(result$estimate / se))
    result$weights <- attr(est, "weights")
    result$object <- est
    result$success <- TRUE
    
  }, error = function(e) {
    result$estimate <- NA
    result$se <- NA
    result$error <- e$message
    result$success <- FALSE
  })
  
  result
}

# Run ASCM
run_ascm <- function(panel_data, treated_name, name = "") {
  result <- list(name = name)
  
  # Prepare data
  ascm_data <- panel_data %>%
    mutate(
      treat = as.integer(unit_id == treated_name & month >= treatment_month),
      time = as.integer(factor(month))
    )
  
  tryCatch({
    fit <- augsynth(
      crime_count ~ treat,
      unit = unit_id,
      time = time,
      data = ascm_data,
      progfunc = "Ridge",
      scm = TRUE,
      CV = FALSE,
      lambda = 1e-3
    )
    
    summ <- summary(fit)
    
    result$estimate <- summ$average_att$Estimate
    result$se <- summ$average_att$Std.Error
    result$ci_lower <- result$estimate - 1.96 * result$se
    result$ci_upper <- result$estimate + 1.96 * result$se
    result$p_value <- summ$average_att$p_value
    result$object <- fit
    result$success <- TRUE
    
  }, error = function(e) {
    result$estimate <- NA
    result$se <- NA
    result$error <- e$message
    result$success <- FALSE
  })
  
  result
}

#==============================================================================
# 6. MAIN ANALYSIS: ZONE EFFECT
#==============================================================================

cat("\n")
cat("================================================================\n")
cat("MAIN ANALYSIS: TREATMENT EFFECT ON ZONE\n")
cat("================================================================\n")

# Prepare data
main_panel <- donor_panel %>%
  filter(unit_id %in% top_donors) %>%
  bind_rows(zone_panel)

main_matrix <- prepare_sdid_matrix(zone_panel, donor_panel, top_donors)

# --- SDID ---
cat("\nFitting SDID (primary method)...\n")
sdid_main <- run_sdid(main_matrix$Y, main_matrix$N0, main_matrix$T0, "Zone - SDID")

if (sdid_main$success) {
  cat(sprintf("  Estimate: %.2f (SE: %.2f)\n", sdid_main$estimate, sdid_main$se))
  cat(sprintf("  95%% CI: [%.2f, %.2f]\n", sdid_main$ci_lower, sdid_main$ci_upper))
  cat(sprintf("  p-value: %.4f\n", sdid_main$p_value))
  
  # Create plot, display, and save
  p_sdid_main <- synthdid_plot(sdid_main$object, se.method = "placebo")
  print(p_sdid_main)  # Display in RStudio
  
  ggsave(here("output", "causal_final", "01_sdid_zone.png"), p_sdid_main,
         width = 10, height = 6, dpi = 300)
}

# --- ASCM ---
cat("\nFitting ASCM (backup method)...\n")
ascm_main <- run_ascm(main_panel, "ZONE", "Zone - ASCM")

if (ascm_main$success) {
  cat(sprintf("  Estimate: %.2f (SE: %.2f)\n", ascm_main$estimate, ascm_main$se))
  cat(sprintf("  95%% CI: [%.2f, %.2f]\n", ascm_main$ci_lower, ascm_main$ci_upper))
  
  p_ascm_main <- plot(ascm_main$object) +
    labs(title = "Augmented SCM: Zone Effect") +
    theme_minimal()
  p_ascm_main
  ggsave(here("output", "causal_final", "02_ascm_zone.png"), p_ascm_main,
         width = 10, height = 6, dpi = 300)
}

#==============================================================================
# 7. SPATIAL SPILLOVER ANALYSIS
#==============================================================================

cat("\n")
cat("================================================================\n")
cat("SPATIAL SPILLOVER: BUFFER ZONE AROUND ORR\n")
cat("================================================================\n")

if (!is.null(spatial_spillover_panel) && nrow(spatial_spillover_panel) > 0) {
  
  # Use same donors (they don't overlap with spillover zone)
  spatial_panel <- donor_panel %>%
    filter(unit_id %in% top_donors) %>%
    bind_rows(spatial_spillover_panel)
  
  spatial_matrix <- prepare_sdid_matrix(spatial_spillover_panel, donor_panel, top_donors)
  
  # --- SDID ---
  cat("\nFitting SDID for spatial spillover...\n")
  sdid_spatial <- run_sdid(spatial_matrix$Y, spatial_matrix$N0, spatial_matrix$T0, 
                           "Spatial Spillover - SDID")
  
  if (sdid_spatial$success) {
    cat(sprintf("  Estimate: %.2f (SE: %.2f)\n", sdid_spatial$estimate, sdid_spatial$se))
    cat(sprintf("  95%% CI: [%.2f, %.2f]\n", sdid_spatial$ci_lower, sdid_spatial$ci_upper))
    cat(sprintf("  p-value: %.4f\n", sdid_spatial$p_value))
    
    # Interpretation
    if (sdid_spatial$estimate > 0 && sdid_spatial$p_value < 0.1) {
      cat("  >> Evidence of spatial DISPLACEMENT (crime increased in buffer)\n")
    } else if (sdid_spatial$estimate < 0 && sdid_spatial$p_value < 0.1) {
      cat("  >> Evidence of DIFFUSION OF BENEFITS (crime decreased in buffer)\n")
    } else {
      cat("  >> No significant spatial spillover detected\n")
    }
    
    png(here("output", "causal_final", "03_sdid_spatial_spillover.png"),
        width = 10, height = 6, units = "in", res = 300)
    synthdid_plot(sdid_spatial$object, se.method = "placebo")
    dev.off()
  }
  
  # --- ASCM ---
  cat("\nFitting ASCM for spatial spillover...\n")
  ascm_spatial <- run_ascm(spatial_panel, "SPATIAL_SPILLOVER", "Spatial Spillover - ASCM")
  
  if (ascm_spatial$success) {
    cat(sprintf("  Estimate: %.2f (SE: %.2f)\n", ascm_spatial$estimate, ascm_spatial$se))
    
    p_ascm_spatial <- plot(ascm_spatial$object) +
      labs(title = "Augmented SCM: Spatial Spillover (Buffer Zone)") +
      theme_minimal()
    ggsave(here("output", "causal_final", "04_ascm_spatial_spillover.png"), p_ascm_spatial,
           width = 10, height = 6, dpi = 300)
  }
}

#==============================================================================
# 8. THEMATIC SPILLOVER ANALYSIS (PROSTITUTION HOT SPOTS)
#==============================================================================

cat("\n")
cat("================================================================\n")
cat("THEMATIC SPILLOVER: PROSTITUTION HOT SPOTS (QUEENS/BROOKLYN)\n")
cat("================================================================\n")

if (!is.null(thematic_spillover_panel) && nrow(thematic_spillover_panel) > 0) {
  
  # Use same donors
  thematic_panel <- donor_panel %>%
    filter(unit_id %in% top_donors) %>%
    bind_rows(thematic_spillover_panel)
  
  thematic_matrix <- prepare_sdid_matrix(thematic_spillover_panel, donor_panel, top_donors)
  
  # --- SDID ---
  cat("\nFitting SDID for thematic spillover...\n")
  sdid_thematic <- run_sdid(thematic_matrix$Y, thematic_matrix$N0, thematic_matrix$T0,
                            "Thematic Spillover - SDID")
  
  if (sdid_thematic$success) {
    cat(sprintf("  Estimate: %.2f (SE: %.2f)\n", sdid_thematic$estimate, sdid_thematic$se))
    cat(sprintf("  95%% CI: [%.2f, %.2f]\n", sdid_thematic$ci_lower, sdid_thematic$ci_upper))
    cat(sprintf("  p-value: %.4f\n", sdid_thematic$p_value))
    
    # Interpretation
    if (sdid_thematic$estimate > 0 && sdid_thematic$p_value < 0.1) {
      cat("  >> Evidence of THEMATIC DISPLACEMENT (crime shifted to other prostitution areas)\n")
    } else if (sdid_thematic$estimate < 0 && sdid_thematic$p_value < 0.1) {
      cat("  >> Crime DECREASED in other prostitution hot spots too\n")
    } else {
      cat("  >> No significant thematic spillover detected\n")
    }
    
    png(here("output", "causal_final", "05_sdid_thematic_spillover.png"),
        width = 10, height = 6, units = "in", res = 300)
    synthdid_plot(sdid_thematic$object, se.method = "placebo")
    dev.off()
  }
  
  # --- ASCM ---
  cat("\nFitting ASCM for thematic spillover...\n")
  ascm_thematic <- run_ascm(thematic_panel, "THEMATIC_SPILLOVER", "Thematic Spillover - ASCM")
  
  if (ascm_thematic$success) {
    cat(sprintf("  Estimate: %.2f (SE: %.2f)\n", ascm_thematic$estimate, ascm_thematic$se))
    
    p_ascm_thematic <- plot(ascm_thematic$object) +
      labs(title = "Augmented SCM: Thematic Spillover (Prostitution Hot Spots)") +
      theme_minimal()
    ggsave(here("output", "causal_final", "06_ascm_thematic_spillover.png"), p_ascm_thematic,
           width = 10, height = 6, dpi = 300)
  }
}

#==============================================================================
# 9. ROBUSTNESS CHECKS
#==============================================================================

cat("\n")
cat("================================================================\n")
cat("ROBUSTNESS CHECKS\n")
cat("================================================================\n")

robustness_results <- list()

# --- 9a. Pre-Spike Matching ---
cat("\n--- Robustness 1: Pre-Spike Matching Window ---\n")
cat("Optimizing weights on Jan 2022 - May 2024 only (before summer spike)\n")

# Restrict to pre-spike period for weight estimation
pre_spike_months <- all_months[all_months < pre_spike_end]
n_pre_spike <- length(pre_spike_months)

tryCatch({
  # Create modified panel with only pre-spike for matching
  # But include all periods for outcome
  
  # This requires manual weight estimation - use correlation on pre-spike
  zone_pre_spike <- zone_panel %>%
    filter(month < pre_spike_end) %>%
    arrange(month) %>%
    pull(crime_count)
  
  donor_cors_prespike <- donor_panel %>%
    filter(month < pre_spike_end) %>%
    group_by(unit_id) %>%
    arrange(month) %>%
    summarize(
      series = list(crime_count),
      .groups = "drop"
    ) %>%
    mutate(
      cor_zone = map_dbl(series, ~cor(.x, zone_pre_spike, use = "complete.obs"))
    ) %>%
    filter(!is.na(cor_zone))
  
  top_donors_prespike <- donor_cors_prespike %>%
    slice_max(cor_zone, n = n_donors_select) %>%
    pull(unit_id)
  
  prespike_matrix <- prepare_sdid_matrix(zone_panel, donor_panel, top_donors_prespike)
  sdid_prespike <- run_sdid(prespike_matrix$Y, prespike_matrix$N0, prespike_matrix$T0,
                            "Pre-Spike Matching")
  
  if (sdid_prespike$success) {
    cat(sprintf("  Estimate: %.2f (SE: %.2f)\n", sdid_prespike$estimate, sdid_prespike$se))
    robustness_results$prespike <- sdid_prespike
  }
  
}, error = function(e) {
  cat(sprintf("  Failed: %s\n", e$message))
})

# --- 9b. Alternative Donor Pool: Queens Only ---
cat("\n--- Robustness 2: Queens-Only Donors ---\n")

tryCatch({
  queens_donors <- nyc_bgs %>%
    filter(boro_name == "Queens", !geoid %in% all_excluded) %>%
    pull(geoid)
  
  queens_good <- donor_stats %>%
    filter(geoid %in% queens_donors,
           total_crime >= min_total_crime,
           nonzero_months >= n_pre * min_nonzero_pct) %>%
    pull(geoid)
  
  cat(sprintf("  Queens donors available: %d\n", length(queens_good)))
  
  if (length(queens_good) >= 50) {
    # Select by correlation
    queens_cors <- donor_cors %>%
      filter(unit_id %in% queens_good) %>%
      slice_max(cor_zone, n = min(100, length(queens_good)))
    
    queens_matrix <- prepare_sdid_matrix(zone_panel, donor_panel, queens_cors$unit_id)
    sdid_queens <- run_sdid(queens_matrix$Y, queens_matrix$N0, queens_matrix$T0,
                            "Queens-Only Donors")
    
    if (sdid_queens$success) {
      cat(sprintf("  Estimate: %.2f (SE: %.2f)\n", sdid_queens$estimate, sdid_queens$se))
      robustness_results$queens_only <- sdid_queens
    }
  }
  
}, error = function(e) {
  cat(sprintf("  Failed: %s\n", e$message))
})

# --- 9c. Alternative Donor Pool: Citywide (exclude Queens) ---
cat("\n--- Robustness 3: Non-Queens Donors ---\n")

tryCatch({
  nonqueens_donors <- nyc_bgs %>%
    filter(boro_name != "Queens", !geoid %in% all_excluded) %>%
    pull(geoid)
  
  nonqueens_good <- donor_stats %>%
    filter(geoid %in% nonqueens_donors,
           total_crime >= min_total_crime,
           nonzero_months >= n_pre * min_nonzero_pct) %>%
    pull(geoid)
  
  cat(sprintf("  Non-Queens donors available: %d\n", length(nonqueens_good)))
  
  if (length(nonqueens_good) >= 50) {
    nonqueens_cors <- donor_cors %>%
      filter(unit_id %in% nonqueens_good) %>%
      slice_max(cor_zone, n = min(100, length(nonqueens_good)))
    
    nonqueens_matrix <- prepare_sdid_matrix(zone_panel, donor_panel, nonqueens_cors$unit_id)
    sdid_nonqueens <- run_sdid(nonqueens_matrix$Y, nonqueens_matrix$N0, nonqueens_matrix$T0,
                               "Non-Queens Donors")
    
    if (sdid_nonqueens$success) {
      cat(sprintf("  Estimate: %.2f (SE: %.2f)\n", sdid_nonqueens$estimate, sdid_nonqueens$se))
      robustness_results$nonqueens <- sdid_nonqueens
    }
  }
  
}, error = function(e) {
  cat(sprintf("  Failed: %s\n", e$message))
})

# --- 9d. Placebo Timing Test ---
cat("\n--- Robustness 4: Placebo Timing (Fake Treatment at 2023-10-01) ---\n")

tryCatch({
  placebo_treatment <- as.Date("2023-10-01")
  
  # Subset data to end before real treatment
  placebo_months <- all_months[all_months < treatment_month]
  
  placebo_zone <- zone_panel %>%
    filter(month < treatment_month)
  
  placebo_donors <- donor_panel %>%
    filter(month < treatment_month, unit_id %in% top_donors)
  
  # Prepare matrix
  placebo_wide <- placebo_donors %>%
    bind_rows(placebo_zone %>% rename(unit_id = unit_id)) %>%
    select(unit_id, month, crime_count) %>%
    pivot_wider(names_from = month, values_from = crime_count, values_fill = 0) %>%
    column_to_rownames("unit_id")
  
  placebo_mat <- as.matrix(placebo_wide)
  zone_row <- placebo_mat["ZONE", , drop = FALSE]
  control_rows <- placebo_mat[rownames(placebo_mat) != "ZONE", , drop = FALSE]
  Y_placebo <- rbind(control_rows, zone_row)
  
  N0_placebo <- nrow(Y_placebo) - 1
  T0_placebo <- sum(as.Date(colnames(Y_placebo)) < placebo_treatment)
  
  sdid_placebo <- run_sdid(Y_placebo, N0_placebo, T0_placebo, "Placebo (Oct 2023)")
  
  if (sdid_placebo$success) {
    cat(sprintf("  Estimate: %.2f (SE: %.2f)\n", sdid_placebo$estimate, sdid_placebo$se))
    cat(sprintf("  p-value: %.4f\n", sdid_placebo$p_value))
    
    if (sdid_placebo$p_value > 0.1) {
      cat("  >> Good: No significant placebo effect (supports causal interpretation)\n")
    } else {
      cat("  >> Warning: Significant placebo effect (pre-trends may be violated)\n")
    }
    
    robustness_results$placebo <- sdid_placebo
  }
  
}, error = function(e) {
  cat(sprintf("  Failed: %s\n", e$message))
})

# --- 9e. Leave-One-Out: Drop Top Weighted Donors ---
cat("\n--- Robustness 5: Leave-One-Out (Drop Top 5 Weighted Donors) ---\n")

tryCatch({
  if (sdid_main$success && !is.null(sdid_main$weights)) {
    # Get donor weights
    omega <- sdid_main$weights$omega
    donor_names <- rownames(main_matrix$Y)[1:main_matrix$N0]
    
    # Find top 5 weighted donors
    weight_df <- tibble(
      unit_id = donor_names,
      weight = as.vector(omega)
    ) %>%
      arrange(desc(weight)) %>%
      slice_head(n = 5)
    
    cat("  Top 5 weighted donors:\n")
    print(weight_df)
    
    # Re-estimate without them
    loo_donors <- setdiff(top_donors, weight_df$unit_id)
    
    loo_matrix <- prepare_sdid_matrix(zone_panel, donor_panel, loo_donors)
    sdid_loo <- run_sdid(loo_matrix$Y, loo_matrix$N0, loo_matrix$T0, "Leave-One-Out")
    
    if (sdid_loo$success) {
      cat(sprintf("  Estimate: %.2f (SE: %.2f)\n", sdid_loo$estimate, sdid_loo$se))
      cat(sprintf("  Change from main: %.2f\n", sdid_loo$estimate - sdid_main$estimate))
      robustness_results$loo <- sdid_loo
    }
  }
  
}, error = function(e) {
  cat(sprintf("  Failed: %s\n", e$message))
})

#==============================================================================
# 10. RESULTS SUMMARY
#==============================================================================

cat("\n")
cat("================================================================\n")
cat("RESULTS SUMMARY\n")
cat("================================================================\n")

# Build results table
results_list <- list(
  list(name = "Zone (Main)", method = "SDID", result = sdid_main),
  list(name = "Zone (Main)", method = "ASCM", result = ascm_main)
)

if (exists("sdid_spatial") && sdid_spatial$success) {
  results_list <- c(results_list, list(
    list(name = "Spatial Spillover", method = "SDID", result = sdid_spatial),
    list(name = "Spatial Spillover", method = "ASCM", result = ascm_spatial)
  ))
}

if (exists("sdid_thematic") && sdid_thematic$success) {
  results_list <- c(results_list, list(
    list(name = "Thematic Spillover", method = "SDID", result = sdid_thematic),
    list(name = "Thematic Spillover", method = "ASCM", result = ascm_thematic)
  ))
}

# Add robustness checks
for (rob_name in names(robustness_results)) {
  results_list <- c(results_list, list(
    list(name = paste0("Robustness: ", rob_name), method = "SDID", 
         result = robustness_results[[rob_name]])
  ))
}

# Create table
results_table <- map_dfr(results_list, function(x) {
  tibble(
    Analysis = x$name,
    Method = x$method,
    Estimate = x$result$estimate,
    SE = x$result$se,
    CI_Lower = x$result$ci_lower,
    CI_Upper = x$result$ci_upper,
    P_Value = x$result$p_value
  )
}) %>%
  mutate(
    Significant = case_when(
      P_Value < 0.01 ~ "***",
      P_Value < 0.05 ~ "**",
      P_Value < 0.1 ~ "*",
      TRUE ~ ""
    )
  )

cat("\n")
print(results_table, n = 30)

write_csv(results_table, here("output", "causal_final", "results_summary.csv"))

#==============================================================================
# 11. VISUALIZATION
#==============================================================================

cat("\n=== GENERATING SUMMARY PLOTS ===\n")

# Forest plot - Main effects
p_forest_main <- results_table %>%
  filter(str_detect(Analysis, "Zone|Spillover"), !str_detect(Analysis, "Robustness")) %>%
  mutate(Analysis = paste(Analysis, Method, sep = " - ")) %>%
  ggplot(aes(x = Estimate, y = reorder(Analysis, Estimate))) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title = "Treatment Effects: Zone and Spillovers",
    subtitle = "Point estimates with 95% CI. Negative = crime reduction.",
    x = "Estimated Effect (monthly crimes)",
    y = ""
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10))

ggsave(here("output", "causal_final", "07_forest_main.png"), p_forest_main,
       width = 10, height = 6, dpi = 300)

# Forest plot - Robustness
p_forest_robust <- results_table %>%
  filter(str_detect(Analysis, "Zone \\(Main\\)|Robustness"), Method == "SDID") %>%
  mutate(Analysis = str_remove(Analysis, "Robustness: ")) %>%
  ggplot(aes(x = Estimate, y = reorder(Analysis, Estimate))) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title = "Robustness Checks (SDID)",
    subtitle = "Main estimate vs. alternative specifications",
    x = "Estimated Effect (monthly crimes)",
    y = ""
  ) +
  theme_minimal()

ggsave(here("output", "causal_final", "08_forest_robustness.png"), p_forest_robust,
       width = 10, height = 6, dpi = 300)

# Time series: Zone vs Spillovers
ts_data <- bind_rows(
  zone_panel %>% mutate(unit = "Treatment Zone"),
  spatial_spillover_panel %>% mutate(unit = "Spatial Spillover (Buffer)"),
  thematic_spillover_panel %>% mutate(unit = "Thematic Spillover (Prostitution BGs)")
)

p_timeseries <- ggplot(ts_data, aes(x = month, y = crime_count, color = unit)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1) +
  geom_vline(xintercept = treatment_month, linetype = "dashed") +
  annotate("rect", xmin = treatment_month, xmax = max(all_months),
           ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "red") +
  scale_color_manual(values = c(
    "Treatment Zone" = "darkred",
    "Spatial Spillover (Buffer)" = "steelblue",
    "Thematic Spillover (Prostitution BGs)" = "forestgreen"
  )) +
  labs(
    title = "Crime Trends: Treatment Zone vs. Spillover Areas",
    subtitle = "Shaded area = post-treatment period",
    x = "Month", y = "Crime Count", color = ""
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(here("output", "causal_final", "09_timeseries_comparison.png"), p_timeseries,
       width = 12, height = 6, dpi = 300)

#==============================================================================
# 12. SAVE RESULTS
#==============================================================================

# Save all results
all_results <- list(
  main = list(sdid = sdid_main, ascm = ascm_main),
  spatial_spillover = if(exists("sdid_spatial")) list(sdid = sdid_spatial, ascm = ascm_spatial) else NULL,
  thematic_spillover = if(exists("sdid_thematic")) list(sdid = sdid_thematic, ascm = ascm_thematic) else NULL,
  robustness = robustness_results,
  config = list(
    n_donors = length(top_donors),
    treatment_month = treatment_month,
    data_start = data_start,
    data_end = data_end,
    n_spatial_spillover_bgs = length(spatial_spillover_bgs),
    n_thematic_spillover_bgs = length(thematic_spillover_bgs)
  )
)

saveRDS(all_results, here("output", "causal_final", "all_results.rds"))

# Save panels
saveRDS(list(
  zone = zone_panel,
  spatial_spillover = spatial_spillover_panel,
  thematic_spillover = thematic_spillover_panel,
  donors = donor_panel %>% filter(unit_id %in% top_donors)
), here("output", "causal_final", "panels.rds"))

# Clean up
plan(sequential)

cat("\n")
cat("================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("================================================================\n")

cat(sprintf("\nOutputs saved to: %s\n", here("output", "causal_final")))
list.files(here("output", "causal_final")) %>%
  paste(" -", .) %>%
  cat(sep = "\n")

cat("\n\n=== KEY FINDINGS ===\n")
if (sdid_main$success) {
  cat(sprintf("Main Effect (SDID): %.2f crimes/month (p = %.4f)\n", 
              sdid_main$estimate, sdid_main$p_value))
}
if (exists("sdid_spatial") && sdid_spatial$success) {
  cat(sprintf("Spatial Spillover: %.2f crimes/month (p = %.4f)\n",
              sdid_spatial$estimate, sdid_spatial$p_value))
}
if (exists("sdid_thematic") && sdid_thematic$success) {
  cat(sprintf("Thematic Spillover: %.2f crimes/month (p = %.4f)\n",
              sdid_thematic$estimate, sdid_thematic$p_value))
}