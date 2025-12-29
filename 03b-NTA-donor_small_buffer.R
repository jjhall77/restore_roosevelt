#smaller buffers

#==============================================================================
# Operation Restore Roosevelt - FINAL Causal Analysis
# 
# Primary Methods: gsynth + Synthetic DiD (NTA-level)
# Spillover: Spatial buffer + Prostitution hot spots
# 
# Based on exploratory analysis findings:
#   - NTA-level provides better scale match than BGs
#   - gsynth and SDID both show significant negative effects
#   - ASCM confirms direction but with wider uncertainty
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

# Parallel setup
n_workers <- min(parallel::detectCores() - 1, 8)
plan(multisession, workers = n_workers)

# Output directory
dir.create(here("output", "causal_final"), showWarnings = FALSE, recursive = TRUE)

#==============================================================================
# 1. SETUP AND DATA LOADING
#==============================================================================

cat("\n=== SETUP ===\n")

# Load EDA objects
eda_objects <- readRDS(here("output", "eda_objects.rds"))
list2env(eda_objects, envir = .GlobalEnv)

# Load NTA shapefile
nynta <- nynta %>% st_transform(2263)
nta_id_col <- intersect(c("ntacode", "nta2020", "NTACode", "NTA2020"), names(nynta))[1]


# Key dates
op_start <- as.Date("2024-10-15")
op_end <- op_start + 89  # Jan 12, 2025
treatment_month <- floor_date(op_start, "month")

# Analysis window
data_start <- as.Date("2022-01-01")
data_end <- as.Date("2025-09-30")

all_months <- seq(data_start, data_end, by = "month")
n_months <- length(all_months)
n_pre <- sum(all_months < treatment_month)
n_post <- n_months - n_pre

cat(sprintf("Treatment: %s (Op ran %s to %s)\n", treatment_month, op_start, op_end))
cat(sprintf("Periods: %d total (%d pre, %d post)\n", n_months, n_pre, n_post))

#==============================================================================
# 2. DEFINE GEOGRAPHIC UNITS
#==============================================================================

cat("\n=== DEFINING GEOGRAPHIC UNITS ===\n")

# --- Treatment Zone ---
# Crimes within rosie_buffer (unchanged)
cat("Treatment zone: rosie_buffer\n")

# --- Spatial Spillover: Direct buffer donut (NOT NTA-based) ---
spillover_distance <- 2640  # 0.5 mile in feet
spatial_buffer <- rosie_buffer %>%
  st_buffer(spillover_distance) %>%
  st_difference(rosie_buffer)  # Donut shape

cat(sprintf("Spatial spillover: %.2f sq mi buffer donut\n", 
            as.numeric(st_area(spatial_buffer)) / 27878400))  # sq ft to sq mi

# --- Thematic Spillover: High-prostitution BGs (NOT NTAs) ---
# Prostitution arrests by BG
prostitution_arrests <- arrests_sf %>%
  filter(date >= data_start, date < treatment_month) %>%
  filter(arrest_boro %in% c("Q", "K")) %>%
  filter(str_detect(law_code, "PL 230"))

pros_by_bg <- prostitution_arrests %>%
  st_join(nyc_bgs %>% select(geoid)) %>%
  st_drop_geometry() %>%
  filter(!is.na(geoid)) %>%
  count(geoid, name = "pros_arrests") %>%
  arrange(desc(pros_arrests))

cat(sprintf("Prostitution arrests in Q/K: %d\n", nrow(prostitution_arrests)))
cat("\nTop 10 prostitution BGs:\n")
print(head(pros_by_bg, 10))

# Identify BGs that overlap with zone or spatial buffer (to exclude)
zone_bgs <- nyc_bgs %>%
  st_filter(rosie_buffer, .predicate = st_intersects) %>%
  pull(geoid)

buffer_bgs <- nyc_bgs %>%
  st_filter(spatial_buffer, .predicate = st_intersects) %>%
  pull(geoid)

# High-prostitution BGs (5+ arrests), excluding zone and buffer
thematic_spillover_bgs <- pros_by_bg %>%
  filter(pros_arrests >= 5) %>%
  filter(!geoid %in% c(zone_bgs, buffer_bgs)) %>%
  pull(geoid)

cat(sprintf("Thematic spillover BGs: %d high-prostitution BGs\n", length(thematic_spillover_bgs)))

# Create combined geometry for thematic spillover
thematic_spillover_area <- nyc_bgs %>%
  filter(geoid %in% thematic_spillover_bgs) %>%
  st_union()

# --- Identify NTAs that intersect zone (for reference) ---
zone_ntas <- nynta %>%
  st_filter(rosie_buffer, .predicate = st_intersects) %>%
  pull(nta_id)

cat(sprintf("NTAs intersecting zone (for reference): %d\n", length(zone_ntas)))

# --- Donor NTAs ---
# All NTAs are potential donors; we'll subtract crimes from excluded areas
donor_ntas <- nynta %>%
  pull(nta_id)

cat(sprintf("Potential donor NTAs: %d (will subtract excluded area crimes)\n", length(donor_ntas)))

#==============================================================================
# 3. BUILD PANEL DATA
#==============================================================================

cat("\n=== BUILDING PANEL DATA ===\n")

# --- Zone panel: crimes within rosie_buffer ---
zone_panel <- street_violent_crime %>%
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

cat(sprintf("Zone monthly crime: mean=%.1f, range=%d-%d\n",
            mean(zone_panel$crime_count),
            min(zone_panel$crime_count),
            max(zone_panel$crime_count)))

# --- Spatial spillover panel: crimes within buffer donut ---
spatial_panel <- street_violent_crime %>%
  filter(date >= data_start) %>%
  st_filter(spatial_buffer, .predicate = st_intersects) %>%
  st_drop_geometry() %>%
  mutate(month = floor_date(date, "month")) %>%
  filter(month <= data_end) %>%
  count(month, name = "crime_count") %>%
  right_join(tibble(month = all_months), by = "month") %>%
  replace_na(list(crime_count = 0)) %>%
  mutate(nta_id = "SPATIAL_SPILLOVER") %>%
  arrange(month)

cat(sprintf("Spatial spillover monthly crime: mean=%.1f, range=%d-%d\n",
            mean(spatial_panel$crime_count),
            min(spatial_panel$crime_count),
            max(spatial_panel$crime_count)))

# --- Thematic spillover panel: crimes within high-prostitution BGs ---
thematic_panel <- street_violent_crime %>%
  filter(date >= data_start) %>%
  st_filter(thematic_spillover_area, .predicate = st_intersects) %>%
  st_drop_geometry() %>%
  mutate(month = floor_date(date, "month")) %>%
  filter(month <= data_end) %>%
  count(month, name = "crime_count") %>%
  right_join(tibble(month = all_months), by = "month") %>%
  replace_na(list(crime_count = 0)) %>%
  mutate(nta_id = "THEMATIC_SPILLOVER") %>%
  arrange(month)

cat(sprintf("Thematic spillover monthly crime: mean=%.1f, range=%d-%d\n",
            mean(thematic_panel$crime_count),
            min(thematic_panel$crime_count),
            max(thematic_panel$crime_count)))

# --- Identify crimes to exclude from donor pool ---
cat("\nIdentifying crimes in excluded areas...\n")

crimes_in_zone <- street_violent_crime %>%
  filter(date >= data_start) %>%
  st_filter(rosie_buffer, .predicate = st_intersects) %>%
  pull(cmplnt_num)

crimes_in_buffer <- street_violent_crime %>%
  filter(date >= data_start) %>%
  st_filter(spatial_buffer, .predicate = st_intersects) %>%
  pull(cmplnt_num)

crimes_in_thematic <- street_violent_crime %>%
  filter(date >= data_start) %>%
  st_filter(thematic_spillover_area, .predicate = st_intersects) %>%
  pull(cmplnt_num)

excluded_crimes <- unique(c(crimes_in_zone, crimes_in_buffer, crimes_in_thematic))

cat(sprintf("Crimes excluded from donor pool:\n"))
cat(sprintf("  Zone: %d\n", length(crimes_in_zone)))
cat(sprintf("  Spatial buffer: %d\n", length(crimes_in_buffer)))
cat(sprintf("  Thematic BGs: %d\n", length(crimes_in_thematic)))
cat(sprintf("  Total unique excluded: %d\n", length(excluded_crimes)))

# --- Build donor panel: NTA-level crimes EXCLUDING zone, buffer, and thematic areas ---
# Tag each crime with its NTA
crime_with_nta <- street_violent_crime %>%
  filter(date >= data_start) %>%
  st_join(nynta %>% select(nta_id)) %>%
  mutate(month = floor_date(date, "month")) %>%
  filter(month <= data_end)

# Count crimes by NTA, excluding those in zone/buffer/thematic areas
crime_by_nta <- crime_with_nta %>%
  st_drop_geometry() %>%
  filter(!is.na(nta_id)) %>%
  filter(!cmplnt_num %in% excluded_crimes) %>%
  count(nta_id, month, name = "crime_count")

# Create complete panel
donor_panel <- expand_grid(nta_id = donor_ntas, month = all_months) %>%
  left_join(crime_by_nta, by = c("nta_id", "month")) %>%
  replace_na(list(crime_count = 0))

# Filter donors to those with sufficient activity
donor_stats <- donor_panel %>%
  filter(month < treatment_month) %>%
  group_by(nta_id) %>%
  summarize(
    total_crime = sum(crime_count),
    nonzero_months = sum(crime_count > 0),
    mean_crime = mean(crime_count),
    .groups = "drop"
  )

good_donors <- donor_stats %>%
  filter(total_crime >= 20, nonzero_months >= n_pre * 0.5) %>%
  pull(nta_id)

cat(sprintf("\nFiltered donors: %d NTAs (from %d)\n", length(good_donors), length(donor_ntas)))

# Summary stats
cat(sprintf("\nMonthly crime means:\n"))
cat(sprintf("  Zone: %.1f\n", mean(zone_panel$crime_count)))
cat(sprintf("  Spatial spillover (buffer): %.1f\n", mean(spatial_panel$crime_count)))
cat(sprintf("  Thematic spillover (pros BGs): %.1f\n", mean(thematic_panel$crime_count)))
cat(sprintf("  Donor NTAs (median): %.1f\n", 
            median(donor_stats %>% filter(nta_id %in% good_donors) %>% pull(mean_crime))))
#=======================================================
# 4. HELPER FUNCTIONS
#==============================================================================

# Prepare SDID matrix
prepare_sdid_matrix <- function(treated_panel, donor_panel, donor_ids) {
  combined <- donor_panel %>%
    filter(nta_id %in% donor_ids) %>%
    bind_rows(treated_panel)
  
  wide <- combined %>%
    select(nta_id, month, crime_count) %>%
    pivot_wider(names_from = month, values_from = crime_count, values_fill = 0) %>%
    column_to_rownames("nta_id")
  
  mat <- as.matrix(wide)
  treated_name <- unique(treated_panel$nta_id)
  treated_row <- mat[treated_name, , drop = FALSE]
  control_rows <- mat[rownames(mat) != treated_name, , drop = FALSE]
  
  Y <- rbind(control_rows, treated_row)
  N0 <- nrow(Y) - 1
  T0 <- sum(as.Date(colnames(Y)) < treatment_month)
  
  list(Y = Y, N0 = N0, T0 = T0)
}

# Run gsynth
run_gsynth <- function(panel_data, treated_name, name = "") {
  panel <- panel_data %>%
    mutate(
      treat = as.integer(nta_id == treated_name & month >= treatment_month),
      log_crime = log(crime_count + 0.5)
    )
  
  result <- list(name = name)
  
  tryCatch({
    gs <- gsynth(
      crime_count ~ treat,
      data = panel,
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
    
    se <- if (!is.null(gs$est.avg) && is.matrix(gs$est.avg)) {
      gs$est.avg[1, "S.E."]
    } else NA_real_
    
    result$estimate <- gs$att.avg
    result$se <- se
    result$ci_lower <- result$estimate - 1.96 * se
    result$ci_upper <- result$estimate + 1.96 * se
    result$p_value <- if (!is.na(se)) 2 * pnorm(-abs(result$estimate / se)) else NA
    result$factors <- gs$r.cv
    result$object <- gs
    result$success <- TRUE
    
  }, error = function(e) {
    result$success <- FALSE
    result$error <- e$message
  })
  
  result
}

# Run SDID
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
    result$object <- est
    result$success <- TRUE
    
  }, error = function(e) {
    result$success <- FALSE
    result$error <- e$message
  })
  
  result
}

library(ggplotify)

#==============================================================================
# 5. MAIN ANALYSIS: TREATMENT EFFECT ON ZONE
#==============================================================================

cat("\n")
cat("================================================================\n")
cat("MAIN ANALYSIS: TREATMENT EFFECT ON ZONE\n")
cat("================================================================\n")

# Combine zone with donors
main_panel <- donor_panel %>%
  filter(nta_id %in% good_donors) %>%
  bind_rows(zone_panel)

# --- Gsynth ---
cat("\nFitting gsynth...\n")
gs_zone <- run_gsynth(main_panel, "ZONE", "Zone - Gsynth")

if (gs_zone$success) {
  cat(sprintf("  ATT: %.2f (SE: %.2f)\n", gs_zone$estimate, gs_zone$se))
  cat(sprintf("  95%% CI: [%.2f, %.2f]\n", gs_zone$ci_lower, gs_zone$ci_upper))
  cat(sprintf("  p-value: %.4f\n", gs_zone$p_value))
  cat(sprintf("  Factors: %d\n", gs_zone$factors))
  
  # Gsynth counterfactual plot - assign output directly
  p_gs_zone_cf <- plot(gs_zone$object, type = "counterfactual", main = "Gsynth: Zone vs Counterfactual")
  print(p_gs_zone_cf)
  ggsave(here("output", "causal_final", "gs_zone_counterfactual.png"), p_gs_zone_cf,
         width = 10, height = 6, dpi = 300)
  
  # Gsynth gap plot
  p_gs_zone_gap <- plot(gs_zone$object, type = "gap", main = "Gsynth: Treatment Effect (Zone)")
  print(p_gs_zone_gap)
  ggsave(here("output", "causal_final", "gs_zone_gap.png"), p_gs_zone_gap,
         width = 10, height = 6, dpi = 300)
  
  cat("  Gsynth plots saved.\n")
}

# --- SDID ---
cat("\nFitting SDID...\n")
main_matrix <- prepare_sdid_matrix(zone_panel, donor_panel, good_donors)
sdid_zone <- run_sdid(main_matrix$Y, main_matrix$N0, main_matrix$T0, "Zone - SDID")

if (sdid_zone$success) {
  cat(sprintf("  ATT: %.2f (SE: %.2f)\n", sdid_zone$estimate, sdid_zone$se))
  cat(sprintf("  95%% CI: [%.2f, %.2f]\n", sdid_zone$ci_lower, sdid_zone$ci_upper))
  cat(sprintf("  p-value: %.4f\n", sdid_zone$p_value))
  
  # SDID plot
  p_sdid_zone <- synthdid_plot(sdid_zone$object, se.method = "placebo")
  print(p_sdid_zone)
  ggsave(here("output", "causal_final", "sdid_zone.png"), p_sdid_zone,
         width = 10, height = 6, dpi = 300)
  
  cat("  SDID plot saved.\n")
}

#==============================================================================
# 6. SPATIAL SPILLOVER ANALYSIS
#==============================================================================

cat("\n")
cat("================================================================\n")
cat("SPATIAL SPILLOVER: BUFFER AROUND ZONE\n")
cat("================================================================\n")

spatial_main_panel <- donor_panel %>%
  filter(nta_id %in% good_donors) %>%
  bind_rows(spatial_panel)

# --- Gsynth ---
cat("\nFitting gsynth for spatial spillover...\n")
gs_spatial <- run_gsynth(spatial_main_panel, "SPATIAL_SPILLOVER", "Spatial - Gsynth")

if (gs_spatial$success) {
  cat(sprintf("  ATT: %.2f (SE: %.2f)\n", gs_spatial$estimate, gs_spatial$se))
  cat(sprintf("  95%% CI: [%.2f, %.2f]\n", gs_spatial$ci_lower, gs_spatial$ci_upper))
  cat(sprintf("  p-value: %.4f\n", gs_spatial$p_value))
  
  if (gs_spatial$estimate > 0 && gs_spatial$p_value < 0.1) {
    cat("  >> DISPLACEMENT detected (crime increased in buffer)\n")
  } else if (gs_spatial$estimate < 0 && gs_spatial$p_value < 0.1) {
    cat("  >> DIFFUSION OF BENEFITS (crime decreased in buffer)\n")
  } else {
    cat("  >> No significant spatial spillover\n")
  }
  
  # Gsynth plots
  p_gs_spatial_cf <- plot(gs_spatial$object, type = "counterfactual", main = "Gsynth: Spatial Spillover vs Counterfactual")
  print(p_gs_spatial_cf)
  ggsave(here("output", "causal_final", "gs_spatial_counterfactual.png"), p_gs_spatial_cf,
         width = 10, height = 6, dpi = 300)
  
  p_gs_spatial_gap <- plot(gs_spatial$object, type = "gap", main = "Gsynth: Effect Over Time (Spatial)")
  print(p_gs_spatial_gap)
  ggsave(here("output", "causal_final", "gs_spatial_gap.png"), p_gs_spatial_gap,
         width = 10, height = 6, dpi = 300)
  
  cat("  Gsynth plots saved.\n")
}

# --- SDID ---
cat("\nFitting SDID for spatial spillover...\n")
spatial_matrix <- prepare_sdid_matrix(spatial_panel, donor_panel, good_donors)
sdid_spatial <- run_sdid(spatial_matrix$Y, spatial_matrix$N0, spatial_matrix$T0, "Spatial - SDID")

if (sdid_spatial$success) {
  cat(sprintf("  ATT: %.2f (SE: %.2f)\n", sdid_spatial$estimate, sdid_spatial$se))
  cat(sprintf("  95%% CI: [%.2f, %.2f]\n", sdid_spatial$ci_lower, sdid_spatial$ci_upper))
  cat(sprintf("  p-value: %.4f\n", sdid_spatial$p_value))
  
  p_sdid_spatial <- synthdid_plot(sdid_spatial$object, se.method = "placebo")
  print(p_sdid_spatial)
  ggsave(here("output", "causal_final", "sdid_spatial.png"), p_sdid_spatial,
         width = 10, height = 6, dpi = 300)
  
  cat("  SDID plot saved.\n")
}

#==============================================================================
# 7. THEMATIC SPILLOVER ANALYSIS (PROSTITUTION HOT SPOTS)
#==============================================================================

cat("\n")
cat("================================================================\n")
cat("THEMATIC SPILLOVER: PROSTITUTION HOT SPOTS\n")
cat("================================================================\n")

thematic_main_panel <- donor_panel %>%
  filter(nta_id %in% good_donors) %>%
  bind_rows(thematic_panel)

# --- Gsynth ---
cat("\nFitting gsynth for thematic spillover...\n")
gs_thematic <- run_gsynth(thematic_main_panel, "THEMATIC_SPILLOVER", "Thematic - Gsynth")

if (gs_thematic$success) {
  cat(sprintf("  ATT: %.2f (SE: %.2f)\n", gs_thematic$estimate, gs_thematic$se))
  cat(sprintf("  95%% CI: [%.2f, %.2f]\n", gs_thematic$ci_lower, gs_thematic$ci_upper))
  cat(sprintf("  p-value: %.4f\n", gs_thematic$p_value))
  
  if (gs_thematic$estimate > 0 && gs_thematic$p_value < 0.1) {
    cat("  >> THEMATIC DISPLACEMENT (crime shifted to other prostitution areas)\n")
  } else if (gs_thematic$estimate < 0 && gs_thematic$p_value < 0.1) {
    cat("  >> Crime also DECREASED in other prostitution areas\n")
  } else {
    cat("  >> No significant thematic spillover\n")
  }
  
  # Gsynth plots
  p_gs_thematic_cf <- plot(gs_thematic$object, type = "counterfactual", main = "Gsynth: Thematic Spillover vs Counterfactual")
  print(p_gs_thematic_cf)
  ggsave(here("output", "causal_final", "gs_thematic_counterfactual.png"), p_gs_thematic_cf,
         width = 10, height = 6, dpi = 300)
  
  p_gs_thematic_gap <- plot(gs_thematic$object, type = "gap", main = "Gsynth: Effect Over Time (Thematic)")
  print(p_gs_thematic_gap)
  ggsave(here("output", "causal_final", "gs_thematic_gap.png"), p_gs_thematic_gap,
         width = 10, height = 6, dpi = 300)
  
  cat("  Gsynth plots saved.\n")
}

# --- SDID ---
cat("\nFitting SDID for thematic spillover...\n")
thematic_matrix <- prepare_sdid_matrix(thematic_panel, donor_panel, good_donors)
sdid_thematic <- run_sdid(thematic_matrix$Y, thematic_matrix$N0, thematic_matrix$T0, "Thematic - SDID")

if (sdid_thematic$success) {
  cat(sprintf("  ATT: %.2f (SE: %.2f)\n", sdid_thematic$estimate, sdid_thematic$se))
  cat(sprintf("  95%% CI: [%.2f, %.2f]\n", sdid_thematic$ci_lower, sdid_thematic$ci_upper))
  cat(sprintf("  p-value: %.4f\n", sdid_thematic$p_value))
  
  p_sdid_thematic <- synthdid_plot(sdid_thematic$object, se.method = "placebo")
  print(p_sdid_thematic)
  ggsave(here("output", "causal_final", "sdid_thematic.png"), p_sdid_thematic,
         width = 10, height = 6, dpi = 300)
  
  cat("  SDID plot saved.\n")
}
#==============================================================================
# 8. RESULTS TABLE
#==============================================================================

cat("\n")
cat("================================================================\n")
cat("RESULTS SUMMARY\n")
cat("================================================================\n")

results_list <- list(
  gs_zone = gs_zone,
  sdid_zone = sdid_zone,
  gs_spatial = gs_spatial,
  sdid_spatial = sdid_spatial,
  gs_thematic = gs_thematic,
  sdid_thematic = sdid_thematic
)

results_table <- tibble(
  Analysis = c("Zone", "Zone", "Spatial Spillover", "Spatial Spillover", 
               "Thematic Spillover", "Thematic Spillover"),
  Method = rep(c("Gsynth", "SDID"), 3),
  Estimate = c(gs_zone$estimate, sdid_zone$estimate,
               gs_spatial$estimate, sdid_spatial$estimate,
               gs_thematic$estimate, sdid_thematic$estimate),
  SE = c(gs_zone$se, sdid_zone$se,
         gs_spatial$se, sdid_spatial$se,
         gs_thematic$se, sdid_thematic$se),
  CI_Lower = c(gs_zone$ci_lower, sdid_zone$ci_lower,
               gs_spatial$ci_lower, sdid_spatial$ci_lower,
               gs_thematic$ci_lower, sdid_thematic$ci_lower),
  CI_Upper = c(gs_zone$ci_upper, sdid_zone$ci_upper,
               gs_spatial$ci_upper, sdid_spatial$ci_upper,
               gs_thematic$ci_upper, sdid_thematic$ci_upper),
  P_Value = c(gs_zone$p_value, sdid_zone$p_value,
              gs_spatial$p_value, sdid_spatial$p_value,
              gs_thematic$p_value, sdid_thematic$p_value)
) %>%
  mutate(
    Significant = case_when(
      P_Value < 0.01 ~ "***",
      P_Value < 0.05 ~ "**",
      P_Value < 0.1 ~ "*",
      TRUE ~ ""
    )
  )

print(results_table)

write_csv(results_table, here("output", "causal_final", "results_summary.csv"))

#==============================================================================
# 9. VISUALIZATIONS
#==============================================================================

cat("\n=== GENERATING PLOTS ===\n")

# --- 9a. Time series comparison ---
ts_data <- bind_rows(
  zone_panel %>% mutate(Unit = "Treatment Zone"),
  spatial_panel %>% mutate(Unit = "Spatial Spillover"),
  thematic_panel %>% mutate(Unit = "Thematic Spillover")
)

p_timeseries <- ggplot(ts_data, aes(x = month, y = crime_count, color = Unit)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  geom_vline(xintercept = treatment_month, linetype = "dashed") +
  annotate("rect", xmin = op_start, xmax = op_end, ymin = -Inf, ymax = Inf,
           alpha = 0.15, fill = "red") +
  annotate("text", x = op_start + 45, y = max(ts_data$crime_count) * 0.95,
           label = "Operation\nPeriod", size = 3, hjust = 0.5) +
  scale_color_manual(values = c(
    "Treatment Zone" = "darkred",
    "Spatial Spillover" = "steelblue",
    "Thematic Spillover" = "forestgreen"
  )) +
  labs(
    title = "Monthly Street Violence: Treatment Zone vs. Spillover Areas",
    x = "Month", y = "Crime Count", color = ""
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p_timeseries)
ggsave(here("output", "causal_final", "01_timeseries.png"), p_timeseries,
       width = 12, height = 6, dpi = 300)

# --- 9b. Forest plot ---
p_forest <- results_table %>%
  mutate(Label = paste(Analysis, "-", Method)) %>%
  ggplot(aes(x = Estimate, y = reorder(Label, Estimate))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0.2) +
  geom_point(aes(color = Analysis), size = 3) +
  scale_color_manual(values = c(
    "Zone" = "darkred",
    "Spatial Spillover" = "steelblue",
    "Thematic Spillover" = "forestgreen"
  )) +
  labs(
    title = "Treatment Effect Estimates",
    subtitle = "Point estimates with 95% CI. Negative = crime reduction.",
    x = "Effect (change in monthly crimes)", y = ""
  ) +
  theme_minimal() +
  theme(legend.position = "none")

print(p_forest)
ggsave(here("output", "causal_final", "02_forest_plot.png"), p_forest,
       width = 10, height = 6, dpi = 300)

# --- 9c. SDID plot for zone ---
if (sdid_zone$success) {
  p_sdid_zone <- synthdid_plot(sdid_zone$object, se.method = "placebo")
  print(p_sdid_zone)
  ggsave(here("output", "causal_final", "03_sdid_zone.png"), p_sdid_zone,
         width = 10, height = 6, dpi = 300)
}

# --- 9d. Gsynth counterfactual ---
if (gs_zone$success) {
  png(here("output", "causal_final", "04_gsynth_zone.png"),
      width = 10, height = 6, units = "in", res = 300)
  plot(gs_zone$object, type = "counterfactual", main = "Gsynth: Zone vs Counterfactual")
  dev.off()
  
  # Gap plot
  png(here("output", "causal_final", "05_gsynth_gap.png"),
      width = 10, height = 6, units = "in", res = 300)
  plot(gs_zone$object, type = "gap", main = "Gsynth: Treatment Effect Over Time")
  dev.off()
}

# --- 9e. ATT by period (pre-treatment fit diagnostic) ---
if (gs_zone$success) {
  gs <- gs_zone$object
  
  att_df <- tibble(
    relative_period = gs$time,
    ATT = gs$att,
    month = sort(unique(main_panel$month))[which(sort(unique(main_panel$month)) == treatment_month) + gs$time]
  ) %>%
    mutate(period = if_else(relative_period < 0, "Pre-treatment", "Post-treatment"))
  
  if (!is.null(gs$est.att)) {
    att_df <- att_df %>%
      mutate(
        SE = gs$est.att[, "S.E."],
        CI_lower = gs$est.att[, "CI.lower"],
        CI_upper = gs$est.att[, "CI.upper"]
      )
  }
  
  p_att <- ggplot(att_df, aes(x = month, y = ATT)) +
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
      title = "Gsynth: Period-by-Period Treatment Effect",
      subtitle = "Pre-treatment ATT should be ~0 if model fits well",
      x = "Month", y = "ATT", color = ""
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  print(p_att)
  ggsave(here("output", "causal_final", "06_att_by_period.png"), p_att,
         width = 10, height = 6, dpi = 300)
  
  # Pre-treatment fit diagnostics
  pre_att <- att_df %>% filter(period == "Pre-treatment")
  cat(sprintf("\nPre-treatment fit:\n"))
  cat(sprintf("  Mean ATT: %.2f (should be ~0)\n", mean(pre_att$ATT)))
  cat(sprintf("  RMSE: %.2f\n", sqrt(mean(pre_att$ATT^2))))
}

#==============================================================================
# 10. INTERPRETATION SUMMARY
#==============================================================================

cat("\n")
cat("================================================================\n")
cat("INTERPRETATION SUMMARY\n")
cat("================================================================\n")

cat("\n--- MAIN EFFECT (ZONE) ---\n")
if (gs_zone$success && sdid_zone$success) {
  avg_effect <- mean(c(gs_zone$estimate, sdid_zone$estimate))
  zone_mean <- mean(zone_panel$crime_count)
  pct_reduction <- abs(avg_effect) / zone_mean * 100
  
  cat(sprintf("Gsynth estimate: %.1f crimes/month (p=%.4f)\n", gs_zone$estimate, gs_zone$p_value))
  cat(sprintf("SDID estimate: %.1f crimes/month (p=%.4f)\n", sdid_zone$estimate, sdid_zone$p_value))
  cat(sprintf("Average effect: %.1f crimes/month (~%.0f%% reduction)\n", avg_effect, pct_reduction))
  
  if (gs_zone$p_value < 0.05 || sdid_zone$p_value < 0.05) {
    cat(">> CONCLUSION: Statistically significant crime reduction in zone\n")
  }
}

cat("\n--- SPATIAL SPILLOVER ---\n")
if (gs_spatial$success && sdid_spatial$success) {
  cat(sprintf("Gsynth estimate: %.1f (p=%.4f)\n", gs_spatial$estimate, gs_spatial$p_value))
  cat(sprintf("SDID estimate: %.1f (p=%.4f)\n", sdid_spatial$estimate, sdid_spatial$p_value))
  
  if (gs_spatial$p_value > 0.1 && sdid_spatial$p_value > 0.1) {
    cat(">> CONCLUSION: No significant spatial displacement detected\n")
  }
}

cat("\n--- THEMATIC SPILLOVER ---\n")
if (gs_thematic$success && sdid_thematic$success) {
  cat(sprintf("Gsynth estimate: %.1f (p=%.4f)\n", gs_thematic$estimate, gs_thematic$p_value))
  cat(sprintf("SDID estimate: %.1f (p=%.4f)\n", sdid_thematic$estimate, sdid_thematic$p_value))
  
  if (gs_thematic$p_value > 0.1 && sdid_thematic$p_value > 0.1) {
    cat(">> CONCLUSION: No significant displacement to other prostitution areas\n")
  }
}

#==============================================================================
# 11. SAVE RESULTS
#==============================================================================

saveRDS(results_list, here("output", "causal_final", "all_results.rds"))
saveRDS(list(
  zone = zone_panel,
  spatial = spatial_panel,
  thematic = thematic_panel,
  donors = donor_panel %>% filter(nta_id %in% good_donors)
), here("output", "causal_final", "panels.rds"))

plan(sequential)

cat("\n================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("================================================================\n")

list.files(here("output", "causal_final")) %>%
  paste(" -", .) %>%
  cat(sep = "\n")

#==============================================================================
# 12. SUGGESTED ROBUSTNESS CHECKS
#==============================================================================

cat("\n")
cat("================================================================\n")
cat("SUGGESTED ROBUSTNESS CHECKS (run separately)\n")
cat("================================================================\n")

cat("
1. PRE-SPIKE MATCHING
   - Re-estimate gsynth/SDID matching only on Jan 2022 - May 2024
   - Avoids the summer 2024 spike influencing weights
   
2. PLACEBO TIMING TEST
   - Assign fake treatment at Oct 2023 (1 year earlier)
   - Should find NO effect if model is well-specified
   
3. LEAVE-ONE-OUT DONORS
   - Remove top 5 most-weighted donors, re-estimate
   - Check sensitivity to influential donors
   
4. ALTERNATIVE DONOR POOLS
   - Queens-only donors (same borough)
   - High-crime NTAs only (better scale match)
   
5. DIFFERENT OUTCOME MEASURES
   - Specific crime types (robbery, assault, etc.)
   - Log transformation robustness
   
6. EXTENDED POST-PERIOD
   - As more 2025 data becomes available
   - Check if effect persists or fades

7. INTERRUPTED TIME SERIES
   - Complement with ITS on zone alone (no donors)
   - Provides different identifying assumptions
")




#==============================================================================
# 13. INTERACTIVE MAP: INSPECT SYNTHETIC CONTROL WEIGHTS
#==============================================================================

library(leaflet)
library(htmltools)

cat("\n=== CREATING INTERACTIVE MAP ===\n")

# Transform to WGS84 for leaflet
nynta_wgs <- nynta %>% st_transform(4326)
nyc_bgs_wgs <- nyc_bgs %>% st_transform(4326)
zone_wgs <- rosie_buffer %>% st_transform(4326)
spatial_buffer_wgs <- spatial_buffer %>% st_transform(4326)
thematic_bgs_wgs <- nyc_bgs %>% 
  filter(geoid %in% thematic_spillover_bgs) %>% 
  st_transform(4326)

# --- Extract SDID weights ---
if (sdid_zone$success) {
  omega <- attr(sdid_zone$object, "weights")$omega
  donor_names <- rownames(main_matrix$Y)[1:main_matrix$N0]
  
  sdid_weights <- tibble(
    nta_id = donor_names,
    sdid_weight = as.vector(omega)
  ) %>%
    arrange(desc(sdid_weight))
  
  cat("Top 10 SDID donor weights:\n")
  print(head(sdid_weights, 10))
}

# --- Extract Gsynth weights (if available) ---
# --- Extract Gsynth weights (if available) ---
gsynth_weights <- tibble(nta_id = character(), gsynth_weight = numeric())

if (gs_zone$success) {
  gs <- gs_zone$object
  
  # Check what's available
  cat("Gsynth weight objects available:\n")
  cat("  wgt.implied:", !is.null(gs$wgt.implied), "\n")
  
  if (!is.null(gs$wgt.implied)) {
    # Check structure
    cat("  Class:", class(gs$wgt.implied), "\n")
    cat("  Length:", length(gs$wgt.implied), "\n")
    
    # Get names - might be names() or rownames() depending on structure
    weight_names <- if (!is.null(rownames(gs$wgt.implied))) {
      rownames(gs$wgt.implied)
    } else if (!is.null(names(gs$wgt.implied))) {
      names(gs$wgt.implied)
    } else {
      NULL
    }
    
    if (!is.null(weight_names)) {
      gsynth_weights <- tibble(
        nta_id = weight_names,
        gsynth_weight = as.vector(gs$wgt.implied)
      )
      
      gsynth_weights <- gsynth_weights %>%
        filter(nta_id != "ZONE") %>%
        arrange(desc(gsynth_weight))
      
      cat("\nTop 10 Gsynth implied weights:\n")
      print(head(gsynth_weights, 10))
    } else {
      cat("  Could not extract weight names from gsynth object\n")
    }
  }
}

# --- Combine weights with NTA geometry ---
nta_map_data <- nynta_wgs %>%
  left_join(sdid_weights, by = "nta_id") %>%
  left_join(gsynth_weights, by = "nta_id") %>%
  mutate(
    # Classify each NTA (now just donors or not)
    unit_type = case_when(
      nta_id %in% good_donors ~ "Donor",
      TRUE ~ "Excluded"
    ),
    # Weight category for donors
    weight_cat = case_when(
      unit_type != "Donor" ~ "Excluded",
      sdid_weight >= 0.05 ~ "High Weight (≥5%)",
      sdid_weight >= 0.01 ~ "Medium Weight (1-5%)",
      sdid_weight > 0 ~ "Low Weight (<1%)",
      TRUE ~ "Zero Weight"
    ),
    # Create popup label
    popup_label = paste0(
      "<b>", nta_name, "</b><br>",
      "NTA: ", nta_id, "<br>",
      "Borough: ", boro_name, "<br>",
      if_else(!is.na(sdid_weight) & sdid_weight > 0, 
              paste0("SDID Weight: ", scales::percent(sdid_weight, accuracy = 0.1), "<br>"),
              ""),
      if_else(!is.na(gsynth_weight) & gsynth_weight > 0,
              paste0("Gsynth Weight: ", scales::percent(gsynth_weight, accuracy = 0.1)),
              "")
    )
  )

# --- Add prostitution arrest info to thematic BGs ---
thematic_bgs_map <- thematic_bgs_wgs %>%
  left_join(pros_by_bg, by = "geoid") %>%
  mutate(
    popup_label = paste0(
      "<b>Thematic Spillover BG</b><br>",
      "GEOID: ", geoid, "<br>",
      "Prostitution Arrests: ", pros_arrests
    )
  )

# --- Color palettes ---
pal_donors <- colorFactor(
  palette = c("orange", "gold", "lightyellow", "gray90", "gray95"),
  levels = c("High Weight (≥5%)", "Medium Weight (1-5%)", "Low Weight (<1%)", "Zero Weight", "Excluded")
)

# --- Create leaflet map ---
map <- leaflet() %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  
  # Add donor NTAs colored by weight
  addPolygons(
    data = nta_map_data %>% filter(unit_type == "Donor"),
    fillColor = ~pal_donors(weight_cat),
    fillOpacity = 0.5,
    color = "white",
    weight = 1,
    popup = ~popup_label,
    label = ~nta_name,
    group = "Donor NTAs"
  ) %>%
  
  # Add treatment zone
  addPolygons(
    data = zone_wgs,
    fillColor = "darkred",
    fillOpacity = 0.5,
    color = "darkred",
    weight = 3,
    label = "Treatment Zone (rosie_buffer)",
    popup = paste0(
      "<b>Treatment Zone</b><br>",
      "Mean monthly crime: ", round(mean(zone_panel$crime_count), 1)
    ),
    group = "Treatment Zone"
  ) %>%
  
  # Add spatial spillover buffer (donut)
  addPolygons(
    data = spatial_buffer_wgs,
    fillColor = "steelblue",
    fillOpacity = 0.4,
    color = "steelblue",
    weight = 2,
    label = "Spatial Spillover Buffer (0.5 mi)",
    popup = paste0(
      "<b>Spatial Spillover Buffer</b><br>",
      "0.5 mile donut around zone<br>",
      "Mean monthly crime: ", round(mean(spatial_panel$crime_count), 1)
    ),
    group = "Spatial Spillover"
  ) %>%
  
  # Add thematic spillover BGs
  addPolygons(
    data = thematic_bgs_map,
    fillColor = "forestgreen",
    fillOpacity = 0.5,
    color = "forestgreen",
    weight = 2,
    popup = ~popup_label,
    label = ~paste("Prostitution BG:", geoid),
    group = "Thematic Spillover (Prostitution BGs)"
  ) %>%
  
  # Add legend
  addLegend(
    position = "bottomright",
    colors = c("darkred", "steelblue", "forestgreen", "orange", "gold", "lightyellow", "gray90"),
    labels = c("Treatment Zone", "Spatial Spillover Buffer", "Thematic Spillover BGs",
               "Donor: High Weight (≥5%)", "Donor: Medium (1-5%)", "Donor: Low (<1%)", "Donor: Zero Weight"),
    title = "Map Legend",
    opacity = 0.8
  ) %>%
  
  # Layer controls
  addLayersControl(
    overlayGroups = c("Donor NTAs", "Treatment Zone", "Spatial Spillover", "Thematic Spillover (Prostitution BGs)"),
    options = layersControlOptions(collapsed = FALSE)
  ) %>%
  
  # Set view to Queens/Roosevelt Ave area
  
  setView(lng = -73.85, lat = 40.75, zoom = 12)

# Display map
map

# Save as HTML
htmlwidgets::saveWidget(map, here("output", "causal_final", "07_synthetic_control_map.html"),
                        selfcontained = TRUE)

cat("Interactive map saved: 07_synthetic_control_map.html\n")

# --- Table of top donors ---
cat("\n=== TOP DONORS BY WEIGHT ===\n")

top_donors_table <- nta_map_data %>%
  st_drop_geometry() %>%
  filter(unit_type == "Donor", sdid_weight > 0) %>%
  select(nta_id, nta_name, boro_name, sdid_weight, gsynth_weight) %>%
  arrange(desc(sdid_weight)) %>%
  head(20)

print(top_donors_table)

write_csv(top_donors_table, here("output", "causal_final", "top_donors.csv"))

# --- Summary of geographic distribution ---
cat("\n=== DONOR WEIGHT BY BOROUGH ===\n")

boro_weights <- nta_map_data %>%
  st_drop_geometry() %>%
  filter(unit_type == "Donor") %>%
  group_by(boro_name) %>%
  summarize(
    n_donors = n(),
    n_nonzero = sum(sdid_weight > 0, na.rm = TRUE),
    total_weight = sum(sdid_weight, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(total_weight))

print(boro_weights)

# --- Thematic spillover BGs summary ---
cat("\n=== THEMATIC SPILLOVER BGs ===\n")

# Use countyfp directly (005=Bronx, 047=Brooklyn, 061=Manhattan, 081=Queens, 085=Staten Island)
thematic_summary <- pros_by_bg %>%
  filter(geoid %in% thematic_spillover_bgs) %>%
  left_join(
    nyc_bgs %>% st_drop_geometry() %>% select(geoid, countyfp),
    by = "geoid"
  ) %>%
  mutate(
    boro_name = case_when(
      countyfp == "005" ~ "Bronx",
      countyfp == "047" ~ "Brooklyn",
      countyfp == "061" ~ "Manhattan",
      countyfp == "081" ~ "Queens",
      countyfp == "085" ~ "Staten Island",
      TRUE ~ "Unknown"
    )
  )

cat(sprintf("Total thematic spillover BGs: %d\n", length(thematic_spillover_bgs)))
cat(sprintf("Total prostitution arrests in these BGs: %d\n", sum(thematic_summary$pros_arrests)))

cat("\nBy borough:\n")
thematic_summary %>%
  group_by(boro_name) %>%
  summarize(
    n_bgs = n(),
    total_arrests = sum(pros_arrests),
    .groups = "drop"
  ) %>%
  arrange(desc(total_arrests)) %>%
  print()

