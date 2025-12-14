#==============================================================================
# Operation Restore Roosevelt - Supplementary EDA
# 
# Additional analyses for methods planning:
# - Pre-intervention spike characterization
# - Donor pool identification
# - Treatment intensity measures
# - Vending enforcement (separate from prostitution)
#==============================================================================

library(here)
library(tidyverse)
library(janitor)
library(lubridate)
library(sf)
library(scales)
library(patchwork)

# Load EDA objects if not in environment
if (!exists("zone")) {
  eda_objects <- readRDS(here("output", "eda_objects.rds"))
  list2env(eda_objects, envir = .GlobalEnv)
}

#==============================================================================
# 1. SPIKE CHARACTERIZATION
#==============================================================================

cat("\n=== PRE-INTERVENTION SPIKE ANALYSIS ===\n")

# Get longer historical view
historical_weekly <- street_violent_zone %>%
  st_drop_geometry() %>%
  filter(date >= as.Date("2019-01-01")) %>%
  mutate(week = floor_date(date, "week")) %>%
  count(week)

# Calculate rolling statistics
historical_weekly <- historical_weekly %>%
  arrange(week) %>%
  mutate(
    rolling_mean_8wk = zoo::rollmean(n, k = 8, fill = NA, align = "right"),
    rolling_sd_8wk = zoo::rollapply(n, width = 8, FUN = sd, fill = NA, align = "right"),
    z_score = (n - rolling_mean_8wk) / rolling_sd_8wk
  )

# Identify anomalous weeks (>2 SD above rolling mean)
anomalous_weeks <- historical_weekly %>%
  filter(z_score > 2, !is.na(z_score))

cat(sprintf("Weeks with z-score > 2: %d\n", nrow(anomalous_weeks)))
cat("Most recent anomalous weeks:\n")
print(tail(anomalous_weeks %>% select(week, n, rolling_mean_8wk, z_score), 10))

# Plot with anomaly highlighting
p_anomaly <- ggplot(historical_weekly, aes(x = week, y = n)) +
  geom_line(alpha = 0.5) +
  geom_line(aes(y = rolling_mean_8wk), color = "blue", linewidth = 0.8) +
  geom_ribbon(aes(ymin = rolling_mean_8wk - 2*rolling_sd_8wk,
                  ymax = rolling_mean_8wk + 2*rolling_sd_8wk),
              alpha = 0.2, fill = "blue") +
  geom_point(data = anomalous_weeks, aes(x = week, y = n), 
             color = "red", size = 2) +
  geom_vline(xintercept = op_start, linetype = "dashed", color = "black") +
  annotate("rect", xmin = op_start, xmax = op_end, ymin = -Inf, ymax = Inf,
           alpha = 0.1, fill = "red") +
  labs(
    title = "Weekly Street Violence with Anomaly Detection",
    subtitle = "Red points = weeks >2 SD above 8-week rolling mean",
    x = "Week", y = "Incidents"
  ) +
  theme_minimal()

p_anomaly

ggsave(here("output", "eda_plots", "20_anomaly_detection.png"), p_anomaly, 
       width = 12, height = 6, dpi = 300)

# When did the spike start?
# Look for sustained elevation in the months before operation
pre_op_6mo <- historical_weekly %>%
  filter(week >= op_start - months(6), week < op_start)

spike_start <- pre_op_6mo %>%
  filter(z_score > 1.5) %>%
  slice_min(week, n = 1)

if (nrow(spike_start) > 0) {
  cat(sprintf("\nSpike appears to start around: %s\n", spike_start$week))
  cat(sprintf("  Count that week: %d (z-score: %.2f)\n", spike_start$n, spike_start$z_score))
}

#==============================================================================
# 2. HISTORICAL SPIKE PATTERNS (FOR COUNTERFACTUAL MODELING)
#==============================================================================

# Look for historical spikes and recovery patterns
# This helps us understand "what does natural recovery from a spike look like?"

historical_spikes <- historical_weekly %>%
  filter(z_score > 2, !is.na(z_score)) %>%
  filter(week < op_start - months(3))  # Exclude recent spike

cat("\n=== HISTORICAL SPIKE RECOVERY PATTERNS ===\n")

if (nrow(historical_spikes) > 0) {
  # For each historical spike, look at the 12 weeks after
  recovery_patterns <- map_dfr(historical_spikes$week, function(spike_week) {
    historical_weekly %>%
      filter(week >= spike_week, week < spike_week + weeks(12)) %>%
      mutate(
        spike_id = as.character(spike_week),
        weeks_post = as.integer(difftime(week, spike_week, units = "weeks"))
      )
  })
  
  # Average recovery pattern
  avg_recovery <- recovery_patterns %>%
    group_by(weeks_post) %>%
    summarize(
      mean_count = mean(n),
      sd_count = sd(n),
      n_spikes = n(),
      .groups = "drop"
    )
  
  cat("Average recovery pattern after historical spikes:\n")
  print(avg_recovery)
  
  p_recovery <- ggplot(recovery_patterns, aes(x = weeks_post, y = n, group = spike_id)) +
    geom_line(alpha = 0.3) +
    geom_line(data = avg_recovery, aes(x = weeks_post, y = mean_count, group = 1),
              color = "red", linewidth = 1.5) +
    labs(
      title = "Recovery Patterns After Historical Spikes",
      subtitle = "Red line = average across all historical spikes",
      x = "Weeks After Spike", y = "Weekly Incidents"
    ) +
    theme_minimal()
  
  ggsave(here("output", "eda_plots", "21_spike_recovery_patterns.png"), p_recovery, 
         width = 10, height = 6, dpi = 300)
}

#==============================================================================
# 3. DONOR POOL IDENTIFICATION
#==============================================================================

cat("\n=== DONOR POOL ANALYSIS ===\n")

# We need areas that:
# 1. Have similar baseline crime levels
# 2. Did NOT experience the pre-intervention spike
# 3. Are not adjacent to the zone (no displacement contamination)

# Get zone's average pre-spike level
zone_baseline <- historical_weekly %>%
  filter(week >= as.Date("2023-01-01"), week < op_start - months(3)) %>%
  summarize(mean_weekly = mean(n), sd_weekly = sd(n))

cat(sprintf("Zone baseline (2023 to 3mo pre-op): mean=%.2f, sd=%.2f\n",
            zone_baseline$mean_weekly, zone_baseline$sd_weekly))

# Get all census tracts' weekly crime series
ct_weekly <- street_violent_crime %>%
  filter(date >= as.Date("2023-01-01")) %>%
  st_join(nyct %>% dplyr::select(boro_ct2020)) %>%
  st_drop_geometry() %>%
  filter(!is.na(boro_ct2020)) %>%
  mutate(week = floor_date(date, "week")) %>%
  count(boro_ct2020, week)

# Fill in zero-count weeks
all_weeks <- tibble(week = seq(as.Date("2023-01-01"), max(ct_weekly$week), by = "week"))
all_cts <- unique(ct_weekly$boro_ct2020)

ct_weekly_complete <- expand_grid(boro_ct2020 = all_cts, week = all_weeks$week) %>%
  left_join(ct_weekly, by = c("boro_ct2020", "week")) %>%
  replace_na(list(n = 0))

# Calculate tract-level stats
ct_stats <- ct_weekly_complete %>%
  filter(week < op_start - months(3)) %>%  # Pre-spike period
  group_by(boro_ct2020) %>%
  summarize(
    mean_weekly = mean(n),
    sd_weekly = sd(n),
    total = sum(n),
    .groups = "drop"
  )

# Find tracts with similar baseline (within 1 SD of zone)
similar_tracts <- ct_stats %>%
  filter(
    mean_weekly >= zone_baseline$mean_weekly - zone_baseline$sd_weekly,
    mean_weekly <= zone_baseline$mean_weekly + zone_baseline$sd_weekly,
    total >= 20  # Minimum data requirement
  ) %>%
  filter(!boro_ct2020 %in% zone_cts)  # Exclude zone tracts

cat(sprintf("\nTracts with similar baseline: %d\n", nrow(similar_tracts)))

# Check which of these also spiked (to exclude them)
# Look at recent 3 months before operation
spike_period_cts <- ct_weekly_complete %>%
  filter(week >= op_start - months(3), week < op_start) %>%
  group_by(boro_ct2020) %>%
  summarize(
    spike_mean = mean(n),
    spike_max = max(n),
    .groups = "drop"
  )

# Join and identify non-spiking similar tracts
potential_donors <- similar_tracts %>%
  left_join(spike_period_cts, by = "boro_ct2020") %>%
  mutate(spike_ratio = spike_mean / mean_weekly) %>%
  filter(spike_ratio < 1.5)  # Did not spike more than 50% above baseline

cat(sprintf("Potential donor tracts (similar baseline, no spike): %d\n", nrow(potential_donors)))

# Exclude tracts adjacent to zone
zone_union <- zone %>% st_union() %>% st_buffer(500)  # 500ft buffer

adjacent_cts <- nyct %>%
  st_filter(zone_union, .predicate = st_intersects) %>%
  pull(boro_ct2020)

final_donors <- potential_donors %>%
  filter(!boro_ct2020 %in% adjacent_cts)

cat(sprintf("Final donor pool (excluding adjacent): %d tracts\n", nrow(final_donors)))

# Save donor pool
write_csv(final_donors, here("output", "eda_plots", "potential_donor_tracts.csv"))

# Visualize donor pool comparison
donor_weekly <- ct_weekly_complete %>%
  filter(boro_ct2020 %in% final_donors$boro_ct2020) %>%
  group_by(week) %>%
  summarize(
    donor_mean = mean(n),
    donor_se = sd(n) / sqrt(n()),
    .groups = "drop"
  )

zone_weekly_agg <- historical_weekly %>%
  select(week, n) %>%
  rename(zone_count = n)

comparison_df <- donor_weekly %>%
  left_join(zone_weekly_agg, by = "week") %>%
  filter(!is.na(zone_count))

p_donor_comparison <- ggplot(comparison_df, aes(x = week)) +
  geom_ribbon(aes(ymin = donor_mean - 1.96*donor_se, 
                  ymax = donor_mean + 1.96*donor_se),
              alpha = 0.3, fill = "blue") +
  geom_line(aes(y = donor_mean, color = "Donor Pool Mean"), linewidth = 0.8) +
  geom_line(aes(y = zone_count, color = "Zone"), linewidth = 0.8) +
  geom_vline(xintercept = op_start, linetype = "dashed") +
  annotate("rect", xmin = op_start, xmax = op_end, ymin = -Inf, ymax = Inf,
           alpha = 0.1, fill = "red") +
  labs(
    title = "Zone vs. Potential Donor Pool",
    subtitle = "Blue ribbon = 95% CI around donor mean",
    x = "Week", y = "Weekly Street Violence", color = ""
  ) +
  scale_color_manual(values = c("Zone" = "red", "Donor Pool Mean" = "blue")) +
  theme_minimal() +
  theme(legend.position = "bottom")

p_donor_comparison

ggsave(here("output", "eda_plots", "22_donor_pool_comparison.png"), p_donor_comparison, 
       width = 12, height = 6, dpi = 300)

#==============================================================================
# 4. TREATMENT INTENSITY MEASURES
#==============================================================================

cat("\n=== TREATMENT INTENSITY ===\n")

# Combine all enforcement measures for a treatment intensity index
# Arrests + Summonses + Directed Patrols

weekly_enforcement <- bind_rows(
  arrests_zone %>%
    st_drop_geometry() %>%
    mutate(week = floor_date(date, "week")) %>%
    count(week) %>%
    mutate(measure = "Arrests"),
  
  summonses_zone %>%
    st_drop_geometry() %>%
    mutate(week = floor_date(date, "week")) %>%
    count(week) %>%
    mutate(measure = "Summonses"),
  
  directed_patrols_zone %>%
    st_drop_geometry() %>%
    mutate(week = floor_date(date, "week")) %>%
    count(week) %>%
    mutate(measure = "Directed Patrols")
) %>%
  filter(week >= as.Date("2022-01-01"))

# Pivot wider for combined index
weekly_enforcement_wide <- weekly_enforcement %>%
  pivot_wider(names_from = measure, values_from = n, values_fill = 0) %>%
  mutate(
    # Standardize each measure and create composite
    arrests_z = scale(Arrests)[,1],
    summonses_z = scale(Summonses)[,1],
    patrols_z = scale(`Directed Patrols`)[,1],
    composite_intensity = (arrests_z + summonses_z + patrols_z) / 3
  )

# Plot treatment intensity
p_intensity <- ggplot(weekly_enforcement_wide, aes(x = week)) +
  geom_line(aes(y = composite_intensity), color = "darkgreen", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = op_start, linetype = "dashed", color = "red") +
  annotate("rect", xmin = op_start, xmax = op_end, ymin = -Inf, ymax = Inf,
           alpha = 0.1, fill = "red") +
  labs(
    title = "Composite Enforcement Intensity in Zone",
    subtitle = "Standardized average of arrests, summonses, and directed patrols",
    x = "Week", y = "Z-Score (Composite)"
  ) +
  theme_minimal()

ggsave(here("output", "eda_plots", "23_treatment_intensity.png"), p_intensity, 
       width = 10, height = 6, dpi = 300)

# Enforcement during vs before operation
enforcement_comparison <- weekly_enforcement_wide %>%
  mutate(period = case_when(
    week >= op_start & week <= op_end ~ "During Operation",
    week >= op_start - weeks(13) & week < op_start ~ "13 Weeks Before",
    TRUE ~ "Historical"
  )) %>%
  filter(period != "Historical") %>%
  group_by(period) %>%
  summarize(
    across(c(Arrests, Summonses, `Directed Patrols`), 
           list(mean = mean, total = sum)),
    .groups = "drop"
  )

cat("\nEnforcement intensity comparison:\n")
print(enforcement_comparison)

write_csv(weekly_enforcement_wide, here("output", "eda_plots", "weekly_enforcement.csv"))

#==============================================================================
# 5. VENDING ENFORCEMENT (SEPARATE FROM PROSTITUTION)
#==============================================================================

cat("\n=== VENDING ENFORCEMENT ===\n")

# Look for vending-related enforcement in summonses
# Common codes: unlicensed general vendor, food vendor without permit, etc.

#vending_keywords <- c("VENDOR", "VENDING", "PEDDL", "FOOD CART", "SIDEWALK OBST")

vending_summonses <- all_summonses %>%
  filter(str_detect(str_to_upper(offense_description), 
                    paste(vending_tickets, collapse = "|")))

cat(sprintf("Total vending-related summonses (all time): %d\n", nrow(vending_summonses)))

# In zone
vending_zone <- vending_summonses %>%
  filter(date >= as.Date("2022-01-01")) %>%
  filter_to_zone()

cat(sprintf("Vending summonses in zone (2022+): %d\n", nrow(vending_zone)))

# Monthly trend
monthly_vending <- vending_zone %>%
  st_drop_geometry() %>%
  mutate(month = floor_date(date, "month")) %>%
  count(month)

if (nrow(monthly_vending) > 0) {
  p_vending <- ggplot(monthly_vending, aes(x = month, y = n)) +
    geom_col(fill = "orange", alpha = 0.8) +
    geom_vline(xintercept = op_start, linetype = "dashed", color = "black") +
    annotate("rect", xmin = op_start, xmax = op_end, ymin = -Inf, ymax = Inf,
             alpha = 0.1, fill = "red") +
    labs(
      title = "Monthly Vending Enforcement in Zone",
      subtitle = "2022 to Present",
      x = "Month", y = "Summonses"
    ) +
    theme_minimal()
  
  ggsave(here("output", "eda_plots", "24_vending_enforcement.png"), p_vending, 
         width = 10, height = 6, dpi = 300)
}
p_vending
#==============================================================================
# 6. CRIME TYPE BREAKDOWN DURING SPIKE
#==============================================================================

cat("\n=== CRIME COMPOSITION DURING SPIKE ===\n")

# What types of crime drove the spike?
spike_period <- street_violent_zone %>%
  st_drop_geometry() %>%
  filter(date >= op_start - months(3), date < op_start)

baseline_period <- street_violent_zone %>%
  st_drop_geometry() %>%
  filter(date >= op_start - months(12), date < op_start - months(3))

spike_composition <- spike_period %>%
  count(ofns_desc, name = "spike_count") %>%
  mutate(spike_pct = 100 * spike_count / sum(spike_count))

baseline_composition <- baseline_period %>%
  count(ofns_desc, name = "baseline_count") %>%
  mutate(baseline_pct = 100 * baseline_count / sum(baseline_count))

composition_comparison <- spike_composition %>%
  full_join(baseline_composition, by = "ofns_desc") %>%
  replace_na(list(spike_count = 0, baseline_count = 0, spike_pct = 0, baseline_pct = 0)) %>%
  mutate(pct_change = spike_pct - baseline_pct) %>%
  arrange(desc(abs(pct_change)))

cat("Crime type shifts during spike period:\n")
print(composition_comparison)

write_csv(composition_comparison, here("output", "eda_plots", "spike_crime_composition.csv"))

#==============================================================================
# 7. SUMMARY STATISTICS TABLE
#==============================================================================

# Create a comprehensive summary table for methods section

summary_stats <- tibble(
  Metric = c(
    "Zone area (sq miles)",
    "Zone as % of PCT 115",
    "Operation duration (days)",
    "Street violent crimes in zone (2019-present)",
    "Weekly street violence - Mean",
    "Weekly street violence - SD",
    "Weekly street violence - Max",
    "Arrests in zone (2022-present)",
    "Summonses in zone (2022-present)",
    "Directed patrols in zone (2022-present)",
    "Arrests during operation",
    "Arrests same period prior year",
    "Arrest change (%)",
    "Potential donor tracts identified",
    "Census tracts in zone",
    "Block groups in zone"
  ),
  Value = c(
    sprintf("%.4f", as.numeric(st_area(zone)) / 27878400),
    sprintf("%.2f%%", 100 * as.numeric(st_area(st_intersection(zone, pct_115))) / as.numeric(st_area(pct_115))),
    as.character(as.integer(op_end - op_start) + 1),
    as.character(nrow(street_violent_zone)),
    sprintf("%.2f", mean(historical_weekly$n, na.rm = TRUE)),
    sprintf("%.2f", sd(historical_weekly$n, na.rm = TRUE)),
    as.character(max(historical_weekly$n, na.rm = TRUE)),
    as.character(nrow(arrests_zone)),
    as.character(nrow(summonses_zone)),
    as.character(nrow(directed_patrols_zone)),
    as.character(nrow(arrests_zone %>% st_drop_geometry() %>% filter(date >= op_start, date <= op_end))),
    as.character(nrow(arrests_zone %>% st_drop_geometry() %>% filter(date >= comp_start, date <= comp_end))),
    sprintf("%.1f%%", 100 * (nrow(arrests_zone %>% st_drop_geometry() %>% filter(date >= op_start, date <= op_end)) - 
                               nrow(arrests_zone %>% st_drop_geometry() %>% filter(date >= comp_start, date <= comp_end))) /
              max(nrow(arrests_zone %>% st_drop_geometry() %>% filter(date >= comp_start, date <= comp_end)), 1)),
    as.character(nrow(final_donors)),
    as.character(length(zone_cts)),
    as.character(length(zone_bgs))
  )
)

write_csv(summary_stats, here("output", "eda_plots", "summary_statistics.csv"))

cat("\n")
cat("============================================================\n")
cat("SUPPLEMENTARY EDA COMPLETE\n")
cat("============================================================\n")
cat("\nNew outputs:\n")
list.files(here("output", "eda_plots"), pattern = "^2[0-4]") %>%
  paste(" -", .) %>%
  cat(sep = "\n")
