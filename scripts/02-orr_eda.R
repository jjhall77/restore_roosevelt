#==============================================================================
# Operation Restore Roosevelt - Exploratory Data Analysis
# 
# This script assumes all objects from the data loading script are already in 
# the environment. Run the data loading script first.
#==============================================================================

library(here)
library(tidyverse)
library(janitor)
library(lubridate)
library(sf)
library(scales)
library(patchwork)

# Create output directory for plots
dir.create(here("output", "eda_plots"), showWarnings = FALSE, recursive = TRUE)

#==============================================================================
# 1. DEFINE KEY PARAMETERS
#==============================================================================

# Operation window
op_start <- as.Date("2024-10-15")
op_end   <- op_start + 89

# Comparison period (same dates, prior year)
comp_start <- op_start - years(1)
comp_end   <- op_end - years(1)

# Pre-intervention period for arrest type comparison
pre_90_start <- op_start - 90
pre_90_end   <- op_start - 1

# Data start dates
data_start_2019 <- as.Date("2019-01-01")
data_start_2022 <- as.Date("2022-01-01")

# Make sure rosie_buffer is loaded
if (!exists("rosie_buffer")) {
  rosie_buffer <- readRDS(here("output", "rosie_buffer.rds"))
}

if (!exists("roosevelt_roads_limited")) {
  roosevelt_roads_limited <- readRDS(here("output", "roosevelt_roads_limited.rds"))
}

#==============================================================================
# 2. CREATE SPATIAL ZONES
#==============================================================================

# Zone = rosie_buffer (500 ft buffer around Roosevelt Ave segment)
zone <- rosie_buffer %>% st_as_sf() %>% st_set_crs(2263)

# 250 ft buffer around zone (displacement buffer) - excludes zone itself
buffer_250 <- zone %>%
  st_buffer(250) %>%
  st_difference(zone) %>%
  st_make_valid()

# Get precinct 115 polygon
pct_115 <- nypp %>% filter(precinct == 115)
pct_110 <- nypp %>% filter(precinct == 110)

# Verify zone is mostly in 115 (it spans 110 and 115)
zone_in_115 <- st_intersection(zone, pct_115) %>% st_area()
zone_total <- st_area(zone)

cat("\n=== ZONE STATISTICS ===\n")
cat(sprintf("Zone total area: %.2f sq ft (%.4f sq miles)\n", 
            as.numeric(zone_total), 
            as.numeric(zone_total) / 27878400))
cat(sprintf("Zone area in PCT 115: %.2f sq ft (%.1f%%)\n",
            as.numeric(zone_in_115),
            100 * as.numeric(zone_in_115) / as.numeric(zone_total)))

# Zone as % of precinct 115
pct_115_area <- st_area(pct_115)
cat(sprintf("PCT 115 total area: %.2f sq ft (%.4f sq miles)\n",
            as.numeric(pct_115_area),
            as.numeric(pct_115_area) / 27878400))
cat(sprintf("Zone as %% of PCT 115: %.2f%%\n",
            100 * as.numeric(zone_in_115) / as.numeric(pct_115_area)))

#==============================================================================
# 3. SPATIAL FILTERING FUNCTIONS
#==============================================================================

# Function to filter sf objects to zone
filter_to_zone <- function(sf_obj) {
  sf_obj %>%
    st_filter(zone, .predicate = st_intersects)
}

# Function to filter sf objects to 250 ft buffer (excluding zone)
filter_to_buffer <- function(sf_obj) {
  sf_obj %>%
    st_filter(buffer_250, .predicate = st_intersects)
}

# Function to filter to precinct 115 boundary (within 250 ft)
filter_to_pct115_buffer <- function(sf_obj) {
  pct_115_buffered <- pct_115 %>% st_buffer(250)
  sf_obj %>%
    st_filter(pct_115_buffered, .predicate = st_intersects)
}

#==============================================================================
# 4. FILTER DATA TO ZONES
#==============================================================================

cat("\n=== FILTERING DATA TO ZONES ===\n")

# Complaints in zone (from 2019)
complaints_zone <- complaints_sf %>%
  mutate(date = mdy(rpt_dt)) %>%
  filter(date >= data_start_2019) %>%
  filter_to_zone()
cat(sprintf("Complaints in zone (2019+): %d\n", nrow(complaints_zone)))

# Complaints in 250ft buffer
complaints_buffer <- complaints_sf %>%
  mutate(date = mdy(rpt_dt)) %>%
  filter(date >= data_start_2019) %>%
  filter_to_buffer()
cat(sprintf("Complaints in buffer (2019+): %d\n", nrow(complaints_buffer)))

# Arrests in zone (from 2022)
arrests_zone <- arrests_sf %>%
  mutate(date = mdy(arrest_date)) %>%
  filter(date >= data_start_2022) %>%
  filter_to_zone()
cat(sprintf("Arrests in zone (2022+): %d\n", nrow(arrests_zone)))

# Arrests in buffer
arrests_buffer <- arrests_sf %>%
  mutate(date = mdy(arrest_date)) %>%
  filter(date >= data_start_2022) %>%
  filter_to_buffer()
cat(sprintf("Arrests in buffer (2022+): %d\n", nrow(arrests_buffer)))

# Directed patrols in zone (from 2022)
directed_patrols_zone <- directed_patrols_sf %>%
  mutate(date = mdy(incident_date)) %>%
  filter(date >= data_start_2022) %>%
  filter_to_zone()
cat(sprintf("Directed patrols in zone (2022+): %d\n", nrow(directed_patrols_zone)))

# Summonses in zone (from 2022)
summonses_zone <- all_summonses %>%
  filter(date >= data_start_2022) %>%
  filter_to_zone()
cat(sprintf("Summonses in zone (2022+): %d\n", nrow(summonses_zone)))

# Street violent crime in zone
street_violent_zone <- street_violent_crime %>%
  filter(date >= data_start_2019) %>%
  filter_to_zone()
cat(sprintf("Street violent crimes in zone (2019+): %d\n", nrow(street_violent_zone)))

#==============================================================================
# 5. CRIME PLOTS BY TYPE IN ZONE (FROM 2019)
#==============================================================================

# Define crime categories
violent_ky  <- c(101, 104, 105, 106, 344)
property_ky <- c(107, 109, 110)

# Create monthly crime counts by type
monthly_crime_zone <- complaints_zone %>%
  st_drop_geometry() %>%
  mutate(
    month = floor_date(date, "month"),
    crime_type = case_when(
      ky_cd %in% violent_ky ~ "Violent",
      ky_cd %in% property_ky ~ "Property",
      TRUE ~ "Other"
    )
  ) %>%
  count(month, crime_type) %>%
  filter(crime_type != "Other")

# Plot: Monthly crime by type in zone
p_crime_zone <- ggplot(monthly_crime_zone, aes(x = month, y = n, color = crime_type)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = op_start, linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_vline(xintercept = op_end, linetype = "dashed", color = "red", linewidth = 0.8) +
  annotate("rect", xmin = op_start, xmax = op_end, ymin = -Inf, ymax = Inf,
           alpha = 0.1, fill = "red") +
  labs(
    title = "Monthly Crime Counts in Roosevelt Ave Zone",
    subtitle = "Shaded area = Operation Restore Roosevelt (Oct 15, 2024 - Jan 12, 2025)",
    x = "Month",
    y = "Count",
    color = "Crime Type"
  ) +
  scale_color_brewer(palette = "Set1") +
  theme_minimal() +
  theme(legend.position = "bottom")

p_crime_zone

ggsave(here("output", "eda_plots", "01_crime_by_type_zone.png"), p_crime_zone, 
       width = 10, height = 6, dpi = 300)


# More granular: specific offense types
monthly_crime_detail <- complaints_zone %>%
  st_drop_geometry() %>%
  mutate(month = floor_date(date, "month")) %>%
  filter(ky_cd %in% c(violent_ky, property_ky)) %>%
  count(month, ofns_desc) %>%
  group_by(ofns_desc) %>%
  filter(sum(n) >= 50) %>%  # Only show categories with enough data

ungroup()

p_crime_detail <- ggplot(monthly_crime_detail, aes(x = month, y = n, color = ofns_desc)) +
  geom_line(linewidth = 0.6, alpha = 0.8) +
  geom_vline(xintercept = op_start, linetype = "dashed", color = "red") +
  annotate("rect", xmin = op_start, xmax = op_end, ymin = -Inf, ymax = Inf,
           alpha = 0.1, fill = "red") +
  labs(
    title = "Monthly Crime by Offense Type in Zone",
    x = "Month", y = "Count", color = "Offense"
  ) +
  theme_minimal() +
  theme(legend.position = "right", legend.text = element_text(size = 7))

p_crime_detail
ggsave(here("output", "eda_plots", "02_crime_detail_zone.png"), p_crime_detail, 
       width = 12, height = 6, dpi = 300)

#==============================================================================
# 6. CRIME PLOTS IN 250FT BUFFER (EXCLUDING ZONE)
#==============================================================================

monthly_crime_buffer <- complaints_buffer %>%
  st_drop_geometry() %>%
  mutate(
    month = floor_date(date, "month"),
    crime_type = case_when(
      ky_cd %in% violent_ky ~ "Violent",
      ky_cd %in% property_ky ~ "Property",
      TRUE ~ "Other"
    )
  ) %>%
  count(month, crime_type) %>%
  filter(crime_type != "Other")

p_crime_buffer <- ggplot(monthly_crime_buffer, aes(x = month, y = n, color = crime_type)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = op_start, linetype = "dashed", color = "red") +
  annotate("rect", xmin = op_start, xmax = op_end, ymin = -Inf, ymax = Inf,
           alpha = 0.1, fill = "red") +
  labs(
    title = "Monthly Crime in 250ft Buffer Around Zone (Displacement Area)",
    subtitle = "Excludes zone itself",
    x = "Month", y = "Count", color = "Crime Type"
  ) +
  scale_color_brewer(palette = "Set1") +
  theme_minimal() +
  theme(legend.position = "bottom")

p_crime_buffer

ggsave(here("output", "eda_plots", "03_crime_by_type_buffer.png"), p_crime_buffer, 
       width = 10, height = 6, dpi = 300)

# Combined zone vs buffer comparison
combined_crime <- bind_rows(
  monthly_crime_zone %>% mutate(location = "Zone"),
  monthly_crime_buffer %>% mutate(location = "250ft Buffer")
)

p_zone_buffer_compare <- ggplot(combined_crime, 
                                 aes(x = month, y = n, color = crime_type, linetype = location)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = op_start, linetype = "dashed", color = "gray40") +
  annotate("rect", xmin = op_start, xmax = op_end, ymin = -Inf, ymax = Inf,
           alpha = 0.1, fill = "red") +
  facet_wrap(~crime_type, scales = "free_y") +
  labs(
    title = "Crime Trends: Zone vs. 250ft Buffer",
    x = "Month", y = "Count", 
    color = "Crime Type", linetype = "Location"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

p_zone_buffer_compare

ggsave(here("output", "eda_plots", "04_zone_buffer_comparison.png"), p_zone_buffer_compare, 
       width = 12, height = 6, dpi = 300)

#==============================================================================
# 7. ARRESTS IN ZONE (2022 TO PRESENT)
#==============================================================================

monthly_arrests_zone <- arrests_zone %>%
  st_drop_geometry() %>%
  mutate(month = floor_date(date, "month")) %>%
  count(month)

p_arrests_zone <- ggplot(monthly_arrests_zone, aes(x = month, y = n)) +
  geom_line(linewidth = 0.8, color = "steelblue") +
  geom_point(size = 1.5, color = "steelblue") +
  geom_vline(xintercept = op_start, linetype = "dashed", color = "red") +
  annotate("rect", xmin = op_start, xmax = op_end, ymin = -Inf, ymax = Inf,
           alpha = 0.1, fill = "red") +
  labs(
    title = "Monthly Arrests in Roosevelt Ave Zone",
    subtitle = "2022 to Present",
    x = "Month", y = "Arrests"
  ) +
  theme_minimal()

p_arrests_zone

ggsave(here("output", "eda_plots", "05_arrests_zone.png"), p_arrests_zone, 
       width = 10, height = 6, dpi = 300)

#==============================================================================
# 8. ROBBERY ARRESTS IN ZONE
#==============================================================================

robbery_arrests_zone <- arrests_zone %>%
  st_drop_geometry() %>%
  filter(ofns_desc == "ROBBERY") %>%
  mutate(month = floor_date(date, "month")) %>%
  count(month)

p_robbery <- ggplot(robbery_arrests_zone, aes(x = month, y = n)) +
  geom_col(fill = "darkred", alpha = 0.8) +
  geom_vline(xintercept = op_start, linetype = "dashed", color = "black") +
  annotate("rect", xmin = op_start, xmax = op_end, ymin = -Inf, ymax = Inf,
           alpha = 0.1, fill = "red") +
  labs(
    title = "Monthly Robbery Arrests in Zone",
    subtitle = "2022 to Present",
    x = "Month", y = "Arrests"
  ) +
  theme_minimal()

p_robbery

ggsave(here("output", "eda_plots", "06_robbery_arrests_zone.png"), p_robbery, 
       width = 10, height = 6, dpi = 300)

#==============================================================================
# 9. ASSAULT ARRESTS IN ZONE (ky_cd 106, 344)
#==============================================================================

assault_arrests_zone <- arrests_zone %>%
  st_drop_geometry() %>%
  filter(ky_cd %in% c(106, 344)) %>%
  mutate(
    month = floor_date(date, "month"),
    assault_type = case_when(
      ky_cd == 106 ~ "Felony Assault",
      ky_cd == 344 ~ "Assault 3 & Related",
      TRUE ~ "Other"
    )
  ) %>%
  count(month, assault_type)

p_assault <- ggplot(assault_arrests_zone, aes(x = month, y = n, fill = assault_type)) +
  geom_col(position = "stack", alpha = 0.8) +
  geom_vline(xintercept = op_start, linetype = "dashed", color = "black") +
  annotate("rect", xmin = op_start, xmax = op_end, ymin = -Inf, ymax = Inf,
           alpha = 0.1, fill = "red") +
  labs(
    title = "Monthly Assault Arrests in Zone (ky_cd 106, 344)",
    subtitle = "2022 to Present",
    x = "Month", y = "Arrests", fill = "Type"
  ) +
  scale_fill_brewer(palette = "Reds") +
  theme_minimal() +
  theme(legend.position = "bottom")

p_assault

ggsave(here("output", "eda_plots", "07_assault_arrests_zone.png"), p_assault, 
       width = 10, height = 6, dpi = 300)

#==============================================================================
# 10. DIRECTED PATROLS IN ZONE
#==============================================================================

monthly_patrols_zone <- directed_patrols_zone %>%
  st_drop_geometry() %>%
  mutate(month = floor_date(date, "month")) %>%
  count(month)

p_patrols <- ggplot(monthly_patrols_zone, aes(x = month, y = n)) +
  geom_line(linewidth = 0.8, color = "darkgreen") +
  geom_point(size = 1.5, color = "darkgreen") +
  geom_vline(xintercept = op_start, linetype = "dashed", color = "red") +
  annotate("rect", xmin = op_start, xmax = op_end, ymin = -Inf, ymax = Inf,
           alpha = 0.1, fill = "red") +
  labs(
    title = "Monthly Directed Patrols in Zone",
    subtitle = "2022 to Present (Measure of Police Effort)",
    x = "Month", y = "Directed Patrol Events"
  ) +
  theme_minimal()

p_patrols

ggsave(here("output", "eda_plots", "08_directed_patrols_zone.png"), p_patrols, 
       width = 10, height = 6, dpi = 300)

#==============================================================================
# 11. SUMMONSES IN ZONE
#==============================================================================

monthly_summonses_zone <- summonses_zone %>%
  st_drop_geometry() %>%
  mutate(month = floor_date(date, "month")) %>%
  count(month, summons_type)

p_summonses <- ggplot(monthly_summonses_zone, aes(x = month, y = n, fill = summons_type)) +
  geom_col(position = "stack", alpha = 0.8) +
  geom_vline(xintercept = op_start, linetype = "dashed", color = "black") +
  annotate("rect", xmin = op_start, xmax = op_end, ymin = -Inf, ymax = Inf,
           alpha = 0.1, fill = "red") +
  labs(
    title = "Monthly Summonses in Zone",
    subtitle = "Criminal Court + OATH Summonses, 2022 to Present",
    x = "Month", y = "Summonses", fill = "Type"
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  theme(legend.position = "bottom")

p_summonses

ggsave(here("output", "eda_plots", "09_summonses_zone.png"), p_summonses, 
       width = 10, height = 6, dpi = 300)

#==============================================================================
# 12. OPERATION PERIOD VS PRIOR YEAR COMPARISON
#==============================================================================

cat("\n=== OPERATION vs PRIOR YEAR COMPARISON ===\n")

# Arrests comparison
arrests_during <- arrests_zone %>%
  st_drop_geometry() %>%
  filter(date >= op_start, date <= op_end) %>%
  nrow()

arrests_prior <- arrests_zone %>%
  st_drop_geometry() %>%
  filter(date >= comp_start, date <= comp_end) %>%
  nrow()

cat(sprintf("Arrests during operation: %d\n", arrests_during))
cat(sprintf("Arrests same period prior year: %d\n", arrests_prior))
cat(sprintf("Change: %+d (%.1f%%)\n", 
            arrests_during - arrests_prior,
            100 * (arrests_during - arrests_prior) / max(arrests_prior, 1)))

# Summonses comparison
summonses_during <- summonses_zone %>%
  st_drop_geometry() %>%
  filter(date >= op_start, date <= op_end) %>%
  nrow()

summonses_prior <- summonses_zone %>%
  st_drop_geometry() %>%
  filter(date >= comp_start, date <= comp_end) %>%
  nrow()

cat(sprintf("\nSummonses during operation: %d\n", summonses_during))
cat(sprintf("Summonses same period prior year: %d\n", summonses_prior))
cat(sprintf("Change: %+d (%.1f%%)\n", 
            summonses_during - summonses_prior,
            100 * (summonses_during - summonses_prior) / max(summonses_prior, 1)))

# Directed patrols comparison
patrols_during <- directed_patrols_zone %>%
  st_drop_geometry() %>%
  filter(date >= op_start, date <= op_end) %>%
  nrow()

patrols_prior <- directed_patrols_zone %>%
  st_drop_geometry() %>%
  filter(date >= comp_start, date <= comp_end) %>%
  nrow()

cat(sprintf("\nDirected patrols during operation: %d\n", patrols_during))
cat(sprintf("Directed patrols same period prior year: %d\n", patrols_prior))
cat(sprintf("Change: %+d (%.1f%%)\n", 
            patrols_during - patrols_prior,
            100 * (patrols_during - patrols_prior) / max(patrols_prior, 1)))

# Create comparison summary table
comparison_df <- tibble(
  Measure = c("Arrests", "Summonses", "Directed Patrols"),
  During_ORR = c(arrests_during, summonses_during, patrols_during),
  Prior_Year = c(arrests_prior, summonses_prior, patrols_prior),
  Change = c(arrests_during - arrests_prior, 
             summonses_during - summonses_prior,
             patrols_during - patrols_prior),
  Pct_Change = c(
    100 * (arrests_during - arrests_prior) / max(arrests_prior, 1),
    100 * (summonses_during - summonses_prior) / max(summonses_prior, 1),
    100 * (patrols_during - patrols_prior) / max(patrols_prior, 1)
  )
)

write_csv(comparison_df, here("output", "eda_plots", "operation_comparison.csv"))

# Bar plot comparison
comparison_long <- comparison_df %>%
  pivot_longer(cols = c(During_ORR, Prior_Year), 
               names_to = "Period", values_to = "Count") %>%
  mutate(Period = recode(Period, 
                         "During_ORR" = "During Operation",
                         "Prior_Year" = "Prior Year (Same Dates)"))

p_comparison <- ggplot(comparison_long, aes(x = Measure, y = Count, fill = Period)) +
  geom_col(position = "dodge", alpha = 0.8) +
  geom_text(aes(label = Count), position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 3) +
  labs(
    title = "Operation Period vs. Prior Year Comparison",
    subtitle = sprintf("Operation: %s to %s", op_start, op_end),
    x = "", y = "Count", fill = ""
  ) +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  theme(legend.position = "bottom")

p_comparison

ggsave(here("output", "eda_plots", "10_operation_comparison.png"), p_comparison, 
       width = 10, height = 6, dpi = 300)

#==============================================================================
# 13. WEEKLY AND MONTHLY STREET VIOLENCE DISTRIBUTIONS
#==============================================================================

# Weekly counts
weekly_street_violence <- street_violent_zone %>%
  st_drop_geometry() %>%
  mutate(week = floor_date(date, "week")) %>%
  count(week) %>%
  filter(week >= data_start_2019)

# Monthly counts
monthly_street_violence <- street_violent_zone %>%
  st_drop_geometry() %>%
  mutate(month = floor_date(date, "month")) %>%
  count(month) %>%
  filter(month >= data_start_2019)

cat("\n=== STREET VIOLENCE DISTRIBUTIONS IN ZONE ===\n")
cat("\nWeekly Street Violence:\n")
cat(sprintf("  Mean: %.2f\n", mean(weekly_street_violence$n)))
cat(sprintf("  Median: %.0f\n", median(weekly_street_violence$n)))
cat(sprintf("  Min: %d\n", min(weekly_street_violence$n)))
cat(sprintf("  Max: %d\n", max(weekly_street_violence$n)))
cat(sprintf("  SD: %.2f\n", sd(weekly_street_violence$n)))

cat("\nMonthly Street Violence:\n")
cat(sprintf("  Mean: %.2f\n", mean(monthly_street_violence$n)))
cat(sprintf("  Median: %.0f\n", median(monthly_street_violence$n)))
cat(sprintf("  Min: %d\n", min(monthly_street_violence$n)))
cat(sprintf("  Max: %d\n", max(monthly_street_violence$n)))
cat(sprintf("  SD: %.2f\n", sd(monthly_street_violence$n)))

# Weekly time series
p_weekly_ts <- ggplot(weekly_street_violence, aes(x = week, y = n)) +
  geom_line(color = "darkred", alpha = 0.7) +
  geom_smooth(method = "loess", span = 0.2, se = FALSE, color = "blue") +
  geom_vline(xintercept = op_start, linetype = "dashed", color = "black") +
  annotate("rect", xmin = op_start, xmax = op_end, ymin = -Inf, ymax = Inf,
           alpha = 0.1, fill = "red") +
  labs(
    title = "Weekly Street Violence in Zone",
    subtitle = "Blue line = LOESS smoother",
    x = "Week", y = "Incidents"
  ) +
  theme_minimal()

p_weekly_ts

# Distribution plots
p_weekly_dist <- ggplot(weekly_street_violence, aes(x = n)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "white", alpha = 0.8) +
  geom_vline(xintercept = mean(weekly_street_violence$n), 
             linetype = "dashed", color = "red") +
  labs(
    title = "Distribution of Weekly Street Violence Counts",
    subtitle = sprintf("Mean = %.1f (red line)", mean(weekly_street_violence$n)),
    x = "Weekly Count", y = "Frequency"
  ) +
  theme_minimal()

p_weekly_dist
p_weekly_combined <- p_weekly_ts / p_weekly_dist
p_weekly_combined
ggsave(here("output", "eda_plots", "11_weekly_street_violence.png"), p_weekly_combined, 
       width = 10, height = 10, dpi = 300)

# Monthly time series and distribution
p_monthly_ts <- ggplot(monthly_street_violence, aes(x = month, y = n)) +
  geom_line(color = "darkred", linewidth = 0.8) +
  geom_point(color = "darkred", size = 2) +
  geom_vline(xintercept = op_start, linetype = "dashed", color = "black") +
  annotate("rect", xmin = op_start, xmax = op_end, ymin = -Inf, ymax = Inf,
           alpha = 0.1, fill = "red") +
  labs(
    title = "Monthly Street Violence in Zone",
    x = "Month", y = "Incidents"
  ) +
  theme_minimal()

p_monthly_ts

p_monthly_dist <- ggplot(monthly_street_violence, aes(x = n)) +
  geom_histogram(binwidth = 2, fill = "steelblue", color = "white", alpha = 0.8) +
  geom_vline(xintercept = mean(monthly_street_violence$n), 
             linetype = "dashed", color = "red") +
  labs(
    title = "Distribution of Monthly Street Violence Counts",
    subtitle = sprintf("Mean = %.1f", mean(monthly_street_violence$n)),
    x = "Monthly Count", y = "Frequency"
  ) +
  theme_minimal()

p_monthly_dist

p_monthly_combined <- p_monthly_ts / p_monthly_dist

p_monthly_combined
ggsave(here("output", "eda_plots", "12_monthly_street_violence.png"), p_monthly_combined, 
       width = 10, height = 10, dpi = 300)

#==============================================================================
# 14. CENSUS GEOGRAPHY VIOLENT CRIME DISTRIBUTIONS (DONOR POOL ASSESSMENT)
#==============================================================================

cat("\n=== CENSUS GEOGRAPHY CRIME DISTRIBUTIONS ===\n")
cat("(For assessing synthetic control donor pool appropriateness)\n")

# Join street violent crime to block groups
street_violent_bg <- street_violent_crime %>%
  filter(date >= data_start_2022) %>%
  st_join(nyc_bgs %>% dplyr::select(geoid)) %>%
  st_drop_geometry() %>%
  filter(!is.na(geoid))

# Monthly counts by block group
bg_monthly <- street_violent_bg %>%
  mutate(month = floor_date(date, "month")) %>%
  count(geoid, month)

# Summary stats per block group
bg_summary <- bg_monthly %>%
  group_by(geoid) %>%
  summarize(
    total_months = n(),
    total_crimes = sum(n),
    mean_monthly = mean(n),
    median_monthly = median(n),
    max_monthly = max(n),
    .groups = "drop"
  )

cat("\nBlock Group Street Violence Summary (2022+):\n")
cat(sprintf("  Total block groups with any street violence: %d\n", nrow(bg_summary)))
cat(sprintf("  Mean monthly crimes per BG: %.2f\n", mean(bg_summary$mean_monthly)))
cat(sprintf("  Median monthly crimes per BG: %.2f\n", median(bg_summary$mean_monthly)))
cat(sprintf("  Max mean monthly: %.2f\n", max(bg_summary$mean_monthly)))

# What does the zone's block groups look like?
zone_bgs <- nyc_bgs %>%
  st_filter(zone, .predicate = st_intersects) %>%
  pull(geoid)

zone_bg_stats <- bg_summary %>%
  filter(geoid %in% zone_bgs)

cat("\nZone Block Groups:\n")
cat(sprintf("  Number of BGs in zone: %d\n", length(zone_bgs)))
cat(sprintf("  Mean monthly in zone BGs: %.2f\n", mean(zone_bg_stats$mean_monthly)))
cat(sprintf("  Range: %.1f - %.1f\n", min(zone_bg_stats$mean_monthly), max(zone_bg_stats$mean_monthly)))

# Distribution of block group means
p_bg_dist <- ggplot(bg_summary, aes(x = mean_monthly)) +
  geom_histogram(binwidth = 0.5, fill = "steelblue", color = "white", alpha = 0.8) +
  geom_vline(xintercept = mean(zone_bg_stats$mean_monthly), 
             linetype = "dashed", color = "red", linewidth = 1) +
  labs(
    title = "Distribution of Mean Monthly Street Violence by Block Group",
    subtitle = sprintf("Red line = Zone BG average (%.1f)", mean(zone_bg_stats$mean_monthly)),
    x = "Mean Monthly Street Violence", y = "Number of Block Groups"
  ) +
  theme_minimal()

p_bg_dist

ggsave(here("output", "eda_plots", "13_bg_distribution.png"), p_bg_dist, 
       width = 10, height = 6, dpi = 300)

# Census tracts
street_violent_ct <- street_violent_crime %>%
  filter(date >= data_start_2022) %>%
  st_join(nyct %>% dplyr::select(boro_ct2020)) %>%
  st_drop_geometry() %>%
  filter(!is.na(boro_ct2020))

ct_monthly <- street_violent_ct %>%
  mutate(month = floor_date(date, "month")) %>%
  count(boro_ct2020, month)

ct_summary <- ct_monthly %>%
  group_by(boro_ct2020) %>%
  summarize(
    total_months = n(),
    total_crimes = sum(n),
    mean_monthly = mean(n),
    median_monthly = median(n),
    max_monthly = max(n),
    .groups = "drop"
  )

cat("\nCensus Tract Street Violence Summary (2022+):\n")
cat(sprintf("  Total tracts with any street violence: %d\n", nrow(ct_summary)))
cat(sprintf("  Mean monthly crimes per tract: %.2f\n", mean(ct_summary$mean_monthly)))
cat(sprintf("  Median monthly crimes per tract: %.2f\n", median(ct_summary$mean_monthly)))

# Zone tracts
zone_cts <- nyct %>%
  st_filter(zone, .predicate = st_intersects) %>%
  pull(boro_ct2020)

zone_ct_stats <- ct_summary %>%
  filter(boro_ct2020 %in% zone_cts)

cat("\nZone Census Tracts:\n")
cat(sprintf("  Number of tracts in zone: %d\n", length(zone_cts)))
cat(sprintf("  Mean monthly in zone tracts: %.2f\n", mean(zone_ct_stats$mean_monthly)))

p_ct_dist <- ggplot(ct_summary, aes(x = mean_monthly)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "white", alpha = 0.8) +
  geom_vline(xintercept = mean(zone_ct_stats$mean_monthly), 
             linetype = "dashed", color = "red", linewidth = 1) +
  labs(
    title = "Distribution of Mean Monthly Street Violence by Census Tract",
    subtitle = sprintf("Red line = Zone tract average (%.1f)", mean(zone_ct_stats$mean_monthly)),
    x = "Mean Monthly Street Violence", y = "Number of Tracts"
  ) +
  theme_minimal()

p_ct_dist

ggsave(here("output", "eda_plots", "14_ct_distribution.png"), p_ct_dist, 
       width = 10, height = 6, dpi = 300)

# Save summaries for donor pool selection
write_csv(bg_summary, here("output", "eda_plots", "block_group_violence_summary.csv"))
write_csv(ct_summary, here("output", "eda_plots", "census_tract_violence_summary.csv"))

#==============================================================================
# 15. ARREST TRENDS WITHIN 250FT OF PCT 115
#==============================================================================

# Get PCT 115 buffered
pct_115_buffered <- pct_115 %>% st_buffer(250)

arrests_pct115_area <- arrests_sf %>%
  mutate(date = mdy(arrest_date)) %>%
  filter(date >= data_start_2022) %>%
  st_filter(pct_115_buffered, .predicate = st_intersects)

monthly_arrests_pct115 <- arrests_pct115_area %>%
  st_drop_geometry() %>%
  mutate(month = floor_date(date, "month")) %>%
  count(month)

p_pct115_arrests <- ggplot(monthly_arrests_pct115, aes(x = month, y = n)) +
  geom_line(linewidth = 0.8, color = "steelblue") +
  geom_point(size = 1.5, color = "steelblue") +
  geom_vline(xintercept = op_start, linetype = "dashed", color = "red") +
  annotate("rect", xmin = op_start, xmax = op_end, ymin = -Inf, ymax = Inf,
           alpha = 0.1, fill = "red") +
  labs(
    title = "Monthly Arrests Within 250ft of PCT 115 Boundary",
    subtitle = "2022 to Present",
    x = "Month", y = "Arrests"
  ) +
  theme_minimal()

p_pct115_arrests

ggsave(here("output", "eda_plots", "15_pct115_arrests.png"), p_pct115_arrests, 
       width = 10, height = 6, dpi = 300)

cat(sprintf("\n=== PCT 115 AREA ARRESTS ===\n"))
cat(sprintf("Total arrests in PCT 115 area (2022+): %d\n", nrow(arrests_pct115_area)))


# Get PCT 110 buffered
pct_110_buffered <- pct_110 %>% st_buffer(250)

arrests_pct110_area <- arrests_sf %>%
  mutate(date = mdy(arrest_date)) %>%
  filter(date >= data_start_2022) %>%
  st_filter(pct_110_buffered, .predicate = st_intersects)

monthly_arrests_pct110 <- arrests_pct110_area %>%
  st_drop_geometry() %>%
  mutate(month = floor_date(date, "month")) %>%
  count(month)

p_pct110_arrests <- ggplot(monthly_arrests_pct110, aes(x = month, y = n)) +
  geom_line(linewidth = 0.8, color = "steelblue") +
  geom_point(size = 1.5, color = "steelblue") +
  geom_vline(xintercept = op_start, linetype = "dashed", color = "red") +
  annotate("rect", xmin = op_start, xmax = op_end, ymin = -Inf, ymax = Inf,
           alpha = 0.1, fill = "red") +
  labs(
    title = "Monthly Arrests Within 250ft of PCT 110 Boundary",
    subtitle = "2022 to Present",
    x = "Month", y = "Arrests"
  ) +
  theme_minimal()

p_pct110_arrests

ggsave(here("output", "eda_plots", "15_pct110_arrests.png"), p_pct115_arrests, 
       width = 10, height = 6, dpi = 300)

cat(sprintf("\n=== PCT 110 AREA ARRESTS ===\n"))
cat(sprintf("Total arrests in PCT 115 area (2022+): %d\n", nrow(arrests_pct110_area)))



#==============================================================================
# 16. TOP ARREST TYPES: 90 DAYS BEFORE VS DURING OPERATION
#==============================================================================

# Pre-operation (90 days before)
arrests_pre <- arrests_zone %>%
  st_drop_geometry() %>%
  filter(date >= pre_90_start, date <= pre_90_end) %>%
  count(pd_desc, sort = TRUE) %>%
  mutate(period = "Pre-Operation (90 days)")

# During operation
arrests_during_detail <- arrests_zone %>%
  st_drop_geometry() %>%
  filter(date >= op_start, date <= op_end) %>%
  count(pd_desc, sort = TRUE) %>%
  mutate(period = "During Operation")

cat("\n=== TOP ARREST TYPES ===\n")
cat("\nTop 15 - 90 Days Before Operation:\n")
print(head(arrests_pre, 15))

cat("\nTop 15 - During Operation:\n")
print(head(arrests_during_detail, 15))

# Combined comparison of top types
top_types <- union(head(arrests_pre$pd_desc, 10), head(arrests_during_detail$pd_desc, 10))

arrest_comparison <- bind_rows(arrests_pre, arrests_during_detail) %>%
  filter(pd_desc %in% top_types) %>%
  pivot_wider(names_from = period, values_from = n, values_fill = 0) %>%
  mutate(change = `During Operation` - `Pre-Operation (90 days)`) %>%
  arrange(desc(`During Operation`))

write_csv(arrest_comparison, here("output", "eda_plots", "arrest_type_comparison.csv"))

# Plot
arrest_comparison_long <- bind_rows(arrests_pre, arrests_during_detail) %>%
  filter(pd_desc %in% head(top_types, 12))

p_arrest_types <- ggplot(arrest_comparison_long, 
                          aes(x = reorder(pd_desc, n), y = n, fill = period)) +
  geom_col(position = "dodge", alpha = 0.8) +
  coord_flip() +
  labs(
    title = "Top Arrest Types: Pre-Operation vs During Operation",
    x = "", y = "Count", fill = ""
  ) +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  theme(legend.position = "bottom")

p_arrest_types

ggsave(here("output", "eda_plots", "16_arrest_types_comparison.png"), p_arrest_types, 
       width = 12, height = 8, dpi = 300)

library(leaflet)
library(scales)
library(htmltools)


# Top 50 BGs by prostitution arrests
top50_bg <- pros_bg %>% slice_max(n, n = 50, with_ties = FALSE)

# Join counts to BG polygons (nyc_bgs is in 2263)
top50_bg_sf <- nyc_bgs %>%
  select(geoid) %>%
  inner_join(top50_bg, by = "geoid") %>%
  mutate(
    label = paste0("BG: ", geoid, "<br>Prostitution arrests: ", n)
  ) %>%
  st_transform(4326)

pal <- colorNumeric("YlOrRd", domain = top50_bg_sf$n)

leaflet(top50_bg_sf) %>%
  addTiles() %>%
  addPolygons(
    color = ~pal(n),
    weight = 1,
    opacity = 1,
    fillColor = ~pal(n),
    fillOpacity = 0.6,
    label = ~label,
    highlightOptions = highlightOptions(weight = 3, bringToFront = TRUE)
  ) %>%
  addLegend(
    pal = pal,
    values = ~n,
    title = htmltools::HTML("Prostitution arrests<br>(Top 50 block groups)"),
    opacity = 0.8
  )








#==============================================================================
# 17. PROSTITUTION HOT SPOTS (CENSUS GEOGRAPHIES)
#==============================================================================

cat("\n=== PROSTITUTION ARREST HOT SPOTS ===\n")

# Prostitution arrests (PL 230)
prostitution_arrests <- arrests_sf %>%
  mutate(date = mdy(arrest_date)) %>%
  filter(str_detect(law_code, "PL 230"))

cat(sprintf("Total prostitution arrests (all time): %d\n", nrow(prostitution_arrests)))

# By block group
pros_bg <- prostitution_arrests %>%
  st_join(nyc_bgs %>% dplyr::select(geoid)) %>%
  st_drop_geometry() %>%
  filter(!is.na(geoid)) %>%
  count(geoid, sort = TRUE)

cat("\nTop 20 Block Groups by Prostitution Arrests:\n")
print(head(pros_bg, 20))

# By census tract
pros_ct <- prostitution_arrests %>%
  st_join(nyct %>% dplyr::select(boro_ct2020, ct_label)) %>%
  st_drop_geometry() %>%
  filter(!is.na(boro_ct2020)) %>%
  count(boro_ct2020, ct_label, sort = TRUE)

cat("\nTop 20 Census Tracts by Prostitution Arrests:\n")
print(head(pros_ct, 20))

# Save for spillover analysis
write_csv(pros_bg, here("output", "eda_plots", "prostitution_arrests_by_bg.csv"))
write_csv(pros_ct, here("output", "eda_plots", "prostitution_arrests_by_tract.csv"))

# How many are in the zone?
pros_in_zone <- prostitution_arrests %>%
  filter_to_zone() %>%
  nrow()
cat(sprintf("\nProstitution arrests in zone: %d (%.1f%% of total)\n", 
            pros_in_zone, 
            100 * pros_in_zone / nrow(prostitution_arrests)))

#==============================================================================
# 18. PROSTITUTION ARRESTS IN ZONE (2022 ONWARD)
#==============================================================================

prostitution_zone <- prostitution_arrests %>%
  filter(date >= data_start_2022) %>%
  filter_to_zone()

monthly_pros_zone <- prostitution_zone %>%
  st_drop_geometry() %>%
  mutate(month = floor_date(date, "month")) %>%
  count(month)

p_prostitution <- ggplot(monthly_pros_zone, aes(x = month, y = n)) +
  geom_col(fill = "purple", alpha = 0.8) +
  geom_vline(xintercept = op_start, linetype = "dashed", color = "black") +
  annotate("rect", xmin = op_start, xmax = op_end, ymin = -Inf, ymax = Inf,
           alpha = 0.1, fill = "red") +
  labs(
    title = "Monthly Prostitution Arrests in Zone (PL 230)",
    subtitle = "2022 to Present",
    x = "Month", y = "Arrests"
  ) +
  theme_minimal()

p_prostitution

ggsave(here("output", "eda_plots", "17_prostitution_arrests_zone.png"), p_prostitution, 
       width = 10, height = 6, dpi = 300)

cat(sprintf("\nProstitution arrests in zone (2022+): %d\n", nrow(prostitution_zone)))

#==============================================================================
# 19. ADDITIONAL: PRE-INTERVENTION SPIKE ANALYSIS
#==============================================================================

# Look at the 6 months before intervention to understand the spike
spike_window_start <- op_start - months(6)
spike_window_end <- op_start - 1

weekly_pre_spike <- street_violent_zone %>%
  st_drop_geometry() %>%
  filter(date >= spike_window_start, date <= spike_window_end) %>%
  mutate(week = floor_date(date, "week")) %>%
  count(week)

p_pre_spike <- ggplot(weekly_pre_spike, aes(x = week, y = n)) +
  geom_line(color = "darkred", linewidth = 0.8) +
  geom_point(color = "darkred", size = 2) +
  geom_smooth(method = "loess", span = 0.5, se = TRUE, color = "blue", alpha = 0.2) +
  geom_vline(xintercept = op_start, linetype = "dashed", color = "black") +
  labs(
    title = "Weekly Street Violence - 6 Months Pre-Intervention",
    subtitle = "Examining the pre-intervention spike pattern",
    x = "Week", y = "Incidents"
  ) +
  theme_minimal()

p_pre_spike

ggsave(here("output", "eda_plots", "18_pre_spike_analysis.png"), p_pre_spike, 
       width = 10, height = 6, dpi = 300)

# Month-over-month comparison for same period
mom_comparison <- street_violent_zone %>%
  st_drop_geometry() %>%
  filter(month(date) %in% c(4:10)) %>%  # Apr-Oct
mutate(
  year = year(date),
  month_num = month(date)
) %>%
  count(year, month_num) %>%
  filter(year >= 2022)

p_mom <- ggplot(mom_comparison, aes(x = factor(month_num), y = n, 
                                     color = factor(year), group = year)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  labs(
    title = "Monthly Street Violence by Year (Apr-Oct)",
    subtitle = "Comparing seasonal patterns across years",
    x = "Month", y = "Incidents", color = "Year"
  ) +
  scale_x_discrete(labels = month.abb[4:10]) +
  scale_color_brewer(palette = "Set1") +
  theme_minimal() +
  theme(legend.position = "bottom")

p_mom

ggsave(here("output", "eda_plots", "19_year_over_year.png"), p_mom, 
       width = 10, height = 6, dpi = 300)

#==============================================================================
# 20. SUMMARY OUTPUT
#==============================================================================

cat("\n\n")
cat("============================================================\n")
cat("EDA COMPLETE - SUMMARY\n")
cat("============================================================\n")
cat(sprintf("Zone area: %.4f sq miles\n", as.numeric(zone_total) / 27878400))
cat(sprintf("Zone as %% of PCT 115: %.2f%%\n", 
            100 * as.numeric(zone_in_115) / as.numeric(pct_115_area)))
cat(sprintf("\nOperation period: %s to %s (%d days)\n", op_start, op_end, 
            as.integer(op_end - op_start) + 1))
cat(sprintf("\nKey findings saved to: %s\n", here("output", "eda_plots")))
cat("\nPlots generated:\n")
list.files(here("output", "eda_plots"), pattern = "\\.png$") %>%
  paste(" -", .) %>%
  cat(sep = "\n")

cat("\n\nCSV outputs:\n")
list.files(here("output", "eda_plots"), pattern = "\\.csv$") %>%
  paste(" -", .) %>%
  cat(sep = "\n")

#==============================================================================
# SAVE KEY OBJECTS FOR LATER USE
#==============================================================================

saveRDS(list(
  zone = zone,
  buffer_250 = buffer_250,
  pct_115 = pct_115,
  op_start = op_start,
  op_end = op_end,
  zone_bgs = zone_bgs,
  zone_cts = zone_cts
), file = here("output", "eda_objects.rds"))

cat("\nKey spatial objects saved to: output/eda_objects.rds\n")
