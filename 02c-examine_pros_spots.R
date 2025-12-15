library(sf)
library(dplyr)
library(leaflet)

pros_top20_ll <- pros_bgs_ll |>
  st_make_valid() |>
  (\(x) {
    crs <- st_crs(x)
    if (is.na(crs) || is.null(crs$epsg) || crs$epsg != 4326) st_transform(x, 4326) else x
  })() |>
  arrange(desc(n)) |>
  mutate(rank = row_number()) |>
  filter(rank <= 20) |>
  mutate(
    label_pt = st_point_on_surface(geometry),
    popup_html = sprintf(
      "<b>Rank:</b> %s<br/><b>GEOID:</b> %s<br/><b>n:</b> %s",
      rank, .geoid, n
    )
  )

leaflet(options = leafletOptions(preferCanvas = TRUE, zoomControl = TRUE)) |>
  addProviderTiles(providers$CartoDB.Positron) |>
  addPolygons(
    data = pros_top20_ll,
    weight = 2,
    opacity = 1,
    fillOpacity = 0.35,
    popup = ~popup_html
  ) |>
  addLabelOnlyMarkers(
    data = pros_top20_ll,
    lng = ~st_coordinates(label_pt)[, 1],
    lat = ~st_coordinates(label_pt)[, 2],
    label = ~as.character(rank),
    labelOptions = labelOptions(
      noHide = TRUE,
      direction = "center",
      textOnly = TRUE,
      style = list(
        "font-weight" = "700",
        "font-size"   = "14px",
        "color"       = "#111",
        "text-shadow" = "0 0 3px rgba(255,255,255,0.9)"
      )
    )
  )



# ranks you want
ranks_want <- c(9, 14)

# GEOIDs for those ranks
bg_ids <- pros_top20_ll %>%
  st_drop_geometry() %>%
  filter(rank %in% ranks_want) %>%
  pull(.geoid)

# monthly violent crime, COMBINED into one line
bg_monthly_combined <- violent_crime %>%
  filter(date >= as.Date("2022-01-01")) %>%
  mutate(month = floor_date(date, "month")) %>%
  st_join(
    nyc_bgs %>% select(.geoid = geoid),
    left = FALSE
  ) %>%
  st_drop_geometry() %>%
  filter(.geoid %in% bg_ids) %>%
  count(month, name = "crime_count") %>%   # <- key change
  complete(
    month = seq(as.Date("2022-01-01"),
                floor_date(Sys.Date(), "month"),
                by = "month"),
    fill = list(crime_count = 0)
  )

ggplot(bg_monthly_combined, aes(x = month, y = crime_count)) +
  geom_line(linewidth = 1.1, color = "steelblue") +
  geom_point(size = 1.6, color = "steelblue") +
  annotate("rect", xmin = op_start, xmax = op_end, ymin = -Inf, ymax = Inf,
           alpha = 0.1, fill = "red") +
  labs(
    title = "Monthly Street Violent Crime (Combined Block Groups)",
    subtitle = "PROS-ranked block groups: ranks 9 and 14",
    x = "Month",
    y = "Incidents"
  ) +
  theme_minimal()




# define rank groups
rank_groups <- tibble(
  group = c("Ranks 9 & 14", "Ranks 3, 5, 11, 12, 20"),
  ranks = list(
    c(9, 14),
    c(3, 5, 11, 12, 20)
  )
)

# map ranks -> GEOIDs
bg_lookup <- pros_top20_ll %>%
  st_drop_geometry() %>%
  select(.geoid, rank)

# monthly violent crime, aggregated by group
bg_monthly_groups <- violent_crime %>%
  filter(date >= as.Date("2022-01-01")) %>%
  mutate(month = floor_date(date, "month")) %>%
  st_join(
    nyc_bgs %>% select(.geoid = geoid),
    left = FALSE
  ) %>%
  st_drop_geometry() %>%
  inner_join(bg_lookup, by = ".geoid") %>%
  inner_join(
    rank_groups %>% unnest(ranks),
    by = c("rank" = "ranks")
  ) %>%
  count(group, month, name = "crime_count") %>%
  complete(
    group,
    month = seq(
      as.Date("2022-01-01"),
      floor_date(Sys.Date(), "month"),
      by = "month"
    ),
    fill = list(crime_count = 0)
  )

ggplot(bg_monthly_groups, aes(x = month, y = crime_count, color = group)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 1.6) +
  annotate("rect", xmin = op_start, xmax = op_end, ymin = -Inf, ymax = Inf,
           alpha = 0.1, fill = "red") +
  labs(
    title = "Monthly Street Violent Crime (Combined PROS Block Groups)",
    subtitle = "Two aggregated rank groups",
    x = "Month",
    y = "Incidents",
    color = ""
  ) +
  theme_minimal()

