#-------------------------------------------------
# 00 — Load data + make sf versions (EPSG:2263)
#-------------------------------------------------

library(here)
library(tidyverse)
library(janitor)
library(lubridate)
library(sf)
source(here("lib","definitions.R"))

#-------------------------------------------------
# List files
#-------------------------------------------------
list.files(here("data"))

#-------------------------------------------------
# LION 
#-------------------------------------------------
# Path to the geodatabase
lion_gdb <- here("data", "lion", "lion.gdb")

# List available layers inside the .gdb
st_layers(lion_gdb)

lion_gdb <- here("data", "lion", "lion.gdb")

lion <- st_read(lion_gdb, layer = "lion") %>%
  st_transform(2263) %>% 
  #filter(feature_typ %in% c('0','W')) %>%
  filter(!st_geometry_type(.) %in% c("MULTICURVE")) %>%
  clean_names() # NYC coordinates, ft
#-------------------------------------------------
# Spatial boundary files (already 2263; enforce)
#-------------------------------------------------
nycb <- st_read(here("data", "nycb2020_25d"), quiet = TRUE) |>
  st_transform(2263) |>
  clean_names()

nyct <- st_read(here("data", "nyct2020_25d"), quiet = TRUE) |>
  st_transform(2263) |>
  clean_names()

nynta <- st_read(here("data", "nynta2020_25d"), quiet = TRUE) |>
  st_transform(2263) |>
  clean_names()

nypp <- st_read(here("data", "nypp_25d"), quiet = TRUE) |>
  st_transform(2263) |>
  clean_names()

#-------------------------------------------------
# CSV reader
#-------------------------------------------------
read_clean_csv <- function(file) {
  readr::read_csv(here("data", file), show_col_types = FALSE) |>
    clean_names()
}

#-------------------------------------------------
# Load CSVs
#-------------------------------------------------
arrests <- bind_rows(
  read_clean_csv("NYPD_Arrests_Data_(Historic)_20251214.csv"),
  read_clean_csv("NYPD_Arrest_Data_(Year_to_Date)_20251214.csv")
)

directed_patrols <- bind_rows(
  read_clean_csv("NYPD_Calls_for_Service_(Historic)_20251214.csv") |>
    mutate(incident_time = as.character(incident_time)),
  read_clean_csv("NYPD_Calls_for_Service_(Year_to_Date)_20251214.csv") |>
    mutate(incident_time = as.character(incident_time))
)

complaints <- bind_rows(
  read_clean_csv("NYPD_Complaint_Data_Historic_20251214.csv") |>
    mutate(housing_psa = as.character(housing_psa)),
  read_clean_csv("NYPD_Complaint_Data_Current_(Year_To_Date)_20251214.csv") |>
    mutate(housing_psa = as.character(housing_psa))
) %>%
  filter(!pd_cd %in% c(111,113,114))

criminal_court_summons <- read_clean_csv("NYPD_Criminal_Court_Summons_(Historic)_20251214.csv")

oath_summons <- read_clean_csv("NYPD_OATH_Summons_Data_20251214.csv")

precincts <- read_clean_csv("nypd_precinct_locations.csv")

#-------------------------------------------------
# Make sf versions (EPSG:2263)
#-------------------------------------------------

# Arrests: authoritative NYPD X/Y (2263)
arrests_sf <- arrests |>
  filter(!is.na(x_coord_cd), !is.na(y_coord_cd)) |>
  st_as_sf(coords = c("x_coord_cd", "y_coord_cd"), crs = 2263, remove = FALSE)

# Calls for Service: authoritative X/Y (geo_cd_x, geo_cd_y) (2263)
directed_patrols_sf <- directed_patrols |>
  filter(!is.na(geo_cd_x), !is.na(geo_cd_y)) |>
  st_as_sf(coords = c("geo_cd_x", "geo_cd_y"), crs = 2263, remove = FALSE)

# Complaints: authoritative NYPD X/Y (2263)
complaints_sf <- complaints |>
  filter(!is.na(x_coord_cd), !is.na(y_coord_cd)) |>
  st_as_sf(coords = c("x_coord_cd", "y_coord_cd"), crs = 2263, remove = FALSE)

# Criminal Court Summonses: authoritative X/Y (2263)
criminal_court_summons_sf <- criminal_court_summons |>
  filter(!is.na(x_coordinate_cd), !is.na(y_coordinate_cd)) |>
  st_as_sf(coords = c("x_coordinate_cd", "y_coordinate_cd"), crs = 2263, remove = FALSE)

# OATH Summonses: X/Y are character → numeric → sf (2263)
oath_summons_sf <- oath_summons |>
  mutate(
    x_coord_cd = as.numeric(x_coord_cd),
    y_coord_cd = as.numeric(y_coord_cd)
  ) |>
  filter(!is.na(x_coord_cd), !is.na(y_coord_cd)) |>
  st_as_sf(coords = c("x_coord_cd", "y_coord_cd"), crs = 2263, remove = FALSE)

# Precincts: lat/long → sf (4326) → transform to 2263
precincts_sf <- precincts |>
  filter(!is.na(longitude), !is.na(latitude)) |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) |>
  st_transform(2263)

#-------------------------------------------------
# Quick CRS checks
#-------------------------------------------------
st_crs(arrests_sf)
st_crs(directed_patrols_sf)
st_crs(complaints_sf)
st_crs(criminal_court_summons_sf)
st_crs(oath_summons_sf)
st_crs(precincts_sf)



#create zone
# Filter out unwanted street types and intersect with pct110_115 as before
roads_filtered <- lion %>%
  filter(!str_detect(street, "DRIVEWAY|PED |PEDESTRIAN PATH")) %>%
  filter(l_boro == 4) 


# Separate the roads that contain "ROOSEVELT" in the street name
roosevelt_roads <- roads_filtered %>%
  filter(str_detect(street, "ROOSEVELT")) %>%
  filter(!str_detect(street, "BRIDGE"))
#see it
plot(st_geometry(roosevelt_roads))


#rosie buffer
rosie_buffer <- roosevelt_roads %>%
  st_buffer(dist = 500) %>% # need streets midblock off rosie
  st_union(.) %>%
  st_make_valid()

plot(st_geometry(rosie_buffer))


#inspect the buffer

library(leaflet)
# Transform CRS for leaflet (since leaflet uses EPSG 4326)
roads_filtered_4326 <- st_transform(roads_filtered, 4326)
rosie_buffer_4326 <- st_transform(rosie_buffer, 4326)

# Create leaflet map
leaflet() %>%
  addTiles() %>%  # Add the default OpenStreetMap tiles
  #addPolylines(data = roads_filtered_4326, color = "black", weight = 1, opacity = 0.7) %>%  # Add roads in black
  addPolygons(data = rosie_buffer_4326, color = "red", weight = 2, fillOpacity = 0.3) %>%  # Add buffer around Roosevelt Avenue in red
  setView(lng = -73.878, lat = 40.748, zoom = 13)  # Set a reasonable zoom level for the area


# Keep everything in 2263 for geometry ops
roads_2263 <- roads_filtered |>
  st_zm(drop = TRUE, what = "ZM") |>
  st_cast("MULTILINESTRING", warn = FALSE) |>
  st_make_valid()

buffer_2263 <- rosie_buffer |>
  st_zm(drop = TRUE, what = "ZM") |>
  st_make_valid()

# Use st_filter first (fast, avoids heavy intersection unless needed)
roads_within_buffer_2263 <- roads_2263 |>
  st_filter(buffer_2263, .predicate = st_intersects)

# If you truly need clipped geometries (not just “within”), then intersect
roads_within_buffer_2263 <- st_intersection(roads_within_buffer_2263, buffer_2263)

# Now transform for leaflet
roads_within_buffer_4326 <- st_transform(roads_within_buffer_2263, 4326)
rosie_buffer_4326 <- st_transform(buffer_2263, 4326)

# Get geometries for the westernmost and easternmost segments based on segment IDs
west_segment <- roosevelt_roads %>% filter(segment_id == '0076052')
east_segment <- roosevelt_roads %>% filter(segment_id == '0083526')

# Create a bounding box that includes both westernmost and easternmost segments
bounding_box <- st_union(west_segment, east_segment) %>% st_bbox()

# Crop the roosevelt_roads to include only segments within the bounding box
roosevelt_roads_limited <- st_crop(roosevelt_roads, bounding_box)

# Plot the filtered Roosevelt Avenue segments
ggplot() +
  geom_sf(data = roads_filtered, color = "lightgray") +  # Plot all roads in black
  geom_sf(data = roosevelt_roads_limited, color = "red", size = 1.5) +  # Plot limited Roosevelt roads in red, thicker
  theme_minimal()

# Create a buffer around the limited Roosevelt Avenue segments
rosie_buffer <- roosevelt_roads_limited %>%
  st_buffer(dist = 500) %>%
  st_union(.) %>%
  st_make_valid()

# Plot the buffer and limited Roosevelt road segments
# Plot the filtered Roosevelt Avenue segments
ggplot() +
  geom_sf(data = roads_filtered, color = "lightgray") +  # Plot all roads in black
  geom_sf(data = rosie_buffer, color = "red", size = 1.5) +  # Plot limited Roosevelt roads in red, thicker
  theme_minimal()


# Filter for LINESTRING geometries in roosevelt_roads_limited before transforming
roosevelt_roads_limited <- roosevelt_roads_limited %>%
  filter(st_geometry_type(.) == "LINESTRING")

# Transform to EPSG 4326 for leaflet
roosevelt_roads_limited_4326 <- st_transform(roosevelt_roads_limited, 4326)

# Create a leaflet map
leaflet() %>%
  addTiles() %>%  # Add default OpenStreetMap tiles
  addPolylines(
    data = roosevelt_roads_limited_4326,
    color = "red",
    weight = 2,
    opacity = 0.8,
    label = ~segment_id,  # Use segment_id as the label on mouseover
    highlightOptions = highlightOptions(color = "blue", weight = 3)
  ) %>%
  setView(lng = -73.878, lat = 40.748, zoom = 13)  # Adj

#new rosie buffer
rosie_buffer <- roosevelt_roads_limited %>%
  st_buffer(dist = 500) %>%
  st_union(.) %>%
  st_make_valid()


rosie_buffer_4326 <- rosie_buffer %>%
  st_transform(4326)


# Create a leaflet map
# Create leaflet map
leaflet() %>%
  addTiles() %>%  # Add the default OpenStreetMap tiles
  #addPolylines(data = roads_filtered_4326, color = "black", weight = 1, opacity = 0.7) %>%  # Add roads in black
  addPolygons(data = rosie_buffer_4326, color = "red", weight = 2, fillOpacity = 0.3) %>%  # Add buffer around Roosevelt Avenue in red
  setView(lng = -73.878, lat = 40.748, zoom = 13)  # Set a reasonable zoom level for the area

# Save the Rosie buffer as an .rds file
saveRDS(
  rosie_buffer,
  file = here("output", "rosie_buffer.rds")
)

# Save the LINESTRING layer as an .rds file
saveRDS(
  roosevelt_roads_limited,
  file = here("output", "roosevelt_roads_limited.rds")
)


bg_path <- here("output", "nyc_block_groups_2020_2263.rds")

if (file.exists(bg_path)) {
  
  nyc_bgs <- readRDS(bg_path)
  
} else {
  
  library(tigris)
  options(tigris_use_cache = TRUE)
  
  nyc_bgs <- block_groups(
    state = "NY",
    county = c("005", "047", "061", "081", "085"),
    year = 2020,
    cb = TRUE,
    class = "sf"
  ) |>
    clean_names() |>
    st_transform(2263)
  
  dir.create(here("output"), showWarnings = FALSE, recursive = TRUE)
  saveRDS(nyc_bgs, bg_path)
}


#combine oath and c's
#-------------------------------------------------
# Combine OATH + Criminal Court summonses (already loaded)
# Expects objects: oath_summons, criminal_court_summons
# Output: all_summonses (sf, EPSG:2263) + date field
#-------------------------------------------------

# Criminal Court summonses (already in common schema)
rosie_cs <- criminal_court_summons |>
  distinct(summons_key, .keep_all = TRUE) |>
  mutate(summons_type = "criminal court") |>
  dplyr::select(
    summons_key,
    summons_date,
    offense_description,
    age_group,
    sex,
    race,
    x_coordinate_cd,
    y_coordinate_cd,
    summons_type
  )

# OATH summonses (rename into common schema)
oath <- oath_summons |>
  mutate(summons_type = "OATH") |>
  rename(
    summons_key = evnt_key,
    summons_date = occur_date,
    offense_description = law_desc,
    x_coordinate_cd = x_coord_cd,
    y_coordinate_cd = y_coord_cd
  ) |>
  mutate(
    x_coordinate_cd = as.numeric(x_coordinate_cd),
    y_coordinate_cd = as.numeric(y_coordinate_cd)
  ) |>
  dplyr::select(
    summons_key,
    summons_date,
    offense_description,
    age_group,
    sex,
    race,
    x_coordinate_cd,
    y_coordinate_cd,
    summons_type
  )

# Combine + make sf + parse date
all_summonses <- bind_rows(oath, rosie_cs) |>
  filter(!is.na(x_coordinate_cd), !is.na(y_coordinate_cd)) |>
  st_as_sf(coords = c("x_coordinate_cd", "y_coordinate_cd"), crs = 2263, remove = FALSE) |>
  mutate(date = mdy(summons_date))


#create crime data:
violent_ky  <- c(101, 104, 105, 106, 344)
property_ky <- c(107, 109, 110)

outside_loc_keywords <- c(
  "FRONT OF","OPPOSITE OF","OUTSIDE","REAR OF",
  "STREET","IN STREET","SIDEWALK"
)

outside_prem_keywords <- c(
  "PARK","STREET","PUBLIC PLACE","HIGHWAY",
  "BRIDGE","SIDEWALK","VACANT LOT",
  "PUBLIC HOUSING AREA","OUTSIDE"
)

# Flag "street" (outdoor)
compl_street <- complaints_sf %>%
  mutate(
    date = mdy(rpt_dt),
    loc_of_occur_desc = str_to_upper(coalesce(loc_of_occur_desc, "")),
    prem_typ_desc     = str_to_upper(coalesce(prem_typ_desc, "")),
    is_outdoor =
      str_detect(loc_of_occur_desc, str_c(outside_loc_keywords, collapse = "|")) |
      str_detect(prem_typ_desc,     str_c(outside_prem_keywords, collapse = "|"))
  )

# 4 objects
violent_crime <- complaints_sf %>%
  filter(ky_cd %in% violent_ky) %>%
  mutate(date = mdy(rpt_dt))

property_crime <- complaints_sf %>%
  filter(ky_cd %in% property_ky) %>%
  mutate(date = mdy(rpt_dt))

street_violent_crime <- compl_street %>%
  filter(is_outdoor, jurisdiction_code != 1, ky_cd %in% violent_ky)

street_property_crime <- compl_street %>%
  filter(is_outdoor, jurisdiction_code != 1, ky_cd %in% property_ky)



# Core spatial objects
glimpse(lion)
glimpse(nycb)
glimpse(nyct)
glimpse(nynta)
glimpse(nypp)

# Raw tabular data
glimpse(arrests)
glimpse(directed_patrols)
glimpse(complaints)
glimpse(criminal_court_summons)
glimpse(oath_summons)
glimpse(precincts)

# sf versions
glimpse(arrests_sf)
glimpse(directed_patrols_sf)
glimpse(complaints_sf)
glimpse(criminal_court_summons_sf)
glimpse(oath_summons_sf)
glimpse(precincts_sf)

# Roads / zone objects
glimpse(roads_filtered)
glimpse(roosevelt_roads)
glimpse(roosevelt_roads_limited)
glimpse(rosie_buffer)

# Census geography
glimpse(nyc_bgs)

# Summonses combined
glimpse(rosie_cs)
glimpse(oath)
glimpse(all_summonses)

# Crime classification outputs
glimpse(violent_crime)
glimpse(property_crime)
glimpse(street_violent_crime)
glimpse(street_property_crime)


x<- violent_crime %>%
  st_drop_geometry() %>%
  count(date)
