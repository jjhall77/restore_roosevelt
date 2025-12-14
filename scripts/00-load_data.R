#-------------------------------------------------
# 00 — Load data + make sf versions (EPSG:2263)
#-------------------------------------------------

library(here)
library(tidyverse)
library(janitor)
library(lubridate)
library(sf)

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
)

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

plot(st_geometry(roosevelt_roads))


#rosie buffer
rosie_buffer <- roosevelt_roads %>%
  st_buffer(dist = 500) %>%
  st_union(.) %>%
  st_make_valid()

plot(st_geometry(rosie_buffer))


#mapping the buffer

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



# Transform CRS for leaflet (since leaflet uses EPSG 4326)
roads_filtered_4326 <- st_transform(roads_filtered, 4326)
rosie_buffer_4326 <- st_transform(rosie_buffer, 4326)


# Make geometries valid to prevent intersection errors
roads_filtered_4326 <- st_make_valid(roads_filtered_4326)
rosie_buffer_4326 <- st_make_valid(rosie_buffer_4326)

# Perform intersection to get only road segments within rosie_buffer
roads_within_buffer <- st_intersection(roads_filtered_4326, rosie_buffer_4326)

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
saveRDS(rosie_buffer, file = "../output/rosie_buffer.rds")
# Save the LINESTRING layer (e.g., roosevelt_roads_limited) as an .rds file
saveRDS(roosevelt_roads_limited, "../output/roosevelt_roads_limited.rds")


#load spatial data
library(tigris); options(tigris_use_cache = TRUE)
nyc_bgs <- block_groups(state = "NY",
                        county = c("005","047","061","081","085"),
                        year = 2020, cb = TRUE, class = "sf") |> 
  clean_names() %>%
  st_transform(2263)   # NYC’s projected CRS


#load pct data
nypp <- read_sf('../data/nypp_25b') %>%
  clean_names() %>%
  st_set_crs(2263)

#load cb data
nycb <- read_sf('../data/nycb2020_25a') %>%
  clean_names() %>%
  st_set_crs(2263) 


nyct <- read_sf('../data/nyct2020_25b')%>%
  clean_names() %>%
  st_set_crs(2263) 





#open records-----------

list.files('../data/')
#load crimes------
x1 <- read_csv('../data/NYPD_Complaint_Data_Current__Year_To_Date__20250722.csv') %>%
  clean_names() %>%
  select(-housing_psa)
x2 <- read_csv('../data/NYPD_Complaint_Data_Historic_20250622.csv') %>%
  clean_names() 
crimes <- bind_rows(x1,x2) %>%
  distinct(cmplnt_num, .keep_all = T) %>%
  mutate(date = mdy(rpt_dt),
         cmplnt_fr_dt = mdy(cmplnt_fr_dt),
         hour  = hour(hms(cmplnt_fr_tm))) %>%
  upDate()

crimes.sf <- crimes %>%
  filter(!is.na(x_coord_cd)) %>%
  st_as_sf(coords = c("x_coord_cd","y_coord_cd"), crs = 2263) 
  
#load arrests--------
list.files('../data/')
x1 <- read_csv('../data/NYPD_Arrest_Data__Year_to_Date__20250622.csv') %>%
  clean_names() %>%
  select(-housing_psa)
x2 <- read_csv('../data/NYPD_Arrests_Data__Historic__20250622.csv') %>%
  clean_names() 

arrests <- bind_rows(x1,x2) %>%
  distinct(arrest_key, .keep_all = T) %>%
  mutate(date = mdy(arrest_date)) %>%
  upDate()

arrests.sf <- arrests %>%
  st_as_sf(coords = c("x_coord_cd","y_coord_cd"), crs = 2263) 



#summonses
list.files('../data/')
#load c court summonses
ytd_cs <- read_csv('../data/NYPD_Criminal_Court_Summons_Incident_Level_Data__Year_To_Date__20250622.csv') %>%
  clean_names() %>%
  mutate(precinct_of_occur = as.double(precinct_of_occur))
hx_cs <- read_csv('../data/NYPD_Criminal_Court_Summons__Historic__20250622.csv') %>%
  clean_names()
rosie_cs <-  bind_rows(ytd_cs, hx_cs) %>%
  distinct(summons_key,.keep_all = T) %>%
  mutate(summons_type = "criminal court")

#load oath
oath <- read_csv('../data/NYPD_OATH_Summons_Data_20250622.csv') %>%
  clean_names() %>%
  #filter(city_nm %in% c("QUEENS")) %>%
  mutate(summons_type = "OATH")


#combine the two with common fields
rosie_cs <- rosie_cs %>%
  dplyr::select(summons_key, summons_date, offense_description, age_group, sex, race, x_coordinate_cd,
                y_coordinate_cd, summons_type)

oath <- oath %>%
  rename(summons_key = evnt_key,
         summons_date = occur_date,
         offense_description = law_desc,
         x_coordinate_cd = x_coord_cd,
         y_coordinate_cd = y_coord_cd) %>%
  dplyr::select(summons_key, summons_date, offense_description, age_group, sex, race, x_coordinate_cd,
                y_coordinate_cd, summons_type)

all_summonses <- bind_rows(oath, rosie_cs) %>%
  filter(!is.na(x_coordinate_cd)) %>%
  st_as_sf(coords = c("x_coordinate_cd","y_coordinate_cd"), crs = 2263) 

all_summonses <- all_summonses %>%
  mutate(date = mdy(summons_date)) %>%
  upDate()


#vacate orders
vacate_orders <- read_csv('../data/Order_to_Repair_Vacate_Orders_20250622.csv') %>%
  clean_names()
vacate.sf <- vacate_orders %>%
  mutate(date = mdy(vacate_effective_date)) %>%
  filter(!is.na(latitude)) %>%
  st_as_sf(coords = c("longitude","latitude")) %>%
  st_set_crs(4326) %>%
  st_transform(2263)


plot(st_geometry(vacate.sf), add = T)
plot(st_geometry(nypp), add = T, col = "red")


#sla
list.files('../data/')
sla <- read_csv('../data/Current_Liquor_Authority_Active_Licenses_20250622.csv') %>%
  clean_names()

sla.sf <- sla %>% 
  ## replace NA with a legal empty geometry token
  mutate(
    wkt = ifelse(is.na(georeference), "POINT EMPTY", georeference),
    geometry = st_as_sfc(wkt, crs = 4326)
  ) %>% 
  st_as_sf() %>%          # upgrades the data-frame to sf
  select(-wkt) %>%        # drop helper column
  st_make_valid() %>%     # belt-and-suspenders
  st_transform(2263)      # NYC / Long Island feet (optional)












#rosie 311
#https://data.cityofnewyork.us/Social-Services/311-Service-Requests-from-2010-to-Present/erm2-nwe9/explore/query/SELECT%0A%20%20%60unique_key%60%2C%0A%20%20%60created_date%60%2C%0A%20%20%60closed_date%60%2C%0A%20%20%60agency%60%2C%0A%20%20%60agency_name%60%2C%0A%20%20%60complaint_type%60%2C%0A%20%20%60descriptor%60%2C%0A%20%20%60location_type%60%2C%0A%20%20%60incident_zip%60%2C%0A%20%20%60incident_address%60%2C%0A%20%20%60street_name%60%2C%0A%20%20%60cross_street_1%60%2C%0A%20%20%60cross_street_2%60%2C%0A%20%20%60intersection_street_1%60%2C%0A%20%20%60intersection_street_2%60%2C%0A%20%20%60address_type%60%2C%0A%20%20%60city%60%2C%0A%20%20%60landmark%60%2C%0A%20%20%60facility_type%60%2C%0A%20%20%60status%60%2C%0A%20%20%60due_date%60%2C%0A%20%20%60resolution_description%60%2C%0A%20%20%60resolution_action_updated_date%60%2C%0A%20%20%60community_board%60%2C%0A%20%20%60bbl%60%2C%0A%20%20%60borough%60%2C%0A%20%20%60x_coordinate_state_plane%60%2C%0A%20%20%60y_coordinate_state_plane%60%2C%0A%20%20%60open_data_channel_type%60%2C%0A%20%20%60park_facility_name%60%2C%0A%20%20%60park_borough%60%2C%0A%20%20%60vehicle_type%60%2C%0A%20%20%60taxi_company_borough%60%2C%0A%20%20%60taxi_pick_up_location%60%2C%0A%20%20%60bridge_highway_name%60%2C%0A%20%20%60bridge_highway_direction%60%2C%0A%20%20%60road_ramp%60%2C%0A%20%20%60bridge_highway_segment%60%2C%0A%20%20%60latitude%60%2C%0A%20%20%60longitude%60%2C%0A%20%20%60location%60%0AWHERE%0A%20%20%28%60created_date%60%0A%20%20%20%20%20BETWEEN%20%222015-12-31T00%3A00%3A00%22%20%3A%3A%20floating_timestamp%0A%20%20%20%20%20AND%20%222024-10-13T13%3A02%3A58%22%20%3A%3A%20floating_timestamp%29%0A%20%20AND%20caseless_one_of%28%0A%20%20%20%20%60incident_zip%60%2C%0A%20%20%20%20%22%22%2C%0A%20%20%20%20%2211372%22%2C%0A%20%20%20%20%2211373%22%2C%0A%20%20%20%20%2211368%22%2C%0A%20%20%20%20%2211377%22%0A%20%20%29%0AORDER%20BY%20%60created_date%60%20DESC%20NULL%20FIRST/page/filter
calls311 <- read_csv('../data/rosie_311.csv') %>%
  clean_names()

calls311 <- calls311 %>%
  mutate(date = mdy_hm(created_date),
         week_date = ymd(cut(date, "weeks" )),
         month_date = ymd(cut(date, "months")),
         year = year(date),
         week_num = week(date),
         month = month(date),
         day = wday(date, label = T),
         hour = hour(date),
         yearday = yday(date))


#load and combine summonses
list.files('../data/')

#load c court summonses
ytd_cs <- read_csv('../data/NYPD_Criminal_Court_Summons_Incident_Level_Data__Year_To_Date__20241014.csv') %>%
  clean_names() %>%
  mutate(precinct_of_occur = as.double(precinct_of_occur))
hx_cs <- read_csv('../data/NYPD_Criminal_Court_Summons__Historic__20241014.csv') %>%
  clean_names()
rosie_cs <-  bind_rows(ytd_cs, hx_cs) %>%
  distinct(summons_key,.keep_all = T) %>%
  mutate(summons_type = "criminal court")


#load oath
oath <- read_csv('../data/NYPD_OATH_Summons_Data_20241014.csv') %>%
  clean_names() %>%
  filter(city_nm %in% c("QUEENS")) %>%
  mutate(summons_type = "OATH")


#combine the two with common fields
rosie_cs <- rosie_cs %>%
  dplyr::select(summons_key, summons_date, offense_description, age_group, sex, race, x_coordinate_cd,
                y_coordinate_cd, summons_type)

oath <- oath %>%
  rename(summons_key = evnt_key,
         summons_date = occur_date,
         offense_description = law_desc,
         x_coordinate_cd = x_coord_cd,
         y_coordinate_cd = y_coord_cd) %>%
  dplyr::select(summons_key, summons_date, offense_description, age_group, sex, race, x_coordinate_cd,
                y_coordinate_cd, summons_type)

all_summonses <- bind_rows(oath, rosie_cs) %>%
  filter(!is.na(x_coordinate_cd)) %>%
  st_as_sf(coords = c("x_coordinate_cd","y_coordinate_cd"), crs = 2263) %>%
  st_intersection(rosie_buffer)

all_summonses <- all_summonses %>%
  mutate(date = mdy(summons_date)) %>%
  upDate()




#arrests
list.files('../data/')

arrests_old <- read_csv('../data/NYPD_Arrests_Data__Historic__20241014.csv') %>%
  clean_names()
arrests_new <- read_csv('../data/NYPD_Arrest_Data__Year_to_Date__20241014 (1).csv') %>%
  clean_names()

arrests <- bind_rows(arrests_old,arrests_new) %>%
  distinct(arrest_key, .keep_all = T)

arrests.sf <- arrests %>%
  mutate(date = mdy(arrest_date)) %>%
  upDate() %>%
  st_as_sf(coords = c("x_coord_cd","y_coord_cd"), crs = 2263) %>%
  st_intersection(rosie_buffer)

#calls 911
list.files('../data/')
calls911 <- read_csv('../data/NYPD_Calls_for_Service__Year_to_Date__20241014.csv') %>%
  clean_names() %>%
  mutate(incident_time = as.character(incident_time))
old_calls911 <- read_csv('../data/NYPD_Calls_for_Service__Historic__20241014.csv') %>%
  clean_names()

calls911 <- bind_rows(calls911, old_calls911) 

calls911.sf <- calls911 %>%
  distinct(cad_evnt_id, .keep_all = T) %>%
  mutate(date = mdy(create_date)) %>%
  upDate() %>%
  st_as_sf(coords = c("geo_cd_x","geo_cd_y"), crs = 2263) %>%
  st_intersection(rosie_buffer) %>%
  mutate(disp_ts = mdy_hms(disp_ts),
         arrivd_ts = mdy_hms(arrivd_ts),
         disp_arrival = arrivd_ts - disp_ts,
         hour = hour(disp_ts))

