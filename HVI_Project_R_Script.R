##########################################################
# Heat Vulnerability and Cooling Center Accessibility
# Eugene / Springfield, Oregon
# Daniel Friese
##########################################################

# This script builds a heat vulnerability index (HVI) at the census block group
# level, identifies statistically significant clusters using Getis-Ord Gi* hot
# spot analysis, and assesses walking distance access to cooling centers.


# ---- PACKAGES ----

library(tidycensus)  # ACS/Census data
library(tigris)      # Boundaries
library(sf)          # Vector data
library(terra)       # Raster data
library(dplyr)       # Data manipulation
library(spdep)       # Spatial analysis
library(tmap)        # Thematic mapping

options(tigris_use_cache = TRUE)

### census_api_key("............", install = TRUE)

crs <- 26910  # NAD83 UTM Zone 10N

##############################################
# PART 1: ACS DATA AND STUDY AREA
##############################################

# My original project used the ArcGIS Enrich tool to pull ACS variables.
# Here, I get them from the Census API via tidycensus.
#
# Poverty, disability and vehicle access are not available at the block group level
# from the Census API, so those were pulled at the tract level instead.

acs_bg_vars <- c(
  
  total_pop = "B01001_001",
  
  # Age 65+ broken into 12 age/sex categories
  age65_m_6566 = "B01001_020",
  age65_m_6769 = "B01001_021",
  age65_m_7074 = "B01001_022",
  age65_m_7579 = "B01001_023",
  age65_m_8084 = "B01001_024",
  age65_m_85plus = "B01001_025",
  age65_f_6566 = "B01001_044",
  age65_f_6769 = "B01001_045",
  age65_f_7074 = "B01001_046",
  age65_f_7579 = "B01001_047",
  age65_f_8084 = "B01001_048",
  age65_f_85plus = "B01001_049",
  
  tenure_denom = "B25003_001",
  tenure_renter = "B25003_003")

acs_tract_vars <- c(
  
  # Poverty
  poverty_denom = "B17001_001",
  poverty_below = "B17001_002",
  
  # Disability - summing individual age/sex subcategories
  disab_denom = "B18101_001",
  disab_m_u5 = "B18101_004",
  disab_m_5_17 = "B18101_007",
  disab_m_18_34 = "B18101_010",
  disab_m_35_64 = "B18101_013",
  disab_m_65_74 = "B18101_016",
  disab_m_75plus = "B18101_019",
  disab_f_u5 = "B18101_023",
  disab_f_5_17 = "B18101_026",
  disab_f_18_34 = "B18101_029",
  disab_f_35_64 = "B18101_032",
  disab_f_65_74 = "B18101_035",
  disab_f_75plus = "B18101_038",
  
  # Vehicle availability
  vehicle_denom = "B08201_001",
  vehicle_none = "B08201_002")

# Download block group data with geometry
lane_bg <- get_acs(
  geography = "block group",
  variables = acs_bg_vars,
  state = "OR",
  county = "Lane",
  year = 2024,
  geometry = TRUE,
  output = "wide")

# Tract-level data for variables not available at block group
lane_tracts <- get_acs(
  geography = "tract",
  variables = acs_tract_vars,
  state = "OR",
  county = "Lane",
  year = 2024,
  output = "wide")

lane_bg <- st_transform(lane_bg, crs = crs)

# Calculate percentages
lane_bg <- lane_bg %>%
  mutate(
    age65_sum = age65_m_6566E + age65_m_6769E + age65_m_7074E +
      age65_m_7579E + age65_m_8084E + age65_m_85plusE +
      age65_f_6566E + age65_f_6769E + age65_f_7074E +
      age65_f_7579E + age65_f_8084E + age65_f_85plusE,
    pct_age65 = ifelse(total_popE > 0, age65_sum / total_popE * 100, NA),
    pct_renter = ifelse(tenure_denomE > 0, tenure_renterE / tenure_denomE * 100, NA))

lane_tracts <- lane_tracts %>%
  mutate(
    pct_poverty = ifelse(poverty_denomE > 0, poverty_belowE / poverty_denomE * 100, NA),
    disab_sum = disab_m_u5E + disab_m_5_17E + disab_m_18_34E +
      disab_m_35_64E + disab_m_65_74E + disab_m_75plusE +
      disab_f_u5E + disab_f_5_17E + disab_f_18_34E +
      disab_f_35_64E + disab_f_65_74E + disab_f_75plusE,
    pct_disab = ifelse(disab_denomE > 0, disab_sum / disab_denomE * 100, NA),
    pct_no_veh = ifelse(vehicle_denomE > 0, vehicle_noneE / vehicle_denomE * 100, NA)) %>%
  select(tract_GEOID = GEOID, pct_poverty, pct_disab, pct_no_veh)

# Join tract percentages to block groups using the first 11 characters of GEOID
lane_bg <- lane_bg %>%
  mutate(tract_GEOID = substr(GEOID, 1, 11)) %>%
  left_join(lane_tracts, by = "tract_GEOID")

lane_bg <- lane_bg %>%
  select(GEOID, NAME, geometry,
         total_pop = total_popE,
         pct_age65, pct_no_veh, pct_renter,
         pct_poverty, pct_disab) %>%
  filter(!is.na(pct_age65), total_pop > 0)

# Load and merge Eugene and Springfield urban growth boundaries
eugene_ugb <- st_read("HVI/Eugene_UGB/Eugene_Urban_Growth_Boundary_(UGB)_-_HUB.shp") %>%
  st_transform(crs = crs)
springfield_ugb <- st_read("HVI/Springfield_UGB/Urban_Growth_Boundary.shp") %>%
  st_transform(crs = crs)

study_area <- st_union(eugene_ugb, springfield_ugb)

study_bg <- lane_bg[lengths(st_intersects(lane_bg, study_area)) > 0, ]


########################################################
# PART 2: RASTER DATA (LST AND TREE CANOPY)
########################################################

# Land surface temperature (LST)
# Derived from Landsat 9 thermal band (ST_B10) in Google Earth Engine,
# converted to Celsius and exported as a GeoTIFF
lst <- rast("HVI/EugSpr_L9_LST_20250712_UTM10N.tif")
lst <- project(lst, paste0("EPSG:", crs))

# NLCD Tree Canopy Cover
canopy <- rast("HVI/TreeCanopyCover_NLCD.tif")
canopy <- project(canopy, paste0("EPSG:", crs))


########################################################
# PART 3: ZONAL STATISTICS
########################################################

# Calculate mean LST and mean tree canopy cover for each block group
bg_vect <- vect(study_bg)

lst_means <- zonal(lst, bg_vect, fun = "mean", na.rm = TRUE)
canopy_means <- zonal(canopy, bg_vect, fun = "mean", na.rm = TRUE)

study_bg <- study_bg %>%
  mutate(
    mean_lst = lst_means[[1]],
    mean_canopy = canopy_means[[1]])


########################################################
# PART 4: BUILD THE HEAT VULNERABILITY INDEX (HVI)
########################################################

# LST is in Celsius and needs to be rescaled to 0-100 to match the ACS percent
# variables before combining them into the index.
rescale_minmax <- function(x) {
  rng <- range(x, na.rm = TRUE)
  (x - rng[1]) / (rng[2] - rng[1]) * 100}

study_bg <- study_bg %>%
  mutate(
    lst_scaled = rescale_minmax(mean_lst),
    canopy_inv = 100 - mean_canopy)  # Tree canopy is a protective factor, so it's inverted

# Weights based on Kohon et al. (2024) and Voelkel et al. (2018) - socioeconomic
# factors weighted more heavily than environmental.
w_age65 <- 20
w_poverty <- 20
w_disab <- 15
w_no_veh <- 15
w_renter <- 10
w_lst <- 12
w_canopy <- 8

study_bg <- study_bg %>%
  mutate(
    hvi = (pct_age65   * w_age65 +
             pct_poverty * w_poverty +
             pct_disab   * w_disab +
             pct_no_veh  * w_no_veh +
             pct_renter  * w_renter +
             lst_scaled  * w_lst +
             canopy_inv  * w_canopy) / 100)


########################################################
# PART 5: GETIS-ORD GI* HOT SPOT ANALYSIS
########################################################

# A fixed distance band of 1250 meters is used. In the original project, this
# value was chosen based on incremental spatial autocorrelation results, which
# showed the first significant clustering peak at that distance.

coords <- st_coordinates(st_centroid(study_bg))

nb <- dnearneigh(coords, d1 = 0, d2 = 1250)
nb_self <- include.self(nb)  # Gi* requires each feature to be its own neighbor

lw <- nb2listw(nb_self, style = "B", zero.policy = TRUE)

gi_result <- localG(study_bg$hvi, lw, zero.policy = TRUE)

study_bg <- study_bg %>%
  mutate(
    gi_z = as.numeric(gi_result),
    
    gi_class = case_when(
      gi_z >  2.58 ~ "Hot Spot 99%",
      gi_z >  1.96 ~ "Hot Spot 95%",
      gi_z >  1.65 ~ "Hot Spot 90%",
      gi_z < -2.58 ~ "Cold Spot 99%",
      gi_z < -1.96 ~ "Cold Spot 95%",
      gi_z < -1.65 ~ "Cold Spot 90%",
      TRUE          ~ "Not Significant"),
    
    # Ordered factor for a logical legend
    gi_class = factor(gi_class, levels = c(
      "Hot Spot 99%", "Hot Spot 95%", "Hot Spot 90%", "Not Significant",
      "Cold Spot 90%", "Cold Spot 95%", "Cold Spot 99%")))


########################################################
# PART 6: COOLING CENTER ACCESS
########################################################

cooling_centers <- read.csv("HVI/CoolingCenters.csv")
cooling_centers <- cooling_centers %>%
  select(NAME = Place_addr, X, Y) %>%
  st_as_sf(coords = c("X", "Y"), crs = 4326) %>%
  st_transform(crs = crs)

# Euclidean walking distance buffers at 3.1 mph (~83 m/min)
# Straight-line buffers will overestimate the reachable area compared to
# the network analysis used in the original project
buf_10 <- st_union(st_buffer(cooling_centers, 830))
buf_15 <- st_union(st_buffer(cooling_centers, 1245))
buf_20 <- st_union(st_buffer(cooling_centers, 1660))

hot_spots <- study_bg %>%
  filter(grepl("Hot Spot", gi_class))

hs_centroids <- st_centroid(hot_spots)

# Classify each hot spot by walking distance to nearest cooling center
hot_spots <- hot_spots %>%
  mutate(
    in_10 = lengths(st_intersects(hs_centroids, buf_10)) > 0,
    in_15 = lengths(st_intersects(hs_centroids, buf_15)) > 0,
    in_20 = lengths(st_intersects(hs_centroids, buf_20)) > 0,
    
    access = case_when(
      in_10          ~ "Within 10 min",
      in_15 & !in_10 ~ "10-15 min",
      in_20 & !in_15 ~ "15-20 min",
      TRUE           ~ "Over 20 min"))

print(table(hot_spots$access))


########################################################
# PART 7: MAPS
########################################################

tmap_mode("plot")

# Mask polygon for area outside the UGB
bbox_poly <- st_as_sfc(st_bbox(study_bg))
outside <- st_difference(bbox_poly, st_union(study_area))

# ---- Figure 2: HVI by block group -------------------------------------------

fig2 <- tm_shape(study_bg, bbox = st_bbox(study_area)) +
  tm_polygons(
    fill = "hvi",
    fill.scale = tm_scale_intervals(style = "quantile", n = 5, values = "YlOrRd"),
    fill.legend = tm_legend(title = "HVI Score"),
    col = "white",
    col_alpha = 0.4) +
  tm_shape(outside) +
  tm_fill(fill = "white") +
  tm_shape(study_area) +
  tm_borders(lwd = 1) +
  tm_title("Heat Vulnerability Index\nEugene / Springfield, OR",
           position = c(0.3, 0.97)) +
  tm_legend(position = c("right", "top")) +
  tm_compass(position = c("right", "bottom")) +
  tm_scalebar(position = c("left", "bottom"))

# ---- Figure 3: Gi* hot spot analysis ----------------------------------------

spot_colors <- c(
  "Hot Spot 99%"    = "#d73027",
  "Hot Spot 95%"    = "#f46d43",
  "Hot Spot 90%"    = "#fdae61",
  "Not Significant" = "#eeeeee",
  "Cold Spot 90%"   = "#abd9e9",
  "Cold Spot 95%"   = "#4575b4",
  "Cold Spot 99%"   = "#313695")

fig3 <- tm_shape(study_bg, bbox = st_bbox(study_area)) +
  tm_polygons(
    fill = "gi_class",
    fill.scale = tm_scale_categorical(values = spot_colors),
    fill.legend = tm_legend(title = "Confidence Level"),
    col = "white",
    col_alpha = 0.4) +
  tm_shape(outside) +
  tm_fill(fill = "white") +
  tm_shape(study_area) +
  tm_borders(lwd = 1) +
  tm_title("Getis-Ord Gi* Hot Spot Analysis\nEugene / Springfield, OR",
           position = c(0.3, 0.97), size = 1.1) +
  tm_legend(position = c("right", "top")) +
  tm_compass(position = c("right", "bottom")) +
  tm_scalebar(position = c("left", "bottom"))

# ---- Figure 4: Walking distance to cooling centers --------------------------

fig4 <- tm_shape(study_bg, bbox = st_bbox(study_area)) +
  tm_polygons(fill = "grey92", col = "grey70", col_alpha = 0.5) +
  tm_shape(buf_20) + tm_fill(fill = "lightgreen", fill_alpha = 0.6) +
  tm_shape(buf_15) + tm_fill(fill = "turquoise", fill_alpha = 0.6) +
  tm_shape(buf_10) + tm_fill(fill = "blue", fill_alpha = 0.6) +
  tm_shape(outside) +
  tm_fill(fill = "white") +
  tm_shape(study_area) +
  tm_borders(lwd = 1) +
  tm_shape(cooling_centers) +
  tm_dots(fill = "red", size = 0.6, shape = 16) +
  tm_add_legend(
    type = "polygons",
    fill = c("blue", "turquoise", "lightgreen"),
    labels = c("Within 10 min walk", "10-15 min walk", "15-20 min walk"),
    title = "Walking Distance\n(~3.1 mph)") +
  tm_add_legend(
    type = "symbols",
    fill = "red",
    shape = 16,
    labels = "Cooling Center") +
  tm_title("Estimated Walking Distance to Cooling Centers\nEugene / Springfield, OR",
           position = c(0.3, 0.97), size = 1) +
  tm_legend(position = c(.75, .9), frame = FALSE) +
  tm_compass(position = c("right", "bottom")) +
  tm_scalebar(position = c("left", "bottom"))

# ---- Figure 5: Hot spots beyond a 10-minute walk ----------------------------

beyond_10 <- hot_spots %>% filter(!in_10)

fig5 <- tm_shape(study_bg, bbox = st_bbox(study_area)) +
  tm_polygons(fill = "grey90", col = "grey70", col_alpha = 0.5) +
  tm_shape(buf_10) +
  tm_fill(fill = "blue", fill_alpha = 0.4) +
  tm_shape(beyond_10) +
  tm_polygons(fill = "red", col = "darkred", col_alpha = 0.85) +
  tm_shape(outside) +
  tm_fill(fill = "white") +
  tm_shape(study_area) +
  tm_borders(lwd = 1) +
  tm_shape(cooling_centers) +
  tm_dots(fill = "black", size = 0.4, shape = 17) +
  tm_add_legend(
    type = "polygons",
    fill = c("red", "blue"),
    labels = c("Hot spot > 10 min from cooling center", "10-min walk area")) +
  tm_add_legend(
    type = "symbols",
    fill = "black",
    shape = 17,
    labels = "Cooling Center") +
  tm_title("Vulnerability Hot Spots > 10-Min Walk\nfrom a Cooling Center",
           position = c(0.3, 0.97)) +
  tm_legend(position = c("right", "bottom")) +
  tm_compass(position = c("right", "top")) +
  tm_scalebar(position = c("left", "bottom"))

# ---- Save all figures -------------------------------------------------------

tmap_save(fig2, "output/fig2_hvi.png",            width = 8, height = 7, dpi = 200)
tmap_save(fig3, "output/fig3_hotspots.png",       width = 8, height = 7, dpi = 200)
tmap_save(fig4, "output/fig4_cooling_access.png", width = 8, height = 7, dpi = 200)
tmap_save(fig5, "output/fig5_beyond_10min.png",   width = 8, height = 7, dpi = 200)



