##########################################################
# Healthcare Deserts in Oregon
# Daniel Friese
##########################################################

# This script identifies healthcare deserts at the census block group level
# using two definitions: straight-line distance to the nearest acute care
# hospital, and a population-weighted composite index. Spatial autocorrelation
# (Global and Local Moran's I) is used to identify clustering patterns.

# ---- PACKAGES ----

library(sf)          # Vector data
library(ggplot2)     # Plotting
library(units)       # Unit-aware distance calculations
library(tidycensus)  # ACS data
library(tidyverse)   # Data wrangling
library(tigris)      # Oregon state boundary for map overlay
library(spdep)       # Spatial weights and Moran's I
library(patchwork)   # Combining multiple ggplots

crs <- 2992  # Oregon Lambert

### census_api_key("............", install = TRUE)

########################################################
# PART 1: HOSPITAL AND BLOCK GROUP DATA
########################################################

# Hospital point data from the Oregon Health Authority
hospitals_sf <- st_read("Hospital_Access\AcuteCareHospitals.gpkg") %>%
  st_transform(crs = crs)

# Median income and population at the block group level from ACS
or_bg <- get_acs(
  geography = "block group",
  variables = c(median_income = "B19013_001", population = "B01003_001"),
  state     = "OR",
  geometry  = TRUE,
  output    = "wide") %>%
  st_transform(crs = crs) %>%
  select(GEOID, median_income = median_incomeE, population = populationE, geometry)

########################################################
# PART 2: DISTANCE TO NEAREST HOSPITAL
########################################################

# Centroids are used as representative points for each block groups
or_bg_centroids <- st_centroid(or_bg)

or_bg$min_dist_km <- as.numeric(
  st_distance(
    or_bg_centroids,
    hospitals_sf[st_nearest_feature(or_bg_centroids, hospitals_sf), ],
    by_element = TRUE)) / 1000

# Definition 1: Block groups whose centroid is more than 50 km from any hospital
healthcare_deserts <- or_bg %>%
  filter(min_dist_km > 50)

########################################################
# PART 3: POPULATION-WEIGHTED DISTANCE
########################################################

# Straight distance alone flags remote but sparsely populated areas that affect
# very few people. Multiplying distance by population density creates a composite
# index that prioritizes block groups where more people lack nearby access.
or_bg <- or_bg %>%
  mutate(
    area_km2          = as.numeric(set_units(st_area(.), "km^2")),
    pop_density       = population / area_km2,
    pop_weighted_dist = min_dist_km * pop_density)

# Definition 2: Top 10% of block groups by population-weighted distance
healthcare_deserts_weighted <- or_bg %>%
  filter(pop_weighted_dist > quantile(pop_weighted_dist, 0.9, na.rm = TRUE))

########################################################
# PART 4: LIFE EXPECTANCY DATA
########################################################

# Life expectancy from the OHA is at the census tract level (coarser than block
# groups). The GEOID arrives as numeric and lowercase, so it needs to be renamed
# and coerced to character to match the tidycensus format.
life_expectancy <- read_csv("OR_life_expectancy.csv") %>%
  rename(GEOID = geoid) %>%
  mutate(GEOID = as.character(GEOID))

# sf objects support left_join() directly, so geometry is preserved without
# stripping and reattaching it.
tracts_with_life <- get_acs(
  geography = "tract",
  variables = "B01003_001",
  state     = "OR",
  geometry  = TRUE) %>%
  st_transform(crs = crs) %>%
  left_join(life_expectancy, by = "GEOID")


########################################################
# PART 5: SOCIOECONOMIC CORRELATIONS
########################################################

# Pearson correlation between median income and distance (unweighted and weighted).
# use = "complete.obs" excludes block groups with missing income data.
cor_income          <- cor(or_bg$min_dist_km,       or_bg$median_income, use = "complete.obs")
cor_income_weighted <- cor(or_bg$pop_weighted_dist, or_bg$median_income, use = "complete.obs")


########################################################
# PART 6: SPATIAL AUTOCORRELATION
########################################################

# Queen contiguity neighbors with row-standardized weights (style = "W").
# zero.policy = TRUE allows isolated polygons with no neighbors to remain in
# the dataset rather than throwing an error.

## Unweighted (distance-based definition)

healthcare_deserts <- healthcare_deserts %>%
  filter(!st_is_empty(geometry))

nb_hd  <- poly2nb(healthcare_deserts)
lw_hd  <- nb2listw(nb_hd, style = "W", zero.policy = TRUE)

global_moran <- moran.test(healthcare_deserts$min_dist_km, lw_hd, zero.policy = TRUE)

local_moran  <- localmoran(healthcare_deserts$min_dist_km, lw_hd, zero.policy = TRUE)

healthcare_deserts <- healthcare_deserts %>%
  mutate(
    local_I = local_moran[, "Ii"],
    local_p = local_moran[, "Pr(z != E(Ii))"])

# Oregon state boundary used as an outline in the local Moran's I maps
oregon_outline <- states(cb = TRUE) %>%
  filter(STUSPS == "OR") %>%
  st_transform(crs = st_crs(or_bg))

# Population-weighted definition

healthcare_deserts_weighted <- healthcare_deserts_weighted %>%
  filter(!st_is_empty(geometry))

nb_hd_w <- poly2nb(healthcare_deserts_weighted)
lw_hd_w <- nb2listw(nb_hd_w, style = "W", zero.policy = TRUE)

global_moran_weighted <- moran.test(healthcare_deserts_weighted$pop_weighted_dist,
                                    lw_hd_w, zero.policy = TRUE)

local_moran_weighted  <- localmoran(healthcare_deserts_weighted$pop_weighted_dist,
                                    lw_hd_w, zero.policy = TRUE)

healthcare_deserts_weighted <- healthcare_deserts_weighted %>%
  mutate(
    local_I = local_moran_weighted[, "Ii"],
    local_p = local_moran_weighted[, "Pr(z != E(Ii))"])

########################################################
# PART 7: MAPS
########################################################

#### Figure 1: Healthcare Deserts (Unweighted)
fig1 <- ggplot() +
  geom_sf(data = or_bg, fill = "grey90", color = "white", size = 0.2) +
  geom_sf(data = healthcare_deserts, fill = "red", color = "darkred",
          alpha = 0.5, size = 0.3) +
  geom_sf(data = hospitals_sf, color = "blue", shape = 21,
          fill = "blue", size = 2) +
  labs(title = "Healthcare Deserts in Oregon",
       subtitle = "Census block groups with centroids > 50 km from the nearest hospital",
       caption = "Data: Oregon Health Authority, American Community Survey") +
  theme_minimal() +
  theme(
    plot.title    = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(face = "italic", size = 9),
    panel.grid    = element_blank(),
    axis.text     = element_blank(),
    axis.ticks    = element_blank(),
    axis.title    = element_blank(),
    plot.caption  = element_text(size = 8, hjust = 0.5))

#### Figure 2: Life Expectancy by Census Tract
fig2 <- ggplot(tracts_with_life) +
  geom_sf(aes(fill = life_expectancy), color = NA) +
  scale_fill_viridis_c(option = "plasma", name = "Life Expectancy") +
  labs(title = "Oregon Life Expectancy by Census Tract",
       caption = "Data: Oregon Health Authority, American Community Survey") +
  theme_minimal() +
  theme(
    plot.title   = element_text(face = "bold", size = 16),
    panel.grid   = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    axis.title   = element_blank(),
    plot.caption = element_text(size = 8, hjust = 0.5))

#### Figure 3: Population Density
fig3 <- ggplot() +
  geom_sf(data = or_bg, aes(fill = pop_density), color = NA) +
  scale_fill_viridis_c(
    option = "viridis",
    trans  = "log10",
    name   = "Pop Density\n(people/km²)",
    breaks = c(1, 10, 100, 1000),
    labels = c("1", "10", "100", "1000")) +
  labs(title = "Population Density in Oregon",
       caption = "Data: American Community Survey") +
  theme_minimal() +
  theme(
    plot.title   = element_text(face = "bold", size = 16),
    panel.grid   = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    axis.title   = element_blank(),
    plot.caption = element_text(size = 8, hjust = 0.5))

#### Figure 4: Healthcare Deserts (Population-Weighted)
fig4 <- ggplot() +
  geom_sf(data = or_bg, fill = "grey90", color = "white", size = 0.2) +
  geom_sf(data = healthcare_deserts_weighted, fill = "red", color = "darkred",
          alpha = 0.5, size = 0.3) +
  geom_sf(data = hospitals_sf, color = "blue", shape = 21,
          fill = "blue", size = 1) +
  labs(title = "Healthcare Deserts in Oregon (Population-Weighted)",
       subtitle = "Block groups with high population-weighted distance to the nearest hospital",
       caption = "Data: Oregon Health Authority, American Community Survey") +
  theme_minimal() +
  theme(
    plot.title    = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(face = "italic", size = 10),
    panel.grid    = element_blank(),
    axis.text     = element_blank(),
    axis.ticks    = element_blank(),
    axis.title    = element_blank(),
    plot.caption  = element_text(size = 8, hjust = 0.5))

#### Figure 5: Portland Area (Population-Weighted
city_labels <- data.frame(
  city = c("Portland", "Beaverton", "Hillsboro"),
  x    = c(780000, 740000, 705000),
  y    = c(1380000, 1360000, 1375000))

fig5 <- ggplot() +
  geom_sf(data = or_bg, fill = "grey90", color = "white", size = 0.2) +
  geom_sf(data = healthcare_deserts_weighted, fill = "red", color = "darkred",
          alpha = 0.5, size = 0.3) +
  geom_sf(data = hospitals_sf, color = "blue", shape = 21,
          fill = "blue", size = 1.5) +
  geom_text(data = city_labels, aes(x = x, y = y, label = city),
            size = 3, color = "gray40", fontface = "bold") +
  coord_sf(xlim = c(690000, 840000), ylim = c(1300000, 1500000),
           expand = FALSE, clip = "off") +
  labs(title = "Healthcare Deserts — Portland Area",
       subtitle = "Census block groups with high population-weighted travel distance",
       caption = "Data: Oregon Health Authority, American Community Survey") +
  theme_minimal() +
  theme(
    plot.title    = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(face = "italic", size = 10),
    panel.grid    = element_blank(),
    axis.text     = element_blank(),
    axis.ticks    = element_blank(),
    axis.title    = element_blank(),
    plot.caption  = element_text(size = 8, hjust = 0.5))

#### Figure 6: Local Moran's I (Unweighted

fig6 <- ggplot(healthcare_deserts) +
  geom_sf(aes(fill = local_I), color = NA) +
  geom_sf(data = oregon_outline, fill = NA, color = "black", size = 0.5) +
  scale_fill_gradient2(midpoint = 0, low = "blue", mid = "white", high = "red",
                       name = "Local Moran's I") +
  labs(title = "Local Moran's I — Distance to Nearest Hospital",
       subtitle = "Identifying clusters in healthcare accessibility") +
  theme_minimal() +
  theme(
    plot.title    = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(face = "italic", size = 10),
    panel.grid    = element_blank(),
    axis.text     = element_blank(),
    axis.ticks    = element_blank(),
    axis.title    = element_blank())

#### Figure 7: Local Moran's I (Population-Weighted, Statewide

fig7 <- ggplot(healthcare_deserts_weighted) +
  geom_sf(aes(fill = local_I), color = NA) +
  geom_sf(data = oregon_outline, fill = NA, color = "black", size = 0.5) +
  scale_fill_gradient2(midpoint = 0, low = "blue", mid = "white", high = "red",
                       name = "Local Moran's I") +
  labs(title = "Local Moran's I — Population-Weighted Distance",
       subtitle = "Identifying clusters in healthcare accessibility") +
  theme_minimal() +
  theme(
    plot.title    = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(face = "italic", size = 10),
    panel.grid    = element_blank(),
    axis.text     = element_blank(),
    axis.ticks    = element_blank(),
    axis.title    = element_blank())

#### Figure 8: Local Moran's I (Population-Weighted, Portland)

fig8 <- ggplot(healthcare_deserts_weighted) +
  geom_sf(aes(fill = local_I), color = NA) +
  geom_sf(data = oregon_outline, fill = NA, color = "black", size = 0.5) +
  geom_text(data = city_labels, aes(x = x, y = y, label = city),
            size = 3, color = "gray40", fontface = "bold") +
  coord_sf(xlim = c(690000, 840000), ylim = c(1300000, 1500000),
           expand = FALSE, clip = "off") +
  scale_fill_gradient2(midpoint = 0, low = "blue", mid = "white", high = "red",
                       name = "Local Moran's I") +
  labs(title = "Local Moran's I — Portland Area",
       subtitle = "Identifying clusters in healthcare accessibility") +
  theme_minimal() +
  theme(
    plot.title    = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(face = "italic", size = 10),
    panel.grid    = element_blank(),
    axis.text     = element_blank(),
    axis.ticks    = element_blank(),
    axis.title    = element_blank())

#### Scatterplots: Income vs. Distance

fig9 <- ggplot(or_bg, aes(x = median_income, y = min_dist_km)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(title = "Median Income vs. Distance to Nearest Hospital",
       x = "Median Household Income",
       y = "Distance to Nearest Hospital (km)") +
  theme_minimal() +
  theme(plot.title = element_text(size = 10))

fig10 <- ggplot(or_bg, aes(x = median_income, y = pop_weighted_dist)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(title = "Median Income vs. Population-Weighted Distance",
       x = "Median Household Income",
       y = "Population-Weighted Distance") +
  theme_minimal() +
  theme(plot.title = element_text(size = 10))

#### Print all figures

print(fig1)
print(fig2)
fig1 + fig2
print(fig3)
fig1 + fig3
print(fig4)
print(fig5)
print(fig6)
print(fig7)
print(fig8)
print(global_moran)
print(global_moran_weighted)
print(fig9); cor_income
print(fig10); cor_income_weighted
       
       