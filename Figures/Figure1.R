###### MAPS Project: Density Dependence #######
### Script name: Figure1.R
### Author(s): SLJ

########### Objective/Description of Script #####################
# create Figure 1, which contains two panels
# Panel 1:  map showing locations of MAPS banding stations
# Panel 2: predicted relationship between productivity and abundance and the two response variables
#################################################################

### Setup ###

# load packages
library(here)
library(dplyr)
library(ggplot2)
library(tidybayes)
library(bayesplot)
library(ggspatial)
library(sf)
library(terra)
library(colorspace)
library(rnaturalearth)
library(ggdist)
library(patchwork)

#################################################################
### Import data files ###

# import productivity-abundance data
DD_breed_dat <- readRDS(here("Outputs", "DD_breed_dat.rds"))

# import station information 
station_info <- read.csv(here("Data", "MAPS_STATION_location_and_operations.csv"), header=T) %>%
  mutate(STA = factor(STA)) # convert STA to a factor

# import data to make species abundance-productivity curve
epred_df_m2 <- readRDS(here("Outputs", "epred_df_m2.rds"))

#################################################################
### Create Panel A of Figure 1: Maps of MAPS banding stations ###

# use productivity-abundance to generate a list of MAPS stations that were in the final analysis
STAlist <- DD_breed_dat %>%
  select(STA) %>% distinct() 
nrow(STAlist) # 378 MAPS stations

# list of stations with number of species (out of 62 in the analysis) at each station
STAlist <- DD_breed_dat %>%
  select(STA, SPEC) %>% 
  distinct() %>%
  group_by(STA) %>%
  count()
nrow(STAlist) # 378 MAPS stations

# import station information to get geographic coordinates
STAlist_coords <- left_join(STAlist, station_info, by = "STA") %>%
  select(STA, DECLNG, DECLAT, STATE, n) %>%
  filter(DECLNG < 999) # remove any stations with missing coordinates (marked by 999 for longitude)
# one stations is missing coordinates

# make the station locations with their data into spatial points and set CRS as NAD83
STA_sf <- st_as_sf(STAlist_coords, coords = c(2:3), crs="+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs +towgs84=0,0,0")

# get map of USA and Canada from rnaturalearth package
USACanada_map <- ne_countries(
  country = c("Canada", "United States of America"), 
  returnclass = "sf"
)

mapcrs <- crs(USACanada_map) # save crs of the map layer

# plot the map layer
(mapplot <- ggplot(data = USACanada_map) +
    geom_sf(fill = "white", color = "black") +
    theme_minimal())


# convert station coordinates to use same crs as map layer
STA_sf_plot <- st_transform(STA_sf, crs = mapcrs)


# add stations points to the USA map 
MAPS_map <- mapplot + 
  layer_spatial(data = STA_sf_plot, color = "#484554", alpha = 0.5) +
  xlab("Longitude") + ylab("Latitude") +
  theme(axis.title.x = element_text(size = 12, margin = margin(t=5), family = "Arial"),
        axis.title.y = element_text(size = 12, margin = margin(r=5), family = "Arial"),
        axis.text = element_text(size=10, family = "Arial"))

MAPS_map

#################################################################
### Create Panel B of Figure 1: Species curve with extracted response variables ###

# get data for species with a Type II curve

SPEC_dat <- epred_df_m2 %>%
  filter(SPEC == "NOCA") # specify the species to keep

# create plot of single species with a Type II curve

TypeII <- ggplot(data = SPEC_dat, aes(x = Adult, y = .epred)) + 
  ggdist::stat_lineribbon(color = "#1874CD") + # change color here
  scale_fill_manual(values = colorspace::lighten("#1874CD", c(0.95, 0.75, 0.5))) + # change color here
  guides(fill = "none") +
  labs(x = "Adult Abundance", y = "Productivity") +
  theme_classic() +
  scale_x_continuous(limits = c(0, max(SPEC_dat$Adult)), expand = c(0, 0)) + # move the y-axis so it intercepts with 0 on x-axis
  theme(axis.title.x = element_text(size=12, family="Arial", margin = margin(t=5)),
        axis.title.y = element_text(size=12, family="Arial", margin = margin(r=5)),
        axis.text = element_text(size=10, family="Arial"))
        

TypeII

# alternatively: if we want the curve without the errorbars

TypeII <- ggplot(data = SPEC_dat, aes(x = Adult, y = .epred)) + 
  geom_smooth(color = "#484554", linewidth=1) + # change color and linewidth here
  guides(fill = "none") +
  labs(x = "Adult Abundance", y = "Productivity") +
  theme_classic() +
  scale_x_continuous(limits = c(0, max(SPEC_dat$Adult)), expand = c(0, 0)) + # move the y-axis so it intercepts with 0 on x-axis
  theme(axis.title.x = element_text(size=12, family="Arial", margin = margin(t=5)),
        axis.title.y = element_text(size=12, family="Arial", margin = margin(r=5)),
        axis.text = element_text(size=10, family="Arial"))

TypeII


# COLOR CODES FOR ANNOTATING PLOT
# Jordan: I used the DarkMint palette to plot Intensity and the Purples palette for Threshold in Fig 3
# these palettes are from the colorspace package
# So you can match with Fig 3, I have plotted the palettes and identified which hex code I used in the histograms
# but any of the hex codes in the palettes should look good
hcl_palettes(palette="DarkMint", n=9, plot=T) # look at palette used above for Intensity
sequential_hcl(9, "DarkMint") # get hex codes
# Intensity: I used "#3F8489" to fill the histogram and "#0E3F5C" to outline the histogram bars

hcl_palettes(palette="Purples", n=9, plot=T) # look at palette used above for Threshold
sequential_hcl(9, "Purples") # get hex codes
# Threshold: I used "#8B7EBB" to fill the histogram and "#3D1778" to outline the histogram bars



#################################################################
### Combine Panels ###

Fig1 <- MAPS_map + TypeII + 
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size=14, family = "Arial", face="bold"))


Fig1

# export plot
ggsave(here("Figures", "Figure1.pdf"), plot = Fig1, 
       width = 30 , height = 15, units = "cm",
       device = cairo_pdf)

