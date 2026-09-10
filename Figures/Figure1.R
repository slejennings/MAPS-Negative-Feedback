###### MAPS Project: Negative Feedback #######
### Script name: Figure1.R
### Author(s): SLJ, JML

########### Objective/Description of Script #####################
# create Figure 1, which contains two panels
# Panel 1: map showing locations of MAPS banding stations
# Panel 2: predicted relationship between productivity and abundance and the two extracted negative feedback variables
#################################################################

### Setup ###

# load packages
library(here)
library(tidyverse)
library(tidybayes)
library(bayesplot)
library(ggspatial)
library(sf)
library(colorspace)
library(ggdist)
library(ggnewscale)
library(patchwork)

#################################################################
### Import data files ###

# import productivity-abundance data
NF_breed_dat <- readRDS(here("Outputs", "NF_breed_dat.rds"))

# import station information 
station_info <- read.csv(here("Data", "MAPS_STATION_location_and_operations.csv"), header=T) %>%
  mutate(STA = factor(STA)) # convert STA to a factor

# import data to make species abundance-productivity curve
epred_df_m2 <- readRDS(here("Outputs", "epred_df_m2.rds"))

#################################################################
### Create Panel A of Figure 1: Maps of MAPS banding stations ###

# use productivity-abundance to generate a list of MAPS stations that were in the final analysis
STAlist <- NF_breed_dat %>%
  select(STA) %>% distinct() 
nrow(STAlist) # 378 MAPS stations

# list of stations with number of species (out of 62 in the analysis) at each station
STAlist <- NF_breed_dat %>%
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
STA_sf <- st_as_sf(STAlist_coords, coords = c(2:3), crs="EPSG:4269")

# get map of the world and subset USA and Canada
world_map <- map_data("world")
na_map <- subset(world_map, region %in% c("USA", "Canada"))

# plot map layer
mapbase <- ggplot(data = na_map, aes(x = long, y = lat, group = group)) +
  geom_polygon(fill = "white", color = "black") +
  theme_minimal()

# add stations to the map
MAPS_map <- mapbase +
  layer_spatial(data = STA_sf, color = "#29335C", alpha = 0.5, size = 2.5) +
  coord_sf(crs="EPSG:4269", xlim=c(-149,-57), ylim=c(27, 69)) +
  xlab("Longitude") + ylab("Latitude") +
  theme(axis.title.x = element_text(size = 12, margin = margin(t=5), family = "Arial"),
        axis.title.y = element_text(size = 12, margin = margin(r=5), family = "Arial"),
        axis.text = element_text(size=11, family = "Arial"))

MAPS_map

#################################################################
### Create Panel B of Figure 1: Species curve with extracted negative feedback variables ###

# get data for species with a Type 2 curve
# using NOCA
SPEC_dat <- epred_df_m2 %>%
  filter(SPEC == "NOCA") # specify the species to keep

# split the curve into the sections that represent threshold and intensity, respectively
threshold_curve <- SPEC_dat %>% filter(Adult<=5) %>% mutate(group="A")
intensity_curve <- SPEC_dat %>% filter(Adult>4) %>% mutate(group="B")

# get points where we obtained the slope/intensity
slopepoints <- intensity_curve %>% 
  group_by(Adult) %>% 
  summarize(predictprod = mean(.epred)) 

# points with some jitter to place them above the curve
slopepoints_jitter <- intensity_curve %>% 
  group_by(Adult) %>% 
  summarize(predictprod =mean(.epred)) %>% 
  mutate(predictprod = predictprod + 0.02) # add a small amount of jitter to move points up above the curve

# create plot of single species with a Type 2 curve
# using non-jittered points here. To change that, alter the data for geom_point() to use slopepoints_jitter

Type2 <- ggplot() +
  ggdist::stat_lineribbon(data = threshold_curve, aes(x = Adult, y = .epred, color=group)) + 
  scale_fill_manual(values = c("#EDECF9", "#D8D4EC", "#BFBADD")) + # change color of threshold error bars here
  guides(fill = "none", color="none") +
  ggnewscale::new_scale_fill() +
  ggdist::stat_lineribbon(data = intensity_curve, aes(x = Adult, y = .epred, color = group)) +
  scale_fill_manual(values = c("#D1FBD4", "#AFE4C4", "#72B4A6")) +
  scale_color_manual(values=c("#715DAA","#2A6D7A")) + # change color of curve lines here
  guides(fill = "none", color = "none") +
  labs(x = "Adult Abundance", y = "Per Capita Productivity") +
  theme_classic() +
  scale_x_continuous(limits = c(0, max(SPEC_dat$Adult)), expand = c(0, 0)) + # move the y-axis so it intercepts with 0 on x-axis
  theme(axis.title.x = element_text(size=12, family="Arial", margin = margin(t=5)),
        axis.title.y = element_text(size=12, family="Arial", margin = margin(r=5)),
        axis.text = element_text(size=11, family="Arial")) +
  geom_point(data = slopepoints, aes(x=Adult, y=predictprod), size=2.25, color="#0E3F5C") + # add points and set their size, shape, color
  xlim(0, 35) # set limits for the x-axis


# note: will likely give a warning message about missing values because we are trimming the x-axis so some data is not being plotted

#################################################################
### Add annotations for Panel B of Figure 1 ###

# add lines for threshold + intensity, labels
Type2_annotate <- Type2 + annotate("segment", x= 5, xend = 5, y = 0.7, yend = 0, colour = "#583A99", linewidth = 1.2, linetype = 2) + 
  annotate("segment", x = 5, xend = 35, y = 0.7, yend = 0.7, linewidth = 1.2, color = "#0E3F5C", linetype = 7, arrow = arrow(angle = 30, length = unit(0.2, "inches"), ends = "both", type = "open")) + 
  annotate("label", label = "Sampling range for negative \n feedback intensity", x = 20, y = 0.6, colour = "#0E3F5C", family = 'Arial', size = 4.5, label.padding = unit(0.5, "lines")) + 
  annotate("label", label = "Negative feedback \n threshold", x = 7, y = 0.13, colour = "#583A99", family = 'Arial', size = 4.5, label.padding = unit(0.5, "lines"))
Type2_annotate

#################################################################
### Combine Panels ###

(Fig1 <- MAPS_map / Type2_annotate + 
  plot_layout(width = unit(c(11, 13), c('cm', 'cm')), heights= unit(c(8, 9), c('cm', 'cm'))) &
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size=14, family = "Arial", face="bold")))

Fig1

# export plot
ggsave(here("Figures", "Figure1.pdf"), plot = Fig1, 
       width = 15 , height = 22, units = "cm",
       device = cairo_pdf)

ggsave(here("Figures", "Figure1.png"), plot = Fig1, 
       width = 15 , height = 22, units = "cm")
