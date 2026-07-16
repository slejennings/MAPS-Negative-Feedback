###### MAPS Project: Density Dependence #######
### Script name: Figure2.R
### Author(s): SLJ, JML

########### Objective/Description of Script #####################
# create Figure 2 
# Panel 1: distribution for threshold and intensity and the correlation between them
#################################################################

### Setup ###

library(here)
library(ggExtra)
library(ggdist)
library(ggnewscale)
library(colorspace)
library(Polychrome)
library(tidyverse)
library(patchwork)

#########################################################################
### Load required files ###

# import data to make species abundance-productivity curve
epred_df_m2 <- readRDS(here("Outputs", "epred_df_m2.rds"))


#########################################################################
# We are using 3 representative species to use as examples of a Type I, Type II and Type III productivity-abundance curve
# Extract a data frame for each of our 3 chosen species

# Type I = BCCH
BCCH_dat <- epred_df_m2 %>%
  filter(SPEC == "BCCH") # specify the species to keep
  
# Type II - KEWA
KEWA_dat <- epred_df_m2 %>%
  filter(SPEC == "KEWA") # specify the species to keep
  
# Type III - BEWR
BEWR_dat <- epred_df_m2 %>%
  filter(SPEC == "BEWR") # specify the species to keep
  
#########################################################################
# Make a plot for each species

# note: if you change the species, you'll need edit ggplot(data =) and scale_x_continuous(limits =)
typeIplot <- ggplot(data = BCCH_dat, aes(x = Adult, y = .epred)) + 
  ggdist::stat_lineribbon(color = "#DA4167") + # change color here
  scale_fill_manual(values = colorspace::lighten("#DA4167", c(0.95, 0.75, 0.5))) + # change color here
  guides(fill = "none") +
  labs(x = "Adult Abundance", y = "Productivity", title = "Type 1") +
  theme_classic() +
  scale_x_continuous(limits = c(0, max(BCCH_dat$Adult)), expand = c(0, 0)) + # move the y-axis so it intercepts with 0 on x-axis
  theme(axis.title.x = element_text(size=12, family="Arial", margin = margin(t=5)),
        axis.title.y = element_text(size=12, family="Arial", margin = margin(r=5)),
        axis.text = element_text(size=10, family="Arial"), 
        plot.title = element_text(size = 14, hjust = 0.5, family = "Arial"))
typeIplot

# we probably want to add a title or text annotation to this plot that says "Type I"
# also, it could be useful to differentiate the types using colors



typeIIplot <- ggplot(data = KEWA_dat, aes(x = Adult, y = .epred)) + 
  ggdist::stat_lineribbon(color = "#29335C") + # change color here
  scale_fill_manual(values = colorspace::lighten("#29335C", c(0.95, 0.75, 0.5))) + # change color here
  guides(fill = "none") +
  labs(x = "Adult Abundance", y = "Productivity", title = "Type 2") +
  theme_classic() +
  scale_x_continuous(limits = c(0, max(KEWA_dat$Adult)), expand = c(0, 0)) + # move the y-axis so it intercepts with 0 on x-axis
  theme(axis.title.x = element_text(size=12, family="Arial", margin = margin(t=5)),
        axis.title.y = element_text(size=12, family="Arial", margin = margin(r=5)),
        axis.text = element_text(size=10, family="Arial"), plot.title = element_text(size = 14, hjust = 0.5, family = "Arial"))
typeIIplot 

typeIIIplot <- ggplot(data = BEWR_dat, aes(x = Adult, y = .epred)) + 
  ggdist::stat_lineribbon(color = "#F5AF00") + 
  scale_fill_manual(values = colorspace::lighten("#F5AF00", c(0.95, 0.75, 0.5))) + # change color here
  guides(fill = "none") +
  labs(x = "Adult Abundance", y = "Productivity", title = "Type 3") +
  theme_classic() +
  scale_x_continuous(limits = c(0, max(BEWR_dat$Adult)), expand = c(0, 0)) + # move the y-axis so it intercepts with 0 on x-axis
  theme(axis.title.x = element_text(size=12, family="Arial", margin = margin(t=5)),
        axis.title.y = element_text(size=12, family="Arial", margin = margin(r=5)),
        axis.text = element_text(size=10, family="Arial"), plot.title = element_text(size = 14, hjust = 0.5, family = "Arial"))
typeIIIplot 

#########################################################################

# Combine plots

# with plot annotations (A, B, C)
typeIplot + typeIIplot + typeIIIplot + 
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size=14, family = "Arial", face="bold"))

# or without plot annotations
typeIplot + typeIIplot + typeIIIplot + 
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size=14, family = "Arial", face="bold"))

