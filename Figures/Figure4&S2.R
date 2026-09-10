###### MAPS Project: Negative Feedback #######
### Script name: Figure4.R
### Author(s): SLJ

########### Objective/Description of Script #####################
# create Figure 4, which contains 6 panels from the trait-based analyses
# also make Figure S2, which the remaining results from the trait-based analyses
#################################################################

### Setup ###

# load packages
library(here)
library(dplyr)
library(ggplot2)
library(tidybayes)
library(bayesplot)
library(colorspace)
library(ggdist)
library(ggnewscale)
library(patchwork)
library(scales)

#################################################################
### Import data files ###

dat_LH <- readRDS(here("Outputs", "dat_LH.rds"))
dat_morph <- readRDS(here("Outputs", "dat_morph.rds"))
dat_trophic <- readRDS(here("Outputs", "dat_trophic.rds"))
dat_ss <- readRDS(here("Outputs", "dat_ss.rds"))

# import models
slope_LH_mod_logn <- readRDS(here("Models/TraitModels", "slope_LH_mod_logn.rds"))
minadult_Morph_mod_nb <- readRDS(here("Models/TraitModels", "minadult_Morph_mod_nb.rds"))
slope_Trophic_mod_logn <- readRDS(here("Models/TraitModels", "slope_Trophic_mod_logn.rds"))
minadult_Trophic_mod_nb <- readRDS(here("Models/TraitModels", "minadult_Trophic_mod_nb.rds"))
slope_SS_mod_logn <- readRDS(here("Models/TraitModels", "slope_SS_mod_logn.rds"))
minadult_SS_mod_nb <- readRDS(here("Models/TraitModels", "minadult_SS_mod_nb.rds"))


#########################################################################
########################### Figure 4 ###################################

###### Life History ######

### plot of clutch size using intensity as response variable

# find the min and max values for clutch size in the data
min(dat_LH$sclitter_or_clutch_size_n) # -1.35
max(dat_LH$sclitter_or_clutch_size_n)  # 3.4

# we also need the mean values of the other variables
mean(dat_LH$scMaximum.longevity) # 0

# get predicted values for clutch size
clutch_slope_epred <- slope_LH_mod_logn %>% 
  epred_draws(newdata = tibble(sclitter_or_clutch_size_n = seq(-1.4, 3.4, 0.1), # create a sequence from the min to the max using intervals of 0.1
                               scMaximum.longevity = c(0)), # fix at mean
                               re_formula = NA)

# make plot
(slope_clutch_plot <- ggplot(clutch_slope_epred, aes(x = sclitter_or_clutch_size_n, y = .epred)) +
    stat_lineribbon(color = "#4D5382", linetype = "dashed") + # color is for the line
    scale_fill_manual(values = colorspace::lighten("#7B82B2", c(0.95, 0.75, 0.5))) + # fill is for the ribbon 
    guides(fill = "none") +
    labs(x = "Scaled Clutch Size", y = "Negative Feedback Intensity") +
    theme_classic() +
    geom_point(data = dat_LH, aes(x= sclitter_or_clutch_size_n, y = abs(estimate)), pch = 19, color = "#4D5382") +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.01)) +
    coord_cartesian(ylim=c(0, 0.15)) +
    theme(axis.text = element_text(size=11), 
          axis.title.y = element_text(size=12, margin=margin(r=8)),
          axis.title.x = element_text(size=12, margin=margin(t=8))))

#### Trophic #####
### make a plot of trophic level using threshold as response variable

# find the min and max values for trophic level in the data
min(dat_trophic$scTrophicLevel) # -2.3
max(dat_trophic$scTrophicLevel)  # 1.37

# we also need the mean values of the other variables
mean(dat_trophic$scStrata) # 0
mean(dat_trophic$scDiet) # 0

# get predicted values for trophic level
trophiclevel_minadult_epred <- minadult_Trophic_mod_nb %>% 
  epred_draws(newdata = tibble(scTrophicLevel = seq(-2.3, 1.4, 0.1), # create a sequence from the min to the max using intervals of 0.1
                               scStrata = c(0), # fix at mean
                               scDiet = c(0)), re_formula = NA)

# make plot
(minadult_trophic_plot <- ggplot(trophiclevel_minadult_epred, aes(x = scTrophicLevel, y = .epred)) +
    stat_lineribbon(color = "#488960") + 
    scale_fill_manual(values = colorspace::lighten("#57a774", c(0.95, 0.75, 0.5))) + 
    guides(fill = "none") +
    labs(x = "Scaled Trophic Level", y = "Negative Feedback Threshold") +
    theme_classic() +
    geom_point(data = dat_trophic, aes(x= scTrophicLevel, y = min_adult), pch = 19, color = "#488960") +
    coord_cartesian(ylim=c(0, 20)) +
    theme(axis.text = element_text(size=11), 
          axis.title.y = element_text(size=12, margin=margin(r=8)),
          axis.title.x = element_text(size=12, margin=margin(t=8))))

###  plot of foraging strata using intensity as response variable

# find the min and max values for strata in the data
min(dat_trophic$scStrata) # -1.82
max(dat_trophic$scStrata)  # 2.41

# we also need the mean values of the other variables
mean(dat_trophic$scTrophicLevel) # 0
mean(dat_trophic$scDiet) # 0

# get predicted values for Strata
strata_slope_epred <- slope_Trophic_mod_logn %>% 
  epred_draws(newdata = tibble(scStrata = seq(-1.9, 2.5, 0.1), # create a sequence from the min to the max using intervals of 0.1
                               scTrophicLevel = c(0), # fix at mean
                               scDiet = c(0)), re_formula = NA)

# make plot
(slope_strata_plot <- ggplot(strata_slope_epred, aes(x = scStrata, y = .epred)) +
    stat_lineribbon(color = "#488960", linetype="dashed") + 
    scale_fill_manual(values = colorspace::lighten("#57a774", c(0.95, 0.75, 0.5))) + # 
    guides(fill = "none") +
    labs(x = "Scaled Foraging Strata Generalism", y = "Negative Feedback Intensity") +
    theme_classic() +
    geom_point(data = dat_trophic, aes(x= scStrata, y = abs(estimate)), pch = 19, color = "#488960") +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.01)) +
    coord_cartesian(ylim=c(0, 0.18)) +
    theme(axis.text = element_text(size=11), 
          axis.title.y = element_text(size=12, margin=margin(r=8)),
          axis.title.x = element_text(size=12, margin=margin(t=8))))


##### Sexual Selection #####

### plot of dichromatism using threshold as response variable

# find the min and max values for sexual dichromatism in the data
min(dat_ss$scDC) # -1.4
max(dat_ss$scDC)  # 2.62

# we also need the mean values of the other variables
mean(dat_ss$scDM) # 0
mean(dat_ss$sc_ssM) # 0
mean(dat_ss$sc_ssF) # 0

# get predicted values for strata
DC_minadult_epred <- minadult_SS_mod_nb %>% 
  epred_draws(newdata = tibble(scDC = seq(-1.4, 2.7, 0.1), # create a sequence from the min to the max using intervals of 0.1
                               scDM = c(0), # fix at mean
                               sc_ssM = c(0), 
                               sc_ssF = c(0)), re_formula = NA)

# make plot
(minadult_DC_plot <- ggplot(DC_minadult_epred, aes(x = scDC, y = .epred)) +
    stat_lineribbon(color = "#06666F") + 
    scale_fill_manual(values = colorspace::lighten("#0A9AA8", c(0.95, 0.75, 0.5))) + 
    guides(fill = "none") +
    labs(x = "Scaled Plumage Sexual Dichromatism", y = "Negative Feedback Threshold") +
    theme_classic() +
    geom_point(data = dat_ss, aes(x= scDC, y = min_adult), pch = 19, color = "#06666F") +
    coord_cartesian(ylim=c(0, 20)) + 
    theme(axis.text = element_text(size=11), 
          axis.title.y = element_text(size=12, margin=margin(r=8)),
          axis.title.x = element_text(size=12, margin=margin(t=8))))


###  plot of male sexual selection using intensity (slope) as response variable

# find the min and max values for male sexual selection in the data
min(dat_ss$sc_ssM) # -0.84
max(dat_ss$sc_ssM)  # 2.77

# we also need the mean values of the other variables
mean(dat_ss$scDC) # 0
mean(dat_ss$scDM) # 0
mean(dat_ss$sc_ssF) # 0

# get predicted values for strata
ssM_slope_epred <- slope_SS_mod_logn %>% 
  epred_draws(newdata = tibble(sc_ssM = seq(-0.9, 2.8, 0.1), # create a sequence from the min to the max using intervals of 0.1
                               scDC = c(0), # fix at mean
                               scDM = c(0), 
                               sc_ssF = c(0)), re_formula = NA)

# make plot
(slope_ssM_plot <- ggplot(ssM_slope_epred, aes(x = sc_ssM, y = .epred)) +
    stat_lineribbon(color = "#06666F", linetype="dashed") + 
    scale_fill_manual(values = colorspace::lighten("#0A9AA8", c(0.95, 0.75, 0.5))) + 
    guides(fill = "none") +
    labs(x = "Scaled Degree of Polygyny in Males", y = "Negative Feedback Intensity") +
    theme_classic() +
    geom_point(data = dat_ss, aes(x= sc_ssM, y = abs(estimate)), pch = 19, color = "#06666F") +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.01)) +
    coord_cartesian(ylim=c(0, 0.1)) +
    theme(axis.text = element_text(size=11), 
          axis.title.y = element_text(size=12, margin=margin(r=8)),
          axis.title.x = element_text(size=12, margin=margin(t=8))))


##### Morphology #####
### plot of beak PC2 using threshold (min # of adults) as response variable

# find the min and max values for Beak PC2 in the data
min(dat_morph$scBeak_PC2) # -2.45
max(dat_morph$scBeak_PC2)  # 2.32

# we also need the mean values of the other response variables
mean(dat_morph$scHWI) # 0
mean(dat_morph$scLNMass) # 0
mean(dat_morph$scBeak_PC1) # 0

# get predicted values for strata
DM_minadult_epred <- minadult_Morph_mod_nb %>% 
  epred_draws(newdata = tibble(scBeak_PC2 = seq(-2.5, 2.4, 0.1), # create a sequence from the min to the max using intervals of 0.1
                               scBeak_PC1 = c(0), # fix at mean
                               scHWI = c(0), 
                               scLNMass = c(0)), re_formula = NA)

# make plot
(minadult_beakPC2_plot <- ggplot(DM_minadult_epred, aes(x = scBeak_PC2, y = .epred)) +
    stat_lineribbon(color = "#4F6B0B") + 
    scale_fill_manual(values = colorspace::lighten("#6E940F", c(0.95, 0.75, 0.5))) + 
    guides(fill = "none") +
    labs(x = "Scaled Beak Shape", y = "Negative Feedback Threshold") +
    theme_classic() +
    geom_point(data = dat_morph, aes(x= scBeak_PC2, y = min_adult), pch = 19, color = "#4F6B0B") +
    coord_cartesian(ylim=c(0, 20)) +
    theme(axis.text = element_text(size=11), 
          axis.title.y = element_text(size=12, margin=margin(r=8)),
          axis.title.x = element_text(size=12, margin=margin(t=8))))



##### Combine Plots into 6 Panel Figure #####

(Fig4 <- slope_clutch_plot + minadult_trophic_plot + slope_strata_plot +
  minadult_beakPC2_plot + minadult_DC_plot + slope_ssM_plot + 
  plot_layout(ncol=3) &
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size=14, family = "Arial", face="bold")))

# export

ggsave(here("Figures", "Figure4.pdf"), plot = Fig4, 
       width = 28 , height = 20, units = "cm",
       device = cairo_pdf)

ggsave(here("Figures", "Figure4.png"), plot = Fig4, 
       width = 28 , height = 20, units = "cm")


#########################################################################
########################### Figure S2 ###################################


##### Sexual Selection #####

### plot of dichromatism using intensity (slope) as response variable

# find the min and max values for sexual dichromatism in the data
min(dat_ss$scDC) # -1.4
max(dat_ss$scDC)  # 2.62

# we also need the mean values of the other variables
mean(dat_ss$scDM) # 0
mean(dat_ss$sc_ssM) # 0
mean(dat_ss$sc_ssF) # 0

# get predicted values for strata
DC_slope_epred <- slope_SS_mod_logn %>% 
  epred_draws(newdata = tibble(scDC = seq(-1.4, 2.7, 0.1), # create a sequence from the min to the max using intervals of 0.1
                               scDM = c(0), # fix at mean
                               sc_ssM = c(0), 
                               sc_ssF = c(0)), re_formula = NA)

# make plot
(FigS2 <- ggplot(DC_slope_epred, aes(x = scDC, y = .epred)) +
    stat_lineribbon(color = "#06666F", linetype="dashed") + 
    scale_fill_manual(values = colorspace::lighten("#0A9AA8", c(0.95, 0.75, 0.5))) + 
    guides(fill = "none") +
    labs(x = "Scaled Plumage Sexual Dichromatism", y = "Negative Feedback Intensity") +
    theme_classic() +
    geom_point(data = dat_ss, aes(x= scDC, y = abs(estimate)), pch = 19, color = "#06666F") +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.01)) +
    coord_cartesian(ylim=c(0, 0.1)) +
    theme(axis.text = element_text(size=11), 
          axis.title.y = element_text(size=12, margin=margin(r=8)),
          axis.title.x = element_text(size=12, margin=margin(t=8))))

ggsave(here("Figures", "FigureS2.pdf"), plot = FigS2, 
       width = 28 , height = 20, units = "cm",
       device = cairo_pdf)

ggsave(here("Figures", "FigureS2.png"), plot = FigS2, 
       width = 12 , height = 10, units = "cm")
