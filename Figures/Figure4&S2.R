###### MAPS Project: Density Dependence #######
### Script name: Figure4.R
### Author(s): SLJ

########### Objective/Description of Script #####################
# create Figure 4, which contains 6 panels from the trait-based analyses
# also make Figure S2, which contains 5 panels showing the remaining results from the trait-based analyses
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

#################################################################
### Import data files ###

dat_LH <- readRDS(here("Outputs", "dat_LH.rds"))
dat_morph <- readRDS(here("Outputs", "dat_morph.rds"))
dat_trophic <- readRDS(here("Outputs", "dat_trophic.rds"))
dat_ss <- readRDS(here("Outputs", "dat_ss.rds"))

# import models
slope_LH_mod_logn <- readRDS(here("Models/TraitModels", "slope_LH_mod_logn.rds"))
minadult_LH_mod_nb <- readRDS(here("Models/TraitModels", "minadult_LH_mod_nb.rds"))
minadult_Morph_mod_nb <- readRDS(here("Models/TraitModels", "minadult_Morph_mod_nb.rds"))
slope_Trophic_mod_logn <- readRDS(here("Models/TraitModels", "slope_Trophic_mod_logn.rds"))
minadult_Trophic_mod_nb <- readRDS(here("Models/TraitModels", "minadult_Trophic_mod_nb.rds"))
slope_SS_mod_logn <- readRDS(here("Models/TraitModels", "slope_SS_mod_logn.rds"))
minadult_SS_mod_nb <- readRDS(here("Models/TraitModels", "minadult_SS_mod_nb.rds"))


#########################################################################
########################### Figure 4 ###################################

###### Life History ######

###  plot of adult survival using threshold as response variable

# find the min and max values for adult survival in the data
min(dat_LH$scAdult.survival) # -2.97
max(dat_LH$scAdult.survival)  # 3.56

# we also need the mean values of the other response variables
mean(dat_LH$scMaximum.longevity) # 0
mean(dat_LH$sclitter_or_clutch_size_n) # 0


# get predicted values for adult survival
survival_minadult_epred <- minadult_LH_mod_nb %>% 
  epred_draws(newdata = tibble(scAdult.survival = seq(-3, 3.6, 0.1), # create a sequence from the min to the max using intervals of 0.1
                               scMaximum.longevity = c(0), # fix max longevity at mean
                               sclitter_or_clutch_size_n = c(0)), re_formula = NA)

# make plot
(minadult_survival_plot <- ggplot(survival_minadult_epred, aes(x = scAdult.survival, y = .epred)) +
    stat_lineribbon(color = "#238b45") + 
    scale_fill_manual(values = colorspace::lighten("#238b45", c(0.95, 0.75, 0.5))) + 
    guides(fill = "none") +
    labs(x = "Scaled Adult Survival", y = "Negative Feedback Threshold") +
    theme_classic() +
    geom_point(data = dat_LH, aes(x= scAdult.survival, y = min_adult), pch = 19, color = "gray30") +
    coord_cartesian(ylim=c(0, 20)) +
    theme(axis.text = element_text(size=11), 
          axis.title.y = element_text(size=12, margin=margin(r=8)),
          axis.title.x = element_text(size=12, margin=margin(t=8))))

### make a plot of longevity using threshold as response variable

# find the min and max values for maximum longevity in the data
min(dat_LH$scMaximum.longevity) # -1.94
max(dat_LH$scMaximum.longevity)  # 3.6

# we also need the mean values of the other variables
mean(dat_LH$sclitter_or_clutch_size_n) # 0
mean(dat_LH$scAdult.survival) # 0

# get predicted values for clutch size
longevity_minadult_epred <- minadult_LH_mod_nb %>% 
  epred_draws(newdata = tibble(scMaximum.longevity = seq(-2, 3.6, 0.1), # create a sequence from the min to the max using intervals of 0.1
                               sclitter_or_clutch_size_n = c(0), # fix at mean
                               scAdult.survival = c(0)), re_formula = NA)

# make plot
(minadult_longevity_plot <- ggplot(longevity_minadult_epred, aes(x = scMaximum.longevity, y = .epred)) +
    stat_lineribbon(color = "#F7710A") + 
    scale_fill_manual(values = colorspace::lighten("#F7710A", c(0.95, 0.75, 0.5))) + 
    guides(fill = "none") +
    labs(x = "Scaled Maximum Longevity", y = "Negative Feedback Threshold") +
    theme_classic() +
    geom_point(data = dat_LH, aes(x= scMaximum.longevity, y = min_adult), pch = 19, color = "gray30") +
    coord_cartesian(ylim=c(0, 22)) +  
    theme(axis.text = element_text(size=11), 
          axis.title.y = element_text(size=12, margin=margin(r=8)),
          axis.title.x = element_text(size=12, margin=margin(t=8))))


### plot of clutch size using intensity as response variable

# find the min and max values for clutch size in the data
min(dat_LH$sclitter_or_clutch_size_n) # -1.35
max(dat_LH$sclitter_or_clutch_size_n)  # 3.4

# we also need the mean values of the other variables
mean(dat_LH$scMaximum.longevity) # 0
mean(dat_LH$scAdult.survival) # 0

# get predicted values for clutch size
clutch_slope_epred <- slope_LH_mod_logn %>% 
  epred_draws(newdata = tibble(sclitter_or_clutch_size_n = seq(-1.4, 3.4, 0.1), # create a sequence from the min to the max using intervals of 0.1
                               scMaximum.longevity = c(0), # fix at mean
                               scAdult.survival = c(0)), re_formula = NA)

# make plot
(slope_clutch_plot <- ggplot(clutch_slope_epred, aes(x = sclitter_or_clutch_size_n, y = .epred)) +
    stat_lineribbon(color = "#56B4E9") + 
    scale_fill_manual(values = colorspace::lighten("#56B4E9", c(0.95, 0.75, 0.5))) + 
    guides(fill = "none") +
    labs(x = "Scaled Clutch Size", y = "Negative Feedback Intensity") +
    theme_classic() +
    geom_point(data = dat_LH, aes(x= sclitter_or_clutch_size_n, y = abs(estimate)), pch = 19, color = "gray30") +
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
    stat_lineribbon(color = "#6518B8") + 
    scale_fill_manual(values = colorspace::lighten("#6518B8", c(0.95, 0.75, 0.5))) + 
    guides(fill = "none") +
    labs(x = "Scaled Trophic Level", y = "Density Dependence Threshold") +
    theme_classic() +
    geom_point(data = dat_trophic, aes(x= scTrophicLevel, y = min_adult), pch = 19, color = "gray30") +
    coord_cartesian(ylim=c(0, 20)) +
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
    stat_lineribbon(color = "#0A9AA8") + 
    scale_fill_manual(values = colorspace::lighten("#0A9AA8", c(0.95, 0.75, 0.5))) + 
    guides(fill = "none") +
    labs(x = "Scaled Sexual Dichromatism", y = "Negative Feedback Threshold") +
    theme_classic() +
    geom_point(data = dat_ss, aes(x= scDC, y = min_adult), pch = 19, color = "gray30") +
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
    stat_lineribbon(color = "#6E940F") + 
    scale_fill_manual(values = colorspace::lighten("#6E940F", c(0.95, 0.75, 0.5))) + 
    guides(fill = "none") +
    labs(x = "Scaled Male Sexual Selection", y = "Negative Feedback Intensity") +
    theme_classic() +
    geom_point(data = dat_ss, aes(x= sc_ssM, y = abs(estimate)), pch = 19, color = "gray30") +
    coord_cartesian(ylim=c(0, 0.1)) +
    theme(axis.text = element_text(size=11), 
          axis.title.y = element_text(size=12, margin=margin(r=8)),
          axis.title.x = element_text(size=12, margin=margin(t=8))))


##### Combine Plots into 6 Panel Figure #####

(minadult_survival_plot + minadult_longevity_plot + slope_clutch_plot) /
(minadult_trophic_plot + minadult_DC_plot + slope_ssM_plot) + 
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size=14, family = "Arial", face="bold"))


#########################################################################
########################### Figure S2 ###################################

##### Life History #####
### plot of clutch size with threshold as response variable

# find the min and max values for clutch size in the data
min(dat_LH$sclitter_or_clutch_size_n) # -1.35
max(dat_LH$sclitter_or_clutch_size_n)  # 3.39

# we also need the mean values of the other variables
mean(dat_LH$scMaximum.longevity) # 0
mean(dat_LH$scAdult.survival) # 0

# get predicted values for clutch size
clutch_minadult_epred <- minadult_LH_mod_nb %>% 
  epred_draws(newdata = tibble(sclitter_or_clutch_size_n = seq(-1.4, 3.4, 0.1), # create a sequence from the min to the max using intervals of 0.1
                               scMaximum.longevity = c(0), # fix at mean
                               scAdult.survival = c(0)), re_formula = NA)

# make plot
(minadult_clutch_plot <-ggplot(clutch_minadult_epred, aes(x = sclitter_or_clutch_size_n, y = .epred)) +
    stat_lineribbon(color = "#56B4E9") + 
    scale_fill_manual(values = colorspace::lighten("#56B4E9", c(0.95, 0.75, 0.5))) + 
    guides(fill = "none") +
    labs(x = "Scaled Clutch Size", y = "Negative Feedback Threshold") +
    theme_classic() +
    geom_point(data = dat_LH, aes(x= sclitter_or_clutch_size_n, y = min_adult), pch = 19, color = "gray30") +
    coord_cartesian(ylim=c(0, 20)) +
    theme(axis.text = element_text(size=11), 
          axis.title.y = element_text(size=12, margin=margin(r=8)),
          axis.title.x = element_text(size=12, margin=margin(t=8))))



### plot of adult survival with intensity as response variable

# find the min and max values for adult survival in the data
min(dat_LH$scAdult.survival) # -2.97
max(dat_LH$scAdult.survival)  # 3.55

# we also need the mean values of the other variables
mean(dat_LH$scMaximum.longevity) # 0
mean(dat_LH$sclitter_or_clutch_size_n) # 0

# get predicted values for adult survival
survival_slope_epred <- slope_LH_mod_logn %>% 
  epred_draws(newdata = tibble(scAdult.survival = seq(-3, 3.6, 0.1), # create a sequence from the min to the max using intervals of 0.1
                               scMaximum.longevity = c(0), # fix at mean
                               sclitter_or_clutch_size_n = c(0)), re_formula = NA)

# make plot
(slope_survival_plot <- ggplot(survival_slope_epred, aes(x = scAdult.survival, y = .epred)) +
    stat_lineribbon(color = "#238b45") + 
    scale_fill_manual(values = colorspace::lighten("#238b45", c(0.95, 0.75, 0.5))) + 
    guides(fill = "none") +
    labs(x = "Scaled Adult Survival", y = "Negative Feedback Intensity") +
    theme_classic() +
    geom_point(data = dat_LH, aes(x= scAdult.survival, y = abs(estimate)), pch = 19, color = "gray30") +
    coord_cartesian(ylim=c(0, 0.15)) +
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
    stat_lineribbon(color = "#bf0404") + 
    scale_fill_manual(values = colorspace::lighten("#bf0404", c(0.95, 0.75, 0.5))) + 
    guides(fill = "none") +
    labs(x = "Scaled Beak Shape", y = "Negative Feedback Threshold") +
    theme_classic() +
    geom_point(data = dat_morph, aes(x= scBeak_PC2, y = min_adult), pch = 19, color = "gray30") +
    coord_cartesian(ylim=c(0, 20)) +
    theme(axis.text = element_text(size=11), 
          axis.title.y = element_text(size=12, margin=margin(r=8)),
          axis.title.x = element_text(size=12, margin=margin(t=8))))

##### Trophic #####
###  plot of strata using intensity as response variable

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
    stat_lineribbon(color = "#AF088F") + 
    scale_fill_manual(values = colorspace::lighten("#AF088F", c(0.95, 0.75, 0.5))) + 
    guides(fill = "none") +
    labs(x = "Scaled Trophic Strata", y = "Negative Feedback Intensity") +
    theme_classic() +
    geom_point(data = dat_trophic, aes(x= scStrata, y = abs(estimate)), pch = 19, color = "gray30") +
    coord_cartesian(ylim=c(0, 0.18)) +
    theme(axis.text = element_text(size=11), 
          axis.title.y = element_text(size=12, margin=margin(r=8)),
          axis.title.x = element_text(size=12, margin=margin(t=8))))


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
(slope_DC_plot <- ggplot(DC_slope_epred, aes(x = scDC, y = .epred)) +
    stat_lineribbon(color = "#0A9AA8") + 
    scale_fill_manual(values = colorspace::lighten("#0A9AA8", c(0.95, 0.75, 0.5))) + 
    guides(fill = "none") +
    labs(x = "Scaled Sexual Dichromatism", y = "Negative Feedback Intensity") +
    theme_classic() +
    geom_point(data = dat_ss, aes(x= scDC, y = abs(estimate)), pch = 19, color = "gray30") +
    coord_cartesian(ylim=c(0, 0.1)) +
    theme(axis.text = element_text(size=11), 
          axis.title.y = element_text(size=12, margin=margin(r=8)),
          axis.title.x = element_text(size=12, margin=margin(t=8))))


##### Combine Plots into 5 Panel Figure #####

(minadult_clutch_plot + slope_survival_plot + minadult_beakPC2_plot +
 slope_strata_plot + slope_DC_plot) + 
  plot_layout(ncol=3) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size=14, family = "Arial", face="bold"))
