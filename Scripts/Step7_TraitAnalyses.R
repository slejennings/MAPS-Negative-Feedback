###### MAPS Project: Density Dependence #######
### Script name: Step7_TraitAnalyses.R
### Author(s): SLJ, CDF

########### Objective/Description of Script #####################
# model relationships between density dependence metrics and species traits using models that take phylogeny into account
#################################################################


### Setup ###

# load packages
library(tidyverse)
library(ape)
library(geiger)
library(here)
library(marginaleffects)
library(performance)
library(brms)
library(bayesplot)
library(loo)
library(tidybayes)
library(broom.mixed)
library(patchwork)

# set theme for plots
theme_set(theme_bw())


# model options
options(
  mc.cores = 4,  # distribute the 4 chains across 4 cores
  brms.threads = 4,
  brms.backend = "cmdstanr",
  brms.file_refit = "on_change"
)

CHAINS <- 4
ITER <- 5000
WARMUP <- 2000
THIN <- 2
BAYES_SEED <- 963


####################################################
### Import files ###

# Import phylogenetic tree
tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))


# Import trait files from previous step

# morphometric traits
DD_morph <- readRDS(here("Outputs", "mapsDD_morphometrics_m2.rds"))
head(DD_morph)
class(DD_morph)

# life history traits
DD_life <- readRDS(here("Outputs", "mapsDD_lifehistory_m2.rds"))
head(DD_life)

# trophic traits
DD_trophic <- readRDS(here("Outputs", "mapsDD_trophic_m2.rds"))
head(DD_trophic)

# sexual selection traits
DD_SS <- readRDS(here("Outputs", "mapsDD_sexualselection_m2.rds"))
head(DD_SS)

##########################################################################################
#################### Density Dependence Metrics and Life History Traits ##################
##########################################################################################

### Join and prune tree and create covariance matrix ###


# replace blanks in Genus species with underscore
dat_LH <- DD_life %>%
  mutate( # make two columns with different names that are otherwise identical with underscore between genus and species
    rownames = str_replace(Species3_BirdTree, " ", "_"), # this one will be moved to rownames
    BirdTree = str_replace(Species3_BirdTree, " ", "_")) %>% # this one will remain in the data frame
  column_to_rownames(., var = "rownames")

head(dat_LH)

dat_LH$Species <- dat_LH$BirdTree

# link data to tree
LH_tree <- geiger::treedata(tree_out, dat_LH, sort=T)

plot(LH_tree$phy, cex=0.2)

# construct a covariance matrix on the relationship between species 
LH_cov <- ape::vcv.phylo(LH_tree$phy, corr = T)
# corr=T should make the diagonal = 1, which should help with posterior for random effect of tree.

##################################################################################
### Examine the variables ###

# two response variables are estimate (Intensity) and min_adult (Threshold)
# examine distributions for both
hist(abs(DD_life$min_adult))
hist(log(abs(DD_life$min_adult))) # log transformed

hist(DD_life$estimate) # raw
hist(abs(DD_life$estimate)) # take absolute value
hist(log(abs(DD_life$estimate))) # take absolute value and log transform

# standardize predictors that will be used as fixed effects in models
dat_LH <- dat_LH |>
  dplyr::mutate(
    scMaximum.longevity      = as.numeric(scale(Maximum.longevity)),
    sclitter_or_clutch_size_n= as.numeric(scale(litter_or_clutch_size_n)),
    scAdult.survival         = as.numeric(scale(Adult.survival)),
    scbroodvalue             = as.numeric(scale(broodvalue))
  )

##################################################################################
### Life History Model 1 ###

# Response variable: number of adults at onset of DD effects (min_adult) 
# Explanatory variables: longevity, clutch size, adult survival, brood value
# Using negative binomial distribution for model

# center Intercept near the typical count on the log scale
mu0_nb_LH <- mean(dat_LH$min_adult, na.rm = TRUE)
# small guard against log(0)
log_mu0_nb_LH <- log(ifelse(is.finite(mu0_nb_LH) && mu0_nb_LH > 0, mu0_nb_LH, 1))

# specify priors
priors_nb_LH <- c(
  # fixed effects (assumes predictors are standardized)
  set_prior("normal(0, 1)", class = "b"),
  
  # intercept centered on log(mean y) with robust tails
  set_prior(paste0("student_t(3, ", round(log_mu0_nb_LH, 3), ", 2.5)"), class = "Intercept"),
  
  # phylogenetic random intercept SD (correlated by the tree)
  set_prior("student_t(3, 0, 1)",  class = "sd", group = "BirdTree"),
  
  # NB overdispersion (shape >= 0): weak preference for some overdispersion
  set_prior("exponential(0.5)", class = "shape")
)

priors_nb_LH


# run NB model
minadult_LH_mod_nb <- brm(
  min_adult ~ scMaximum.longevity + sclitter_or_clutch_size_n + scAdult.survival + scbroodvalue + (1|gr(BirdTree, cov = LH_cov)),
  data = dat_LH,
  family = negbinomial(),
  data2 = list(LH_cov = LH_cov),
  control = list(adapt_delta = 0.999), # increased to 0.999
  backend = "cmdstanr",
  init = "random",
  save_pars = save_pars(all = TRUE),
  prior = priors_nb_LH,
  chains = CHAINS, iter = ITER, warmup = WARMUP, thin = THIN, seed = BAYES_SEED
)
# no divergent transitions


### Examine the model ###
summary(minadult_LH_mod_nb)
tidy(minadult_LH_mod_nb)
# Rhat and ESS are fine

fixef(minadult_LH_mod_nb, probs=c(0.05, 0.95)) # look at 90% CrI
fixef(minadult_LH_mod_nb, probs=c(0.075, 0.925)) # look at 85% CrI

# negative relationship for adult survival. 95% CrI does not overlap zero
# negative trend for clutch size. 85% CrI does not overlap zero
# positive relationship for longevity. 90% CrI does not overlap zero

# trace and density plots
color_scheme_set("brightblue")
plot(minadult_LH_mod_nb, nvariables = 4, ask = FALSE, theme=theme_minimal())

# posterior predictive checks
pp_check(minadult_LH_mod_nb)
pp_check(minadult_LH_mod_nb, type="hist", ndraws = 11) 
pp_check(minadult_LH_mod_nb, type="ecdf_overlay")
pp_check(minadult_LH_mod_nb, type="intervals", ndraws = 100)

### Plot results ###

plot(conditional_effects(minadult_LH_mod_nb, effects="scAdult.survival"), points = TRUE)
plot(conditional_effects(minadult_LH_mod_nb, effects="sclitter_or_clutch_size_n"), points = TRUE)
plot(conditional_effects(minadult_LH_mod_nb, effects="scMaximum.longevity"), points = TRUE) 

### make a fancy plot of adult survival for life history model using minimum adults (threshold) of DD as response variable

# find the min and max values for adult survival in the data
min(dat_LH$scAdult.survival) # -2.97
max(dat_LH$scAdult.survival)  # 3.56

# we also need the mean values of the other response variables
mean(dat_LH$scMaximum.longevity) # 0
mean(dat_LH$sclitter_or_clutch_size_n) # 0
mean(dat_LH$scbroodvalue) # 0

# get predicted values for adult survival
survival_minadult_epred <- minadult_LH_mod_nb %>% 
  epred_draws(newdata = tibble(scAdult.survival = seq(-3, 3.6, 0.1), # create a sequence from the min to the max using intervals of 0.1
                               scMaximum.longevity = c(0), # fix max longevity at mean
                               sclitter_or_clutch_size_n = c(0),
                               scbroodvalue = c(0)), re_formula = NA)

# make plot
(minadult_survival_plot <- ggplot(survival_minadult_epred, aes(x = scAdult.survival, y = .epred)) +
  stat_lineribbon(color = "#238b45") + 
  scale_fill_manual(values = colorspace::lighten("#238b45", c(0.95, 0.75, 0.5))) + 
  guides(fill = "none") +
  labs(x = "Scaled Adult Survival", y = "Density Dependence Threshold") +
  theme_classic() +
  geom_point(data = dat_LH, aes(x= scAdult.survival, y = min_adult), pch = 19, color = "gray30") +
  coord_cartesian(ylim=c(0, 45)) + # ADDING LIMITS to Y-AXIS TO ZOOM IN 
  theme(axis.text.x = element_text(size=12), 
        axis.text.y = element_text(size=12), 
        axis.title.x = element_text(size=14), 
        axis.title.y = element_text(size=14)))

### make a fancy plot of clutch size for life history model using minimum adults (threshold) of DD as response variable

# find the min and max values for clutch size in the data
min(dat_LH$sclitter_or_clutch_size_n) # -1.35
max(dat_LH$sclitter_or_clutch_size_n)  # 3.39

# we also need the mean values of the other response variables
mean(dat_LH$scMaximum.longevity) # 0
mean(dat_LH$scAdult.survival) # 0
mean(dat_LH$scbroodvalue) # 0

# get predicted values for clutch size
clutch_minadult_epred <- minadult_LH_mod_nb %>% 
  epred_draws(newdata = tibble(sclitter_or_clutch_size_n = seq(-1.4, 3.4, 0.1), # create a sequence from the min to the max using intervals of 0.1
                               scMaximum.longevity = c(0), # fix at mean
                               scAdult.survival = c(0),
                               scbroodvalue = c(0)), re_formula = NA)

# make plot
(minadult_clutch_plot <-ggplot(clutch_minadult_epred, aes(x = sclitter_or_clutch_size_n, y = .epred)) +
  stat_lineribbon(color = "#56B4E9") + 
  scale_fill_manual(values = colorspace::lighten("#56B4E9", c(0.95, 0.75, 0.5))) + 
  guides(fill = "none") +
  labs(x = "Scaled Clutch Size", y = "Density Dependence Threshold") +
  theme_classic() +
  geom_point(data = dat_LH, aes(x= sclitter_or_clutch_size_n, y = min_adult), pch = 19, color = "gray30") +
  coord_cartesian(ylim=c(0, 20)) + # ADDING LIMITS to Y-AXIS TO ZOOM IN 
  theme(axis.text.x = element_text(size=12), 
        axis.text.y = element_text(size=12), 
        axis.title.x = element_text(size=14), 
        axis.title.y = element_text(size=14)))

### make a fancy plot of longevity for life history model using minimum adults (threshold) of DD as response variable

# find the min and max values for maximum longevity in the data
min(dat_LH$scMaximum.longevity) # -1.94
max(dat_LH$scMaximum.longevity)  # 3.6

# we also need the mean values of the other response variables
mean(dat_LH$sclitter_or_clutch_size_n) # 0
mean(dat_LH$scAdult.survival) # 0
mean(dat_LH$scbroodvalue) # 0

# get predicted values for clutch size
longevity_minadult_epred <- minadult_LH_mod_nb %>% 
  epred_draws(newdata = tibble(scMaximum.longevity = seq(-2, 3.6, 0.1), # create a sequence from the min to the max using intervals of 0.1
                               sclitter_or_clutch_size_n = c(0), # fix at mean
                               scAdult.survival = c(0),
                               scbroodvalue = c(0)), re_formula = NA)

# make plot
(minadult_longevity_plot <- ggplot(longevity_minadult_epred, aes(x = scMaximum.longevity, y = .epred)) +
  stat_lineribbon(color = "#F7710A") + 
  scale_fill_manual(values = colorspace::lighten("#F7710A", c(0.95, 0.75, 0.5))) + 
  guides(fill = "none") +
  labs(x = "Scaled Maximum Longevity", y = "Density Dependence Threshold") +
  theme_classic() +
  geom_point(data = dat_LH, aes(x= scMaximum.longevity, y = min_adult), pch = 19, color = "gray30") +
  coord_cartesian(ylim=c(0, 40)) + # ADDING LIMITS to Y-AXIS TO ZOOM IN 
  theme(axis.text.x = element_text(size=12), 
        axis.text.y = element_text(size=12), 
        axis.title.x = element_text(size=14), 
        axis.title.y = element_text(size=14))) 



##################################################################################
### Life History Model 2 ###

# Response variable: estimated intensity (slope) of density dependence
# Explanatory variables: longevity, clutch size, adult survival, brood value


# Take the absolute value of the slope estimate to make all the values positive 
#  Explore lognormal as possible distribution for the model
hist(log(abs(dat_LH$estimate))) # here's what that would look like


# Set intercept location at the log-median of the (positive) outcome
mu0_logn_LH <- log(median(abs(dat_LH$estimate), na.rm = TRUE))
mu0_logn_LH

# specify priors
priors_logn_LH <- c(
  # Fixed effects on log scale: ~N(0,1) ⇒ ~2.7× multiplicative change per 1 SD shift
  prior(normal(0, 1), class = "b"),
  
  # Intercept centered at log-median outcome, weak scale
  prior(student_t(3, -4.279, 2), class = "Intercept"), ## add mu0_lgn_LH by hand here
  
  # Lognormal residual SD (on log scale): weak half-t
  prior(student_t(3, 0, 0.5), class = "sigma"),
  
  # Random-intercept SDs (Species and phylogenetic effect via gr(BirdTree, cov = LH_cov))
  prior(student_t(3, 0, 0.5), class = "sd", group = "BirdTree")
)

priors_logn_LH

# Run lognormal model
# use scaled and centered variables in dat_LH
# for this model, increase adapt delta to 0.99 and tree depth to 15
# set init to "random" to reduce starting at extreme values

slope_LH_mod_logn <- brm(
  abs(estimate) ~ scMaximum.longevity + sclitter_or_clutch_size_n + scAdult.survival+ scbroodvalue + (1|gr(BirdTree, cov = LH_cov)),
  data = dat_LH,
  family = lognormal(link="identity"),
  data2 = list(LH_cov = LH_cov),
  control = list(adapt_delta = 0.999),
  save_pars = save_pars(all = TRUE),
  backend = "cmdstanr",
  init = "random",
  prior = priors_logn_LH, 
  chains = CHAINS, iter = ITER, warmup = WARMUP, thin = THIN, seed = BAYES_SEED
)  

# no divergent transitions


### Examine the model ###

summary(slope_LH_mod_logn)
tidy(slope_LH_mod_logn)
#  ESS and Rhat look fine

fixef(slope_LH_mod_logn, probs=c(0.05, 0.95)) # look at 90% CI

# positive relationship for clutch size. 95% CrI does not overlap zero
# positive relationship for adult survival. 90% CrI does not overlap (90% CI)


# trace and density plots
plot(slope_LH_mod_logn, nvariables = 4, ask = FALSE, theme=theme_minimal()) 
# posterior distribution of SD of tree still skewed, but this seems to reflect little phylogenetic signal


# posterior predictive check
pp_check(slope_LH_mod_logn) + xlim(0, 0.3)
pp_check(slope_LH_mod_logn, type="hist", ndraws = 11)
pp_check(slope_LH_mod_logn, type="ecdf_overlay") 
pp_check(slope_LH_mod_logn, type="intervals", ndraws = 100)

### Plot results ###

# adult survival
plot(conditional_effects(slope_LH_mod_logn, effects="scAdult.survival"), points = TRUE) 

# clutch size
plot(conditional_effects(slope_LH_mod_logn, effects="sclitter_or_clutch_size_n"), points = TRUE) 

### make a fancy plot of adult survival for life history model using intensity of DD as response variable

# find the min and max values for adult survival in the data
min(dat_LH$scAdult.survival) # -2.97
max(dat_LH$scAdult.survival)  # 3.55

# we also need the mean values of the other response variables
mean(dat_LH$scMaximum.longevity) # 0
mean(dat_LH$sclitter_or_clutch_size_n) # 0
mean(dat_LH$scbroodvalue) # 0

# get predicted values for adult survival
survival_slope_epred <- slope_LH_mod_logn %>% 
  epred_draws(newdata = tibble(scAdult.survival = seq(-3, 3.6, 0.1), # create a sequence from the min to the max using intervals of 0.1
                               scMaximum.longevity = c(0), # fix at mean
                               sclitter_or_clutch_size_n = c(0),
                               scbroodvalue = c(0)), re_formula = NA)

# make plot
(slope_survival_plot <- ggplot(survival_slope_epred, aes(x = scAdult.survival, y = .epred)) +
  stat_lineribbon(color = "#238b45") + 
  scale_fill_manual(values = colorspace::lighten("#238b45", c(0.95, 0.75, 0.5))) + 
  guides(fill = "none") +
  labs(x = "Scaled Adult Survival", y = "Intensity of Density Dependence)") +
  theme_classic() +
  geom_point(data = dat_LH, aes(x= scAdult.survival, y = abs(estimate)), pch = 19, color = "gray30") +
  coord_cartesian(ylim=c(0, 0.35)) + # ADDING LIMITS to Y-AXIS TO ZOOM IN 
  theme(axis.text.x = element_text(size=12), 
        axis.text.y = element_text(size=12), 
        axis.title.x = element_text(size=14), 
        axis.title.y = element_text(size=14)))

### make a fancy plot of clutch size for life history model using slope (intensity) of DD as response variable

# find the min and max values for clutch size in the data
min(dat_LH$sclitter_or_clutch_size_n) # -1.35
max(dat_LH$sclitter_or_clutch_size_n)  # 3.4

# we also need the mean values of the other response variables
mean(dat_LH$scMaximum.longevity) # 0
mean(dat_LH$scAdult.survival) # 0
mean(dat_LH$scbroodvalue) # 0

# get predicted values for clutch size
clutch_slope_epred <- slope_LH_mod_logn %>% 
  epred_draws(newdata = tibble(sclitter_or_clutch_size_n = seq(-1.4, 3.4, 0.1), # create a sequence from the min to the max using intervals of 0.1
                               scMaximum.longevity = c(0), # fix at mean
                               scAdult.survival = c(0),
                               scbroodvalue = c(0)), re_formula = NA)

# make plot
(slope_clutch_plot <- ggplot(clutch_slope_epred, aes(x = sclitter_or_clutch_size_n, y = .epred)) +
  stat_lineribbon(color = "#56B4E9") + 
  scale_fill_manual(values = colorspace::lighten("#56B4E9", c(0.95, 0.75, 0.5))) + 
  guides(fill = "none") +
  labs(x = "Scaled Clutch Size", y = "Intensity of Density Dependence") +
  theme_classic() +
  geom_point(data = dat_LH, aes(x= sclitter_or_clutch_size_n, y = abs(estimate)), pch = 19, color = "gray30") +
  coord_cartesian(ylim=c(0, 0.5)) + # ADDING LIMITS to Y-AXIS TO ZOOM IN 
  theme(axis.text.x = element_text(size=12), 
        axis.text.y = element_text(size=12), 
        axis.title.x = element_text(size=14), 
        axis.title.y = element_text(size=14)))


#########################################################################################
#################### Density Dependence Metrics and Morphometric Traits #################
#########################################################################################

### Join and prune tree and create covariance matrix ###

head(DD_morph)

# replace blanks in Genus species with underscore
dat_morph <- DD_morph %>%
  mutate( # make two columns with different names that are otherwise identical with underscore between genus and species
    rownames = str_replace(Species3_BirdTree, " ", "_"), # this one will be moved to rownames
    BirdTree = str_replace(Species3_BirdTree, " ", "_")) %>% # this one will remain in the data frame
  column_to_rownames(., var = "rownames")

head(dat_morph)

# link data to tree
morph_tree <- geiger::treedata(tree_out, dat_morph, sort=T)

plot(morph_tree$phy, cex=0.4)

# construct a covariance matrix on the relationship between species and scale
morph_cov <- ape::vcv.phylo(morph_tree$phy, corr =T)

# standardize predictors that will be used as fixed effects
dat_morph <- dat_morph |>
  dplyr::mutate(
    scHWI = as.numeric(scale(Hand.Wing.Index)),
    scLNMass= as.numeric(scale(log(Mass))),
    scBeak_PC1 = as.numeric(scale(Beak_PC1)),
    scBeak_PC2 = as.numeric(scale(Beak_PC2))
  )

##################################################################################
### Morphometrics Model 1 ###

# Response variable: number of adults at onset of DD effects (min_adult) 
# Explanatory variables: HWI, body mass, Beak PC1 and PC2
# Using negative binomial model again

# center Intercept near the typical count on the log scale
mu0_nb_Morph <- mean(dat_morph$min_adult, na.rm = TRUE)
# small guard against log(0)
log_mu0_nb_Morph <- log(ifelse(is.finite(mu0_nb_Morph) && mu0_nb_Morph > 0, mu0_nb_Morph, 1))

# specify priors
priors_nb_Morph <- c(
  # fixed effects (assumes predictors are standardized)
  set_prior("normal(0, 1)", class = "b"),
  
  # intercept centered on log(mean y) with robust tails
  set_prior(paste0("student_t(3, ", round(log_mu0_nb_Morph, 3), ", 2.5)"), class = "Intercept"),
  
  # phylogenetic random intercept SD (correlated by the tree)
  set_prior("student_t(3, 0, 1)",  class = "sd", group = "BirdTree"),
  
  # NB overdispersion (shape >= 0): weak preference for some overdispersion
  set_prior("exponential(0.5)", class = "shape")
)


# run NB model
minadult_Morph_mod_nb <- brm(
  min_adult ~ scHWI + scLNMass + scBeak_PC1 + scBeak_PC2 + (1|gr(BirdTree, cov = morph_cov)),
  data = dat_morph,
  family = negbinomial(),
  data2 = list(morph_cov = morph_cov),
  control = list(adapt_delta = 0.99), 
  backend = "cmdstanr",
  init = "random",
  save_pars = save_pars(all = TRUE),
  prior = priors_nb_Morph,
  chains = CHAINS, iter = ITER, warmup = WARMUP, thin = THIN, seed = BAYES_SEED
)

### Examine the model ###

summary(minadult_Morph_mod_nb) 
tidy(minadult_Morph_mod_nb)
# Rhat and ESS look okay

fixef(minadult_Morph_mod_nb, probs=c(0.05, 0.95)) # look at 90% CI
fixef(minadult_Morph_mod_nb, probs=c(0.075, 0.925)) # look at 85% CI

# negative trend with beak PC2. 85% CrI does not overlap zero

# trace and density plots
plot(minadult_Morph_mod_nb, nvariables = 4, ask = FALSE, theme=theme_minimal()) 
# these look ok... maybe more phylogenetic signal in this model than the life history one

# posterior predictive checks
pp_check(minadult_Morph_mod_nb)
pp_check(minadult_Morph_mod_nb, type="hist", ndraws = 11)
pp_check(minadult_Morph_mod_nb, type="ecdf_overlay")
pp_check(minadult_Morph_mod_nb, type="intervals", ndraws = 100)

### Plot results ###

# make plot of scBeak_PC2 (weak effect)
plot(conditional_effects(minadult_Morph_mod_nb, effects="scBeak_PC2"), points = TRUE) 

### make a fancy plot of beak PC2 for morphometric model using threshold (min # of adults) of DD as response variable

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
  labs(x = "Scaled Beak Shape (relative to size)", y = "Density Dependence Threshold") +
  theme_classic() +
  geom_point(data = dat_morph, aes(x= scBeak_PC2, y = min_adult), pch = 19, color = "gray30") +
  coord_cartesian(ylim=c(0, 20)) + # ADDING LIMITS to Y-AXIS TO ZOOM IN 
  theme(axis.text.x = element_text(size=12), 
        axis.text.y = element_text(size=12), 
        axis.title.x = element_text(size=14), 
        axis.title.y = element_text(size=14)))

##################################################################################
### Morphometric Model 2 ###
# Response variable: estimated intensity (slope) of density dependence
# Explanatory variables: HWI, body mass, Beak PC1 and PC2

# using lognormal model after taking the absolute value of the intensity (slope)

### scale all predictors
# Set intercept location at the log-median of the (positive) outcome
mu0_logn_Morph <- log(median(abs(dat_morph$estimate), na.rm = TRUE))
mu0_logn_Morph

# specify custom priors
priors_logn_Morph <- c(
  # Slopes: predictors are standardized, so Normal(0,1) is a good weakly informative prior
  prior(normal(0, 1), class = "b"),
  
  # Intercept centered at log-median outcome, weak scale
  prior(student_t(3, -4.279, 2), class = "Intercept"), ## add mu0_logn_Morph by hand here
  
  # Lognormal residual SD (on log scale): weak
  prior(student_t(3, 0, 1), class = "sigma"),
 
  # Phylogenetic (group-level) SD: weak
  prior(student_t(3, 0, 1), class = "sd", group = "BirdTree")
)


# Run lognormal model
# for this model, increased adapt delta and tree depth
# set init to "random" to reduce starting at extreme values

slope_Morph_mod_logn <- brm(
  abs(estimate) ~ scHWI + scLNMass + scBeak_PC1 + scBeak_PC2 + (1|gr(BirdTree, cov = morph_cov)),
  data = dat_morph,
  family = lognormal(link="identity"),
  data2 = list(morph_cov = morph_cov),
  control = list(adapt_delta = 0.999),
  backend = "cmdstanr",
  init = "random",
  save_pars = save_pars(all = TRUE),
  prior = priors_logn_Morph,
  chains = CHAINS, iter = ITER, warmup = WARMUP, thin = THIN, seed = BAYES_SEED
)  



###### examine the model ########
summary(slope_Morph_mod_logn)
# ESS and Rhat for fixed effects look fine
# ESS for sigma is somewhat low but still greater than 400 (4* # of chains)

fixef(slope_Morph_mod_logn, probs=c(0.05, 0.95)) # look at 90% CI
fixef(slope_Morph_mod_logn, probs=c(0.075, 0.925)) # look at 85% CI
# no notable relationships

# trace and density plots
plot(slope_Morph_mod_logn, nvariables = 4, ask = FALSE)

# posterior predictive checks
pp_check(slope_Morph_mod_logn) + xlim(0, 0.3)
pp_check(slope_Morph_mod_logn, type="hist", ndraws = 11)
pp_check(slope_Morph_mod_logn, type="ecdf_overlay")
pp_check(slope_Morph_mod_logn, type="intervals", ndraws = 100)

### Plot results ###

# no results to plot

#########################################################################################
#################### Density Dependence Metrics and Trophic Traits ######################
#########################################################################################

### Join and prune tree and create covariance matrix ###

head(DD_trophic)

# replace blanks in Genus species with underscore
dat_trophic <- DD_trophic %>%
  mutate( # make two columns with different names that are otherwise identical with underscore between genus and species
    rownames = str_replace(Species3_BirdTree, " ", "_"), # this one will be moved to rownames
    BirdTree = str_replace(Species3_BirdTree, " ", "_")) %>% # this one will remain in the data frame
  column_to_rownames(., var = "rownames")

dat_trophic$Species <- dat_trophic$BirdTree

#### link data to tree

trophic_tree <- geiger::treedata(tree_out, dat_trophic, sort=T)

plot(trophic_tree$phy, cex=0.2)

# construct a covariance matrix on the relationship between species
trophic_cov <- ape::vcv.phylo(trophic_tree$phy, corr = T)
# corr=T should make the diagonal = 1, which should help with posterior for random effect of tree

# standardize predictors that will be used as fixed effects
dat_trophic <- dat_trophic |>
  dplyr::mutate(
    scDiet = as.numeric(scale(Diet.score)),
    scStrata = as.numeric(scale(Strata.score)),
    scTrophicLevel = as.numeric(scale(TrophicLevel))
  )

##################################################################################
#### Trophic Model 1 #### 

# Response variable: adults at onset of DD effects (min_adult) 
# Explanatory variables:  trophic level, diet generalism and foraging strata generalism
# Using negative binomial model

# Center Intercept near the typical count on the log scale
mu0_nb_Trophic <- mean(dat_trophic$min_adult, na.rm = TRUE)
# small guard against log(0)
log_mu0_nb_Trophic <- log(ifelse(is.finite(mu0_nb_Trophic) && mu0_nb_Trophic > 0, mu0_nb_Trophic, 1))

# specify priors
priors_nb_Trophic <- c(
  # fixed effects (assumes predictors are standardized)
  set_prior("normal(0, 1)", class = "b"),
  
  # intercept centered on log(mean y) with robust tails
  set_prior(paste0("student_t(3, ", round(mu0_nb_Trophic, 3), ", 2.5)"), class = "Intercept"),
  
  # phylogenetic random intercept SD (correlated by the tree)
  set_prior("student_t(3, 0, 1)",  class = "sd", group = "BirdTree"),
  
  # NB overdispersion (shape >= 0): weak preference for some overdispersion
  set_prior("exponential(0.5)", class = "shape")
)


# run NB model
minadult_Trophic_mod_nb <- brm(
  min_adult ~ scDiet + scStrata + scTrophicLevel + (1|gr(BirdTree, cov = trophic_cov)),
  data = dat_trophic,
  family = negbinomial(),
  data2 = list(trophic_cov = trophic_cov),
  control = list(adapt_delta = 0.99), 
  backend = "cmdstanr",
  init = "random",
  save_pars = save_pars(all = TRUE),
  prior = priors_nb_Trophic,
  chains = CHAINS, iter = ITER, warmup = WARMUP, thin = THIN, seed = BAYES_SEED
)
# no divergent transitions 


### Examine the model ###

summary(minadult_Trophic_mod_nb) 
tidy(minadult_Trophic_mod_nb)
# Rhat and ESS are good

fixef(minadult_Trophic_mod_nb, probs=c(0.05, 0.95)) # look at 90% CI
fixef(minadult_Trophic_mod_nb, probs=c(0.075, 0.925)) # look at 85% CI

# negative relationship for trophic level. 95% CrI does not overlap zero

# trace and density plots
plot(minadult_Trophic_mod_nb, nvariables = 4, ask = FALSE) 

# posterior predictive checks
pp_check(minadult_Trophic_mod_nb)
pp_check(minadult_Trophic_mod_nb, type="hist", ndraws = 11) 
pp_check(minadult_Trophic_mod_nb, type="ecdf_overlay") 
pp_check(minadult_Trophic_mod_nb, type="intervals", ndraws = 100) 


### Plot results ###

# plot for Trophic Level
plot(conditional_effects(minadult_Trophic_mod_nb, effects="scTrophicLevel"), points = TRUE) 


### make a fancy plot of trophic level for trophic model using minimum adults (threshold) of DD as response variable

# find the min and max values for trophic level in the data
min(dat_trophic$scTrophicLevel) # -2.3
max(dat_trophic$scTrophicLevel)  # 1.37

# we also need the mean values of the other response variables
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
  coord_cartesian(ylim=c(0, 30)) + # ADDING LIMITS to Y-AXIS TO ZOOM IN 
  theme(axis.text.x = element_text(size=12), 
        axis.text.y = element_text(size=12), 
        axis.title.x = element_text(size=14), 
        axis.title.y = element_text(size=14)))

##################################################################################
### Trophic Model 2 ###

# Response variable: estimated intensity (slope) of density dependence
# Explanatory variables: trophic level, diet generalism and foraging strata generalism
# Using lognormal model


# Set intercept location at the log-median of the (positive) outcome
mu0_logn_Trophic <- log(median(abs(dat_trophic$estimate), na.rm = TRUE))
mu0_logn_Trophic

# specify priors
# to help with divergent transitions, tighten random effect SD and residual sigma from 1 to 0.5

priors_logn_Trophic <- c(
  
  # Slopes: predictors are standardized, so Normal(0,1) is a good weakly informative prior
  prior(normal(0, 1), class = "b"),
  
  # Intercept centered at log-median outcome, weak scale
  prior(student_t(3, -4.279, 2), class = "Intercept"), ## adding mu0_logn_Trophic by hand here
  
  # Lognormal residual SD (on log scale): weak, using 0.5 rather than 1 to decrease divergent transitions 
  prior(student_t(3, 0, 0.5), class = "sigma"),
  
  # Phylogenetic (group-level) SD: weak, using 0.5 rather than 1 to decrease divergent transitions
  prior(student_t(3, 0, 0.5), class = "sd", group = "BirdTree")
)


# run lognormal model
# for this model, increased adapt delta and tree depth
# set init to "random" to reduce starting at extreme values
slope_Trophic_mod_logn <- brm(
  abs(estimate) ~ scDiet + scStrata + scTrophicLevel + (1|gr(BirdTree, cov = trophic_cov)),
  data = dat_trophic,
  family = lognormal(link="identity"),
  data2 = list(trophic_cov = trophic_cov),
  control = list(adapt_delta = 0.999),
  backend = "cmdstanr",
  init = "random",
  save_pars = save_pars(all = TRUE),
  prior = priors_logn_Trophic, 
  chains = CHAINS, iter = ITER, warmup = WARMUP, thin = THIN, seed = BAYES_SEED
)  

# no divergent transitions

### Examine the model ###
summary(slope_Trophic_mod_logn)
#  ESS and Rhat look fine

fixef(slope_Trophic_mod_logn, probs=c(0.05, 0.95)) # look at 90% CI
fixef(slope_Trophic_mod_logn, probs=c(0.075, 0.925)) # look at 85% CI

# positive trend for strata. 85% CrI does not overlap zero

# trace and density plots
plot(slope_Trophic_mod_logn, nvariables = 4, ask = FALSE, theme=theme_minimal()) ## posterior distributions look better than any model so far... seem to have more phylogenetic signal as well

# posterior predictive check
pp_check(slope_Trophic_mod_logn) + xlim(0, 0.3)
pp_check(slope_Trophic_mod_logn, type="hist", ndraws = 11) 
pp_check(slope_Trophic_mod_logn, type="ecdf_overlay") 
pp_check(slope_Trophic_mod_logn, type="intervals", ndraws = 100)


### Plot results ###

# plot for Strata
plot(conditional_effects(slope_Trophic_mod_logn, effects="scStrata"), points = TRUE) 


### make a fancy plot of strata for trophic model using slope (intensity) of DD as response variable

# find the min and max values for strata in the data
min(dat_trophic$scStrata) # -1.82
max(dat_trophic$scStrata)  # 2.41

# we also need the mean values of the other response variables
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
    labs(x = "Scaled Trophic Strata", y = "Intensity of Density Dependence") +
    theme_classic() +
    geom_point(data = dat_trophic, aes(x= scStrata, y = abs(estimate)), pch = 19, color = "gray30") +
    coord_cartesian(ylim=c(0, 0.25)) + # ADDING LIMITS to Y-AXIS TO ZOOM IN 
    theme(axis.text.x = element_text(size=12), 
          axis.text.y = element_text(size=12), 
          axis.title.x = element_text(size=14), 
          axis.title.y = element_text(size=14)))


#########################################################################################
########### Density Dependence Metrics and Sexual Selection Traits ######################
#########################################################################################

### Join and prune tree and create covariance matrix ###

# replace blanks in Genus species with underscore
dat_ss <- DD_SS %>%
  mutate( # make two columns with different names that are otherwise identical with underscore between genus and species
    rownames = str_replace(Species3_BirdTree, " ", "_"), # this one will be moved to rownames
    BirdTree = str_replace(Species3_BirdTree, " ", "_")) %>% # this one will remain in the data frame
  column_to_rownames(., var = "rownames") %>%
  filter(!is.na(Plumage_DC)) # two species are missing values for DC, drop them
  

head(dat_ss)

dat_ss$Species <- dat_ss$BirdTree

# link data to tree

ss_tree <- geiger::treedata(tree_out, dat_ss, sort=T)

plot(ss_tree$phy, cex=0.4)

# construct a covariance matrix on the relationship between species
# corr=T should make the diagonal = 1, which should help with posterior for random effect of tree.
ss_cov <- ape::vcv.phylo(ss_tree$phy, corr = T)

# standardize predictors that will be used as fixed effects
dat_ss <- dat_ss %>%
  dplyr::mutate(
    scDM  = as.numeric(scale(Wing_DM)),
    scDC = as.numeric(scale(Plumage_DC)),
    sc_ssM = as.numeric(scale(sex.sel.m)),
    sc_ssF = as.numeric(scale(sex.sel.f)))


##################################################################################
### Sexual Selection Model 1 ### 

# Response variable: number of adults at onset of DD effects (min_adult)
# Explanatory variables: Plumage_DC, Wing_DM and Mating system
# Using negative binomial model


# center Intercept near the typical count on the log scale
mu0_nb_SS <- mean(dat_ss$min_adult, na.rm = TRUE)
# small guard against log(0)
log_mu0_nb_SS <- log(ifelse(is.finite(mu0_nb_SS) && mu0_nb_SS > 0, mu0_nb_SS, 1))

# specify priors
priors_nb_SS <- c(
  # fixed effects (assumes predictors are standardized)
  set_prior("normal(0, 1)", class = "b"),
  
  # intercept centered on log(mean y) with robust tails
  set_prior(paste0("student_t(3, ", round(log_mu0_nb_SS, 3), ", 2.5)"), class = "Intercept"),
  
  # phylogenetic random intercept SD (correlated by the tree)
  set_prior("student_t(3, 0, 1)",  class = "sd", group = "BirdTree"),
  
  # NB overdispersion (shape >= 0): weak preference for some overdispersion
  set_prior("exponential(0.5)", class = "shape")
)



# run NB model

minadult_SS_mod_nb <- brm(
  min_adult ~ scDC + scDM + sc_ssM + sc_ssF  + (1|gr(BirdTree, cov = ss_cov)),
  data = dat_ss,
  family = negbinomial(),
  data2 = list(ss_cov = ss_cov),
  control = list(adapt_delta = 0.999), 
  save_pars = save_pars(all = TRUE),
  prior = priors_nb_SS,
  backend = "cmdstanr",
  init = "random",
  chains = CHAINS, iter = ITER, warmup = WARMUP, thin = THIN, seed = BAYES_SEED
)

# no divergent transitions
# reminder: two species were removed from the analysis due to lack of dichromatism - downy woodpecker and red-breasted sapsucker

### Examine the model ###
summary(minadult_SS_mod_nb) 
tidy(minadult_SS_mod_nb)
# ESS and Rhat look good

fixef(minadult_SS_mod_nb, probs=c(0.05, 0.95)) # look at 90% CI
fixef(minadult_SS_mod_nb, probs=c(0.075, 0.925)) # look at 85% CI

# positive relationship for sexual dichromatism. 90% CrI does not overlap zero

# trace and density plots
plot(minadult_SS_mod_nb, nvariables = 4, ask = FALSE, , theme=theme_minimal()) 

# posterior predictive checks
pp_check(minadult_SS_mod_nb)
pp_check(minadult_SS_mod_nb, type="hist", ndraws = 11)
pp_check(minadult_SS_mod_nb, type="ecdf_overlay") 
pp_check(minadult_SS_mod_nb, type="intervals", ndraws = 100) 

### Plot results ###

plot(conditional_effects(minadult_SS_mod_nb, effects="scDC"), points = TRUE) 


### make a fancy plot of dichromatism for sexual selection model using minimum adults (threshold) of DD as response variable

# find the min and max values for sexual dichromatism in the data
min(dat_ss$scDC) # -1.4
max(dat_ss$scDC)  # 2.62

# we also need the mean values of the other response variables
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
  labs(x = "Scaled Sexual Dichromatism", y = "Density Dependence Threshold") +
  theme_classic() +
  geom_point(data = dat_ss, aes(x= scDC, y = min_adult), pch = 19, color = "gray30") +
  coord_cartesian(ylim=c(0, 20)) + # ADDING LIMITS to Y-AXIS TO ZOOM IN 
  theme(axis.text.x = element_text(size=12), 
        axis.text.y = element_text(size=12), 
        axis.title.x = element_text(size=14), 
        axis.title.y = element_text(size=14)))

##################################################################################
### Sexual Selection Model 2 ###

# Response variable: estimated intensity (slope) of density dependence
# Explanatory variables: Plumage_DC, Wing_DM and Mating system
# Using lognormal model

# Set intercept location at the log-median of the (positive) outcome
mu0_logn_SS <- log(median(abs(dat_ss$estimate), na.rm = TRUE))

# specify priors
# to help with divergent transitions, tightening random effect SD and residual sigma from 1 to 0.5

priors_logn_SS <- c(
  # Slopes: predictors are standardized, so Normal(0,1) is a good weakly informative prior
  prior(normal(0, 1), class = "b"),
  
  # Intercept centered at log-median outcome, weak scale
  prior(student_t(3, -4.134, 2), class = "Intercept"), ## adding mu0_logn_SS by hand here as -4.134
  
  # Lognormal residual SD (on log scale): weak, using 0.5 rather than 1 to decrease divergent transitions 
  prior(student_t(3, 0, 0.5), class = "sigma"),
 
   # Phylogenetic (group-level) SD: weak, using 0.5 rather than 1 to decrease divergent transitions
  prior(student_t(3, 0, 0.5), class = "sd", group = "BirdTree")
)


# run lognormal model
# increased adapt delta and tree depth to help with divergent transitions
# set init to "random" to reduce starting at extreme values

slope_SS_mod_logn <- brm(
  abs(estimate) ~ scDC + scDM + sc_ssM + sc_ssF + (1|gr(BirdTree, cov = ss_cov)),
  data = dat_ss,
  family = lognormal(link="identity"),
  data2 = list(ss_cov = ss_cov),
  control = list(adapt_delta = 0.999),
  backend = "cmdstanr",
  init = "random",
  save_pars = save_pars(all = TRUE),
  prior = priors_logn_SS, 
  chains = CHAINS, iter = ITER, warmup = WARMUP, thin = THIN, seed = BAYES_SEED
)  

# no divergent transitions
# same two species were excluded as previous model

### Examine the model ###

summary(slope_SS_mod_logn)
tidy(slope_SS_mod_logn)
#  ESS and Rhat look fine

fixef(slope_SS_mod_logn, probs=c(0.05, 0.95))
fixef(slope_SS_mod_logn, probs=c(0.075, 0.925))

# positive relationship for ssM. 95% CrI does not overlap zero
# negative trend for dichromatism. 85% CrI does not overlap zero

# trace and density plots
plot(slope_SS_mod_logn, nvariables = 4, ask = FALSE,theme=theme_minimal()) 

# posterior predictive checks
pp_check(slope_SS_mod_logn) + xlim(0, 0.3)
pp_check(slope_SS_mod_logn, type="hist", ndraws = 11) 
pp_check(slope_SS_mod_logn, type="ecdf_overlay") 
pp_check(slope_SS_mod_logn, type="intervals", ndraws = 100) 


### Plot results ###

# plots for ssM and DC
plot(conditional_effects(slope_SS_mod_logn, effects="scDC"), points = TRUE) 
plot(conditional_effects(slope_SS_mod_logn, effects="sc_ssM"), points = TRUE)

### make a fancy plot of dichromatism for sexual selection model using intensity (slope) of DD as response variable

# find the min and max values for sexual dichromatism in the data
min(dat_ss$scDC) # -1.4
max(dat_ss$scDC)  # 2.62

# we also need the mean values of the other response variables
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
  labs(x = "Scaled Sexual Dichromatism", y = "Intensity of Density Dependence") +
  theme_classic() +
  geom_point(data = dat_ss, aes(x= scDC, y = abs(estimate)), pch = 19, color = "gray30") +
  coord_cartesian(ylim=c(0, 0.1)) + # ADDING LIMITS to Y-AXIS TO ZOOM IN 
  theme(axis.text.x = element_text(size=12), 
        axis.text.y = element_text(size=12), 
        axis.title.x = element_text(size=14), 
        axis.title.y = element_text(size=14)))

### make a fancy plot of male sexual selection using intensity (slope) of DD as response variable

# find the min and max values for male sexual selection in the data
min(dat_ss$sc_ssM) # -0.84
max(dat_ss$sc_ssM)  # 2.77

# we also need the mean values of the other response variables
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
  labs(x = "Scaled Male Sexual Selection", y = "Intensity of Density Dependence") +
  theme_classic() +
  geom_point(data = dat_ss, aes(x= sc_ssM, y = abs(estimate)), pch = 19, color = "gray30") +
  coord_cartesian(ylim=c(0, 0.1)) + # ADDING LIMITS to Y-AXIS TO ZOOM IN 
  theme(axis.text.x = element_text(size=12), 
        axis.text.y = element_text(size=12), 
        axis.title.x = element_text(size=14), 
        axis.title.y = element_text(size=14)))


##########################################################################
##########################################################################
##########################################################################


# function to tidy a single brms model summary into one data frame 
.tidy_brms_summary <- function(model, model_name = deparse(substitute(model))) {
  # Decide what we have
  is_fit <- inherits(model, "brmsfit")
  smry   <- if (is_fit) summary(model) else model
  
  # Safe family string
  get_family_str <- function(x) {
    if (inherits(x, "brmsfit")) {
      f <- brms::family(x)
      return(paste0(f$family, "_", f$link))
    }
    # summary list may or may not hold a family object/char
    if (is.list(x) && !is.null(x$family)) {
      f <- x$family
      if (is.list(f) && !is.null(f$family) && !is.null(f$link))
        return(paste0(f$family, "_", f$link))
      if (is.character(f)) return(paste(f, collapse = "_"))
    }
    NA_character_
  }
  
  fam_str <- get_family_str(model)
  
  # number of observations
  nobs_val <- tryCatch(
    if (is_fit) stats::nobs(model) else NA_integer_,
    error = function(e) NA_integer_
  )
  
  # Fixed effects (always present in brms summary)
  fixed_tbl <-
    tibble::as_tibble(smry$fixed, rownames = "term") |>
    dplyr::mutate(component = "fixed", group = NA_character_)
  
  # Random-effect SD/cor pieces (may be absent)
  random_tbl <-
    if (!is.null(smry$random) && length(smry$random)) {
      purrr::imap_dfr(smry$random, function(mat, grp) {
        tibble::as_tibble(mat, rownames = "term") |>
          dplyr::mutate(component = "random", group = grp)
      })
    } else {
      tibble::tibble(term = character(),
                     Estimate = numeric(), `Est.Error` = numeric(),
                     `l-95% CI` = numeric(), `u-95% CI` = numeric(),
                     Rhat = numeric(), Bulk_ESS = numeric(), Tail_ESS = numeric(),
                     component = character(), group = character())
    }
  
  out <- dplyr::bind_rows(fixed_tbl, random_tbl) |>
    dplyr::mutate(model = model_name,
                  family = fam_str,
                  nobs = nobs_val) |>
    dplyr:: mutate(across(where(is.numeric), ~round(.x, 3))) |> # ROUND ALL NUMERIC COLUMNS TO HAVE 3 DIGITS AFTER DECIMAL
    dplyr::relocate(model, family, nobs, component, group, term)
  
  # Ensure expected columns exist
  for (col in c("Estimate", "Est.Error", "l-95% CI", "u-95% CI",
                "Rhat", "Bulk_ESS", "Tail_ESS")) {
    if (!col %in% names(out)) out[[col]] <- NA_real_
  }
  
  out
}

# write CSVs for multiple brms models
write_brms_summaries <- function(models, out_dir = "Results", write_combined = TRUE,
                                 combined_filename = "all_traitmodels_summary.csv") {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # If user passed a bare list of models, name them automatically
  if (is.null(names(models)) || any(names(models) == "")) {
    names(models) <- paste0("model_", seq_along(models))
  }
  
  tidy_list <- imap(models, ~ .tidy_brms_summary(.x, .y))
  
  # Write one CSV per model
  walk(tidy_list, function(df) {
    fn <- file.path(out_dir, paste0(unique(df$model), ".csv"))
    write_csv(df, fn)
  })
  
  # Optionally write combined CSV
  combined <- bind_rows(tidy_list) %>%
    arrange(model, component, group, term)
  
  if (write_combined) {
    write_csv(combined, file.path(out_dir, combined_filename))
  }
  
  invisible(combined)
}

modList <- list(
            SlopeLH = summary(slope_LH_mod_logn), 
            SlopeMorph = summary(slope_Morph_mod_logn),
            SlopeTrophic = summary(slope_Trophic_mod_logn),
            SlopeSS = summary(slope_SS_mod_logn),
             MinAdultLH = summary(minadult_LH_mod_nb),
             MinAdultMorph = summary(minadult_Morph_mod_nb),
             MinAdultTrophic = summary(minadult_Trophic_mod_nb),
             MinAdultSS = summary(minadult_SS_mod_nb))

# write files
write_brms_summaries(modList)

# save all models
saveRDS(slope_LH_mod_logn, here("Models/TraitModels", "slope_LH_mod_logn.rds"))
saveRDS(minadult_LH_mod_nb, here("Models/TraitModels", "minadult_LH_mod_nb.rds"))
saveRDS(slope_Morph_mod_logn, here("Models/TraitModels", "slope_Morph_mod_logn.rds"))
saveRDS(minadult_Morph_mod_nb, here("Models/TraitModels", "minadult_Morph_mod_nb.rds"))
saveRDS(slope_Trophic_mod_logn, here("Models/TraitModels", "slope_Trophic_mod_logn.rds"))
saveRDS(minadult_Trophic_mod_nb, here("Models/TraitModels", "minadult_Trophic_mod_nb.rds"))
saveRDS(slope_SS_mod_logn, here("Models/TraitModels", "slope_SS_mod_logn.rds"))
saveRDS(minadult_SS_mod_nb, here("Models/TraitModels", "minadult_SS_mod_nb.rds"))
