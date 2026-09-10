###### MAPS Project: Negative Feedback #######
### Script name: Step3_PriorSensitivityTesting.R
### Author(s): SLJ

########### Objective/Description of Script #####################
# testing the two sets of custom priors using prior sensitivity testing from the priorsense package
# the overarching goal is to settle on one set of priors for the models that seems to work decently across all species
#################################################################

### Setup ###

# load packages
library(tidyverse)
library(tidybayes)
library(broom)
library(broom.mixed) 
library(performance)
library(bayesplot) 
library(here)
library(brms)
library(priorsense)
library(scales)
library(ggforce)
library(MetBrewer)
library(gridExtra)
library(marginaleffects)
library(loo)
library(cmdstanr)

#######################################
# model options
options(
  mc.cores = 4,  # distribute the 4 chains across 4 cores
  brms.threads = 4,
  brms.backend = "cmdstanr",
  brms.file_refit = "on_change"
)

# model specifications
CHAINS <- 4
ITER <- 5000
WARMUP <- 2000
THIN <- 2
BAYES_SEED <- 456

######################################################################################
# specify two groups of priors that we came up with in previous script
priors_1 <- c(prior(student_t(3, 0, 1), class = Intercept, dpar="hu"), # intercept of hurdle model
              prior(student_t(3,-0.1, 1), class = b, dpar="hu"), # slope of the hurdle model
              prior(normal(-0.1, 1), class=Intercept), # intercept of the lognormal model
              prior(normal(-0.1, 1), class=b), # slope of the lognormal model
              prior(student_t(3, 0, 2.5), class=sd, dpar="hu"), # group-level standard deviation for hurdle model -> same as default
              prior(student_t(3, 0, 2.5), class=sd), # group-level standard deviation for lognormal model -> same as default
              prior(student_t(3, 0, 2.5), class=sigma), # sigma (residual standard deviation) for lognormal model -> same as default
              prior(lkj(1), class = cor)) # correlation between random slopes and intercepts -> same as default


priors_2 <- c(prior(student_t(3, -0.4, 1), class = Intercept, dpar="hu"), # intercept of hurdle model
              prior(student_t(3,-0.2, 1), class = b, dpar="hu"), # slope of the hurdle model
              prior(normal(-0.2, 1), class=Intercept), # intercept of the lognormal model
              prior(normal(-0.2, 1), class=b), # slope of the lognormal model
              prior(student_t(3, 0, 2.5), class=sd, dpar="hu"), # group-level standard deviation for hurdle model -> same as default
              prior(student_t(3, 0, 2.5), class=sd), # group-level standard deviation for lognormal model -> same as default
              prior(student_t(3, 0, 2.5), class=sigma), # sigma (residual standard deviation) for lognormal model -> same as default
              prior(lkj(1), class = cor)) # correlation between random slopes and intercepts -> same as default

###################################################################################
# import data 
NF_breed_dat <- readRDS(here("Outputs", "NF_breed_dat.rds"))

######################################################################################

# run models using priors_1
m_allspp_p1 <- NF_breed_dat %>%
  group_by(SPEC) %>% # group each species into its own data frame
  nest() %>%
  mutate(dat1 = map(data, ~as_tibble(droplevels(.x))), # make version of species data frame where unused factor levels for STA and year are dropped
         model = map(dat1, # use this df in the model
                  ~ brms::brm(
                    bf(FY_to_A  ~ Adult + (1 + Adult|STA) + (1|year), 
                       hu ~ Adult + (1 + Adult|STA) + (1|year)), 
                    data = .,
                    family = hurdle_lognormal(),
                    prior = priors_1,
                    save_pars = save_pars(all = TRUE),
                    control = list(adapt_delta = 0.99),
                    chains = CHAINS, iter = ITER, warmup = WARMUP, thin = THIN, seed = BAYES_SEED)))

saveRDS(m_allspp_p1, here("Models", "m_allspp_priors1.rds"))

# use powerscale sensitivity test from priorsense package on each model
sensitivity_p1 <- m_allspp_p1 %>%
  mutate(
    summary = map(model, ~summary), # add model summary
    sensitivity = map(model, ~priorsense::powerscale_sensitivity(.x)))

# extract prior sensitivity results to examine them
# first, make a list of priors that we want to examine
priorslist<- c("b_Intercept", "b_Adult", "b_hu_Intercept", "b_hu_Adult", 
               "sd_STA__Intercept", "sd_STA__Adult", "sd_year__Intercept",
               "sd_STA__hu_Intercept", "sd_STA__hu_Adult", "sd_year__hu_Intercept",
               "cor_STA__hu_Intercept__hu_Adult", "cor_STA__Intercept__Adult", "sigma")

priorsensitivity_p1 <- sensitivity_p1 %>% unnest(sensitivity) %>% 
  select(-data, -model, -summary) %>%
  filter(variable %in% priorslist) %>%
  filter(diagnosis != "-")

saveRDS(priorsensitivity_p1, here("Outputs", "priors1_priorsensitivity.rds"))

# look at conflicts
# there are some prior-data conflicts, which suggests the priors may be too strong
conflict_p1 <- priorsensitivity_p1 %>% filter(diagnosis == "potential prior-data conflict")
conflict_p1count <- conflict_p1 %>% group_by(variable) %>% count() %>% rename(n_p1 = n)

######################################################################################
# run models using priors_2
m_allspp_p2 <- NF_breed_dat %>%
  group_by(SPEC) %>% # group each species into its own data frame
  nest() %>%
  mutate(dat1 = map(data, ~as_tibble(droplevels(.x))), # make version of species data frame where unused factor levels for STA and year are dropped
         model = map(dat1, # use this df in the model
                  ~ brms::brm(
                    bf(FY_to_A  ~ Adult + (1 + Adult|STA) + (1|year), 
                       hu ~ Adult + (1 + Adult|STA) + (1|year)), 
                    data = .,
                    family = hurdle_lognormal(),
                    prior = priors_2,
                    save_pars = save_pars(all = TRUE),
                    control = list(adapt_delta = 0.99),
                    chains = CHAINS, iter = ITER, warmup = WARMUP, thin = THIN, seed = BAYES_SEED)))

saveRDS(m_allspp_p2, here("Models", "m_allspp_priors2.rds"))

# use powerscale sensitivity test from priorsense package on each model
sensitivity_p2 <- m_allspp_p2 %>%
  mutate(
    summary = map(model, ~summary), # add model summary
    sensitivity = map(model, ~priorsense::powerscale_sensitivity(.x)))

# extract prior sensitivity results to examine them
# first, make a list of priors that we want to examine
priorslist <- c("b_Intercept", "b_Adult", "b_hu_Intercept", "b_hu_Adult", 
               "sd_STA__Intercept", "sd_STA__Adult", "sd_year__Intercept",
               "sd_STA__hu_Intercept", "sd_STA__hu_Adult", "sd_year__hu_Intercept",
               "cor_STA__hu_Intercept__hu_Adult", "cor_STA__Intercept__Adult", "sigma")

priorsensitivity_p2 <- sensitivity_p2 %>% unnest(sensitivity) %>% 
  select(-data, -model, -summary) %>%
  filter(variable %in% priorslist) %>%
  filter(diagnosis != "-")

saveRDS(priorsensitivity_p2, here("Outputs", "priors2_priorsensitivity.rds"))

# look at conflicts
# very similar to priors_1, although slightly better
# there are still a decent number of conflicts for the hurdle model slope (fixed effect)
conflict_p2 <- priorsensitivity_p2 %>% filter(diagnosis == "potential prior-data conflict")
conflict_p2count <- conflict_p2 %>% group_by(variable) %>% count() %>% rename(n_p2 = n)


# to investigate further, we will run a batch of models with default flat prior for slope of Adult for hurdle model
# we are making this prior as weak as possible to see if the prior-data conflict warnings are resolved
# all other priors in this set as the same as priors_2 
priors_3 <- c(prior(student_t(3, -0.4, 1), class = Intercept, dpar="hu"), # intercept of hurdle model
              prior(normal(-0.2, 1), class=Intercept), # intercept of the lognormal model
              prior(normal(-0.2, 1), class=b), # slope of the lognormal model
              prior(student_t(3, 0, 2.5), class=sd, dpar="hu"), # group-level standard deviation for hurdle model -> same as default
              prior(student_t(3, 0, 2.5), class=sd), # group-level standard deviation for lognormal model -> same as default
              prior(student_t(3, 0, 2.5), class=sigma), # sigma (residual standard deviation) for lognormal model -> same as default
              prior(lkj(1), class = cor)) # correlation between random slopes and intercepts -> same as default

m_allspp_p3 <- NF_breed_dat %>%
  group_by(SPEC) %>% # group each species into its own data frame
  nest() %>%
  mutate(dat1 = map(data, ~as_tibble(droplevels(.x))), # make version of species data frame where unused factor levels for STA and year are dropped
         model = map(dat1, # use this df in the model
                  ~ brms::brm(
                    bf(FY_to_A  ~ Adult + (1 + Adult|STA) + (1|year), 
                       hu ~ Adult + (1 + Adult|STA) + (1|year)), 
                    data = .,
                    family = hurdle_lognormal(),
                    prior = priors_3,
                    save_pars = save_pars(all = TRUE),
                    control = list(adapt_delta = 0.99),
                    chains = CHAINS, iter = ITER, warmup = WARMUP, thin = THIN, seed = BAYES_SEED)))

saveRDS(m_allspp_p3, here("Models", "m_allspp_priors3.rds"))

# use powerscale sensitivity test from priorsense package on each model
sensitivity_p3 <- m_allspp_p3 %>%
  mutate(
    summary = map(model, ~summary), # add model summary
    sensitivity = map(model, ~priorsense::powerscale_sensitivity(.x)))

# extract prior sensitivity results to examine them
# first, make a list of priors that we want to examine
priorslist<- c("b_Intercept", "b_Adult", "b_hu_Intercept", "b_hu_Adult", 
               "sd_STA__Intercept", "sd_STA__Adult", "sd_year__Intercept",
               "sd_STA__hu_Intercept", "sd_STA__hu_Adult", "sd_year__hu_Intercept",
               "cor_STA__hu_Intercept__hu_Adult", "cor_STA__Intercept__Adult", "sigma")

priorsensitivity_p3 <- sensitivity_p3 %>% unnest(sensitivity) %>% 
  select(-data, -model, -summary) %>%
  filter(variable %in% priorslist) %>%
  filter(diagnosis != "-")

saveRDS(priorsensitivity_p3, here("Outputs", "priors3_priorsensitivity.rds"))

# look at conflicts
conflict_p3 <- priorsensitivity_p3 %>% filter(diagnosis == "potential prior-data conflict")
conflict_p3count <- conflict_p3 %>% group_by(variable) %>% count() %>% rename(n_p3 = n)

# the conflicts remain even with the flat slope
# we will proceed with the priors_2 set of priors for all models

rm(m_allspp_p1)
rm(m_allspp_p2)
rm(m_allspp_p3)
