###### MAPS Project: Density Dependence #######
### Script name: Step4_ModelComparison.R
### Author(s): SLJ

########### Objective/Description of Script #####################
# run 3 models for density dependence with 1 model per species:
# m1: random intercepts for STA and year
# m2: random intercepts and slopes for STA and random intercepts for year
# m3: random intercepts and slopes for both STA and year
# overarching goal: compare these models across all species 
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
# set some colors for plotting
clrs <- MetBrewer::met.brewer("Juarez")
clrs

# model options
options(
  mc.cores = 4,  # distribute the 4 chains across 4 cores
  brms.threads = 4,
  brms.backend = "cmdstanr",
  brms.file_refit = "on_change"
)

# model specifications
CHAINS <- 4
ITER <- 9000
WARMUP <- 3000
THIN <- 3
BAYES_SEED <- 289
# increased iterations, warmup and thin above defaults to reduce autocorrelation and increase ESS

##################################################################################
# specify custom priors (determined in previous scripts)

# for model 1
custompriors_a <- c(prior(student_t(3, -0.4, 1), class = Intercept, dpar="hu"), # intercept of hurdle model
                    prior(student_t(3,-0.2, 1), class = b, dpar="hu"), # slope of the hurdle model
                    prior(normal(-0.2, 1), class=Intercept), # intercept of the lognormal model
                    prior(normal(-0.2, 1), class=b), # slope of the lognormal model
                    prior(student_t(3, 0, 2.5), class=sd, dpar="hu"), # group-level standard deviation for hurdle model -> same as default
                    prior(student_t(3, 0, 2.5), class=sd), # group-level standard deviation for lognormal model -> same as default
                    prior(student_t(3, 0, 2.5), class=sigma)) # sigma (residual standard deviation) for lognormal model -> same as default

# for models 2 and 3 where we need a prior for the correlation between random slopes and intercepts
custompriors_b <- c(prior(student_t(3, -0.4, 1), class = Intercept, dpar="hu"), # intercept of hurdle model
                    prior(student_t(3,-0.2, 1), class = b, dpar="hu"), # slope of the hurdle model
                    prior(normal(-0.2, 1), class=Intercept), # intercept of the lognormal model
                    prior(normal(-0.2, 1), class=b), # slope of the lognormal model
                    prior(student_t(3, 0, 2.5), class=sd, dpar="hu"), # group-level standard deviation for hurdle model -> same as default
                    prior(student_t(3, 0, 2.5), class=sd), # group-level standard deviation for lognormal model -> same as default
                    prior(student_t(3, 0, 2.5), class=sigma), # sigma (residual standard deviation) for lognormal model -> same as default
                    prior(lkj(1), class = cor)) # correlation between random slopes and intercepts -> same as default

########################################################
########### Get Species-specific data ##################

# import data 
DD_breed_dat <- readRDS(here("Outputs", "DD_breed_dat.rds"))

######################## Run the models #############################

# we are using the loo package to perform PSIS-LOO cross validation on each model after it has run

###### Model 1 ########

# run m1 for all species
m1_allspp_loo <- DD_breed_dat %>%
  group_by(SPEC) %>% # group each species into its own data frame
  nest() %>%
  mutate(dat1 = map(data, ~as_tibble(droplevels(.x))), # make version of species data frame where unused factor levels for STA and year are dropped
         m1 = map( dat1, # use this df in the model
                   ~ brms::brm(
                     bf(FY_to_A  ~ Adult + (1|STA) + (1|year), 
                        hu ~ Adult + (1|STA) + (1|year)), 
                     data = .,
                     family = hurdle_lognormal(),
                     prior = custompriors_a,
                     save_pars = save_pars(all = TRUE),
                     control = list(adapt_delta=0.99),
                     chains = CHAINS, iter = ITER, warmup = WARMUP, thin = THIN, seed = BAYES_SEED)),
         loo_m1 = map(m1, ~loo(.x, save_psis = T)))

saveRDS(m1_allspp_loo, here("Models", "m1_allspp_loo.rds"))   
rm(m1_allspp_loo) # remove to free up space

###### Model 2 ########

# run m2 for all species
m2_allspp_loo <- DD_breed_dat %>%
  group_by(SPEC) %>% # group each species into its own data frame
  nest() %>%
  mutate(dat1 = map(data, ~as_tibble(droplevels(.x))), # make version of species data frame where unused factor levels for STA and year are dropped
         m2 = map(dat1, # use this df in the model
                  ~ brms::brm(
                    bf(FY_to_A  ~ Adult + (1 + Adult|STA) + (1|year), 
                       hu ~ Adult + (1 + Adult|STA) + (1|year)), 
                    data = .,
                    family = hurdle_lognormal(),
                    prior = custompriors_b,
                    save_pars = save_pars(all = TRUE),
                    control = list(adapt_delta=0.99),
                    chains = CHAINS, iter = ITER, warmup = WARMUP, thin = THIN, seed = BAYES_SEED)),
         loo_m2 = map(m2, ~loo(.x, save_psis = T)))

saveRDS(m2_allspp_loo, here("Models", "m2_allspp_loo.rds"))
rm(m2_allspp_loo) # remove to free up space

###### Model 3 ########

# run m3 for all species
m3_allspp_loo <- DD_breed_dat %>%
  group_by(SPEC) %>% # group each species into its own data frame
  nest() %>%
  mutate(dat1 = map(data, ~as_tibble(droplevels(.x))), # make version of species data frame where unused factor levels for STA and year are dropped
         m3 = map( dat1, # use this df in the model
                   ~ brms::brm(
                     bf(FY_to_A  ~ Adult + (1 + Adult|STA) + (1 + Adult|year), 
                        hu ~ Adult + (1 + Adult|STA) + (1 + Adult|year)), 
                     data = .,
                     family = hurdle_lognormal(),
                     prior = custompriors_b,
                     save_pars = save_pars(all = TRUE),
                     control = list(adapt_delta=0.99),
                     chains = CHAINS, iter = ITER, warmup = WARMUP, thin = THIN, seed = BAYES_SEED)),
         loo_m3 = map(m3, ~loo(.x, save_psis = T)))

saveRDS(m3_allspp_loo, here("Models", "m3_allspp_loo.rds")) 
rm(m3_allspp_loo) # remove to free up space

######################################################################################
########## Use loo package for model evaluation ################

###### Model 1 ########
# import m1 with loo results
m1_allspp_loo <- readRDS(here("Models", "m1_allspp_loo.rds"))

# examine the results to see how many models have high Pareto K values
# for models with ok/good Pareto K values, we can extract ELPD and its standard error
# for models with high/poor Pareto K values, we need to try some other approaches to get reliable ELPD and SE estimates

psis_loo_m1 <- m1_allspp_loo %>% pull(loo_m1, name=SPEC) # this pulls out the results of loo
psis_loo_m1 # we want to simplify this to only have the number of good, bad and very bad for each species

class(psis_loo_m1) # this is a list
# to get the pareto K values for the first species in the list, we would use this code:
psis_loo_m1$BHGR$diagnostics$pareto_k
pluck(psis_loo_m1, 1, 3, 1) # the same thing can be achieved using the pluck function

# we can also get the estimates for the first species in the list, including elpd and its SE using this:
psis_loo_m1$BHGR$estimates
pluck(psis_loo_m1, 1, 1, 1) # elpd
pluck(psis_loo_m1, 1, 1, 4) # SE of elpd

rm(psis_loo_m1) # free up space by removing this from the environment

# use pluck function to add the pareto k values for each species as a column in the data frame
# create a summary of the pareto k values for each species
m1_paretok_summary <- m1_allspp_loo %>% 
  mutate(pareto_k = map(loo_m1, ~ pluck(., 3, 1))) %>%
  pull(pareto_k, name = SPEC) %>% # pull pareto k values into a list where each item in the list is named by species and contains a vector of pareto k values
  enframe(., value = "paretoK") %>% # force the list into a data frame and put the vector of pareto k values into a column called paretoK
  mutate( # use map and mutate on the vector of paretoK values for each species
    total = map(paretoK, length), # get the total number of data points for the species
    m1_good = map(paretoK, ~ sum(. <=0.7)),  # find the number of data points for each species that have good pareto k values (< or = 0.7)
    m1_bad = map(paretoK, ~ sum(. >0.7))) %>% # find the number of data points for each species that have bad pareto k values (>0.7)
select(-paretoK) %>% # drop the vector of pareto k values
rename(SPEC = name)
# get values from loo for species that are OK (no high pareto K values)

# vector of species that are OK
m1_goodPK <- m1_paretok_summary %>% filter(m1_bad == 0) %>% pull(SPEC)

# get elpd and SE for species that are OK
ELPD_m1_goodPK <-  m1_allspp_loo %>% 
  mutate(elpd_loo = map(loo_m1, ~ pluck(., 1, 1)), # extract elpd
         SE_elpd = map(loo_m1, ~ pluck(., 1, 4))) %>% # extract SE for elpd
  select(SPEC, elpd_loo, SE_elpd) %>%
  filter(SPEC %in% m1_goodPK) # use vector from previous step to filter results

# identify species with problematic PK values to try moment matching
m1_highPK <- m1_paretok_summary %>% 
  filter(m1_bad > 0) %>% pull(SPEC)

# use loo with moment matching
m1_momentmatch <- m1_allspp_loo %>% select(SPEC, m1) %>%
  filter(SPEC %in% m1_highPK) %>%
  mutate(loo_mm_m1 = map(m1, ~loo(.x, moment_match = T)))
# this works for some but not all species 

saveRDS(m1_momentmatch, here("Models", "m1_momentmatch.rds"))

rm(m1_allspp_loo) # free up memory by removing this large object

# get a summary of the results for this batch of species
m1_mm_paretok_summary <- m1_momentmatch %>% 
  mutate(pareto_k = map(loo_mm_m1, ~ pluck(., 3, 1))) %>%
  pull(pareto_k, name = SPEC) %>% # pull pareto k values into a list where each item in the list is named by species and contains a vector of pareto k values
  enframe(., value = "paretoK") %>% # force the list into a data frame and put the vector of pareto k values into a column called paretoK
  mutate( # use map and mutate on the vector of paretoK values for each species
    total = map(paretoK, length), # get the total number of data points for the species
    m1_good = map(paretoK, ~ sum(. <=0.7)),  # find the number of data points for each species that have good pareto k values (< or = 0.7)
    m1_bad = map(paretoK, ~ sum(. >0.7 & .<=1))) %>% # find the number of data points for each species that have bad pareto k values (>0.7)
  select(-paretoK) %>% # drop the vector of pareto k values
  rename(SPEC = name)

# vector of species that are OK
m1_mm_goodPK <- m1_mm_paretok_summary %>% filter(m1_bad == 0) %>% pull(SPEC)

# vector of species that still have high PK values
m1_mm_highPK <- m1_mm_paretok_summary %>% filter(m1_bad > 0) %>% pull(SPEC)

# get elpd and SE for species where moment matching worked
ELPD_m1_mm_goodPK <-  m1_momentmatch %>% 
  mutate(elpd_loo = map(loo_mm_m1, ~ pluck(., 1, 1)), # extract elpd
         SE_elpd = map(loo_mm_m1, ~ pluck(., 1, 4))) %>% # extract SE for elpd
  select(SPEC, elpd_loo, SE_elpd) %>%
  filter(SPEC %in% m1_mm_goodPK) # use vector from previous step to filter results

# finally, use reloo() for species that still have high PK values
# this will compute exact cross-validation for problematic observations where approximate LOO-CV may return incorrect results
m1_to_reloo <- m1_momentmatch %>% 
  select(SPEC, m1, loo_mm_m1) %>%
  filter(SPEC %in% m1_mm_highPK) 

rm(m1_momentmatch) # free up memory before using reloo function

options(future.globals.maxSize = +Inf)
m1_reloo <- m1_to_reloo %>%
  mutate(reloo_m1 = map2(m1, loo_mm_m1, ~reloo(.x, .y)))

# save
saveRDS(m1_reloo, here("Models", "m1_reloo.rds"))

options(future.globals.maxSize = 500 *1024^2) # return to default of 500MB

# get elpd and SE for models that needed reloo
ELPD_m1_reloo_goodPK <-  m1_reloo %>% 
  mutate(elpd_loo = map(reloo_m1, ~ pluck(., 1, 1)), # extract elpd
         SE_elpd = map(reloo_m1, ~ pluck(., 1, 4))) %>% # extract SE for elpd
  select(SPEC, elpd_loo, SE_elpd)

# combine elpd and SE for all species for model 1
ELPD_m1_allspp <- bind_rows(ELPD_m1_goodPK, # species that were good with just loo()
                            ELPD_m1_mm_goodPK, # species where moment matching helped
                            ELPD_m1_reloo_goodPK) %>% # species that needed reloo()
  rename(m1_elpd = elpd_loo, m1_SE = SE_elpd)

nrow(ELPD_m1_allspp) # should be 62

# save into Outputs folder
saveRDS(ELPD_m1_allspp, here("Outputs", "ELPD_m1_allspp.rds"))

# remove object from the environment
rm(m1_reloo)


###### Model 2 ########
### repeat for model 2

# import m2 with loo results
m2_allspp_loo <- readRDS(here("Models", "m2_allspp_loo.rds"))

# use pluck function to add the pareto k values for each species as a column in the data frame
# create a summary of the pareto k values for each species
m2_paretok_summary <- m2_allspp_loo %>% 
  mutate(pareto_k = map(loo_m2, ~ pluck(., 3, 1))) %>%
  pull(pareto_k, name = SPEC) %>% # pull pareto k values into a list where each item in the list is named by species and contains a vector of pareto k values
  enframe(., value = "paretoK") %>% # force the list into a data frame and put the vector of pareto k values into a column called paretoK
  mutate( # use map and mutate on the vector of paretoK values for each species
    total = map(paretoK, length), # get the total number of data points for the species
    m2_good = map(paretoK, ~ sum(. <=0.7)),  # find the number of data points for each species that have good pareto k values (< or = 0.7)
    m2_bad = map(paretoK, ~ sum(. >0.7))) %>% # find the number of data points for each species that have bad pareto k values (>0.7)
  select(-paretoK) %>% # drop the vector of pareto k values
  rename(SPEC = name)
# get values from loo for species that are OK (no high pareto K values)

# vector of species that are OK
m2_goodPK <- m2_paretok_summary %>% filter(m2_bad == 0) %>% pull(SPEC)

# get elpd and SE for species that are OK
ELPD_m2_goodPK <-  m2_allspp_loo %>% 
  mutate(elpd_loo = map(loo_m2, ~ pluck(., 1, 1)), # extract elpd
         SE_elpd = map(loo_m2, ~ pluck(., 1, 4))) %>% # extract SE for elpd
  select(SPEC, elpd_loo, SE_elpd) %>%
  filter(SPEC %in% m2_goodPK) # use vector from previous step to filter results

# identify species with problematic PK values to try moment matching
m2_highPK <- m2_paretok_summary %>% 
  filter(m2_bad > 0) %>% pull(SPEC)

# use loo with moment matching
m2_to_momentmatch <- m2_allspp_loo %>% select(SPEC, m2) %>%
  filter(SPEC %in% m2_highPK) 
rm(m2_allspp_loo) # free up memory by removing this large object

m2_momentmatch <- m2_to_momentmatch %>%
  mutate(loo_mm_m2 = map(m2, ~loo(.x, moment_match = T)))
# this works for some but not all species 

saveRDS(m2_momentmatch, here("Models", "m2_momentmatch.rds"))

# get a summary of the results for this batch of species
m2_mm_paretok_summary <- m2_momentmatch %>% 
  mutate(pareto_k = map(loo_mm_m2, ~ pluck(., 3, 1))) %>%
  pull(pareto_k, name = SPEC) %>% # pull pareto k values into a list where each item in the list is named by species and contains a vector of pareto k values
  enframe(., value = "paretoK") %>% # force the list into a data frame and put the vector of pareto k values into a column called paretoK
  mutate( # use map and mutate on the vector of paretoK values for each species
    total = map(paretoK, length), # get the total number of data points for the species
    m2_good = map(paretoK, ~ sum(. <=0.7)),  # find the number of data points for each species that have good pareto k values (< or = 0.7)
    m2_bad = map(paretoK, ~ sum(. >0.7 & .<=1))) %>% # find the number of data points for each species that have bad pareto k values (>0.7)
  select(-paretoK) %>% # drop the vector of pareto k values
  rename(SPEC = name)

# vector of species that are OK
m2_mm_goodPK <- m2_mm_paretok_summary %>% filter(m2_bad == 0) %>% pull(SPEC)

# vector of species that still have high PK values
m2_mm_highPK <- m2_mm_paretok_summary %>% filter(m2_bad > 0) %>% pull(SPEC)

# get elpd and SE for species where moment matching worked
ELPD_m2_mm_goodPK <-  m2_momentmatch %>% 
  mutate(elpd_loo = map(loo_mm_m2, ~ pluck(., 1, 1)), # extract elpd
         SE_elpd = map(loo_mm_m2, ~ pluck(., 1, 4))) %>% # extract SE for elpd
  select(SPEC, elpd_loo, SE_elpd) %>%
  filter(SPEC %in% m2_mm_goodPK) # use vector from previous step to filter results

# finally, use reloo() for species that still have high PK values
# this will compute exact cross-validation for problematic observations where approximate LOO-CV may return incorrect results
m2_to_reloo <- m2_momentmatch %>% 
  select(SPEC, m2, loo_mm_m2) %>%
  filter(SPEC %in% m2_mm_highPK) 

rm(m2_momentmatch) # free up memory before using reloo function

options(future.globals.maxSize = +Inf)
m2_reloo <- m2_to_reloo %>%
  mutate(reloo_m2 = map2(m2, loo_mm_m2, ~reloo(.x, .y)))

# save
saveRDS(m2_reloo, here("Models", "m2_reloo.rds"))

options(future.globals.maxSize = 500 *1024^2) # return to default of 500MB

# get elpd and SE for models that needed reloo
ELPD_m2_reloo_goodPK <-  m2_reloo %>% 
  mutate(elpd_loo = map(reloo_m2, ~ pluck(., 1, 1)), # extract elpd
         SE_elpd = map(reloo_m2, ~ pluck(., 1, 4))) %>% # extract SE for elpd
  select(SPEC, elpd_loo, SE_elpd)

# combine elpd and SE for all species for model 2
ELPD_m2_allspp <- bind_rows(ELPD_m2_goodPK, # species that were good with just loo()
                            ELPD_m2_mm_goodPK, # species where moment matching helped
                            ELPD_m2_reloo_goodPK) %>% # species that needed reloo()
  rename(m2_elpd = elpd_loo, m2_SE = SE_elpd)

nrow(ELPD_m2_allspp) # should be 62

# save into Outputs folder
saveRDS(ELPD_m2_allspp, here("Outputs", "ELPD_m2_allspp.rds"))

# remove object from the environment
rm(m2_reloo)


###### Model 3 ########
### repeat for model 3

# import m1 with loo results
m3_allspp_loo <- readRDS(here("Models", "m3_allspp_loo.rds"))

# use pluck function to add the pareto k values for each species as a column in the data frame
# create a summary of the pareto k values for each species
m3_paretok_summary <- m3_allspp_loo %>% 
  mutate(pareto_k = map(loo_m3, ~ pluck(., 3, 1))) %>%
  pull(pareto_k, name = SPEC) %>% # pull pareto k values into a list where each item in the list is named by species and contains a vector of pareto k values
  enframe(., value = "paretoK") %>% # force the list into a data frame and put the vector of pareto k values into a column called paretoK
  mutate( # use map and mutate on the vector of paretoK values for each species
    total = map(paretoK, length), # get the total number of data points for the species
    m3_good = map(paretoK, ~ sum(. <=0.7)),  # find the number of data points for each species that have good pareto k values (< or = 0.7)
    m3_bad = map(paretoK, ~ sum(. >0.7))) %>% # find the number of data points for each species that have bad pareto k values (>0.7)
  select(-paretoK) %>% # drop the vector of pareto k values
  rename(SPEC = name)
# get values from loo for species that are OK (no high pareto K values)

# vector of species that are OK
m3_goodPK <- m3_paretok_summary %>% filter(m3_bad == 0) %>% pull(SPEC)

# get elpd and SE for species that are OK
ELPD_m3_goodPK <-  m3_allspp_loo %>% 
  mutate(elpd_loo = map(loo_m3, ~ pluck(., 1, 1)), # extract elpd
         SE_elpd = map(loo_m3, ~ pluck(., 1, 4))) %>% # extract SE for elpd
  select(SPEC, elpd_loo, SE_elpd) %>%
  filter(SPEC %in% m3_goodPK) # use vector from previous step to filter results

# identify species with problematic PK values to try moment matching
m3_highPK <- m3_paretok_summary %>% 
  filter(m3_bad > 0) %>% pull(SPEC)

# use loo with moment matching
m3_to_momentmatch <- m3_allspp_loo %>% select(SPEC, m3) %>%
  filter(SPEC %in% m3_highPK) 
rm(m3_allspp_loo) # free up memory by removing this large object

m3_momentmatch <- m3_to_momentmatch %>%
  mutate(loo_mm_m3 = map(m3, ~loo(.x, moment_match = T)))

saveRDS(m3_momentmatch, here("Models", "m3_momentmatch.rds"))

# get a summary of the results for this batch of species
m3_mm_paretok_summary <- m3_momentmatch %>% 
  mutate(pareto_k = map(loo_mm_m3, ~ pluck(., 3, 1))) %>%
  pull(pareto_k, name = SPEC) %>% # pull pareto k values into a list where each item in the list is named by species and contains a vector of pareto k values
  enframe(., value = "paretoK") %>% # force the list into a data frame and put the vector of pareto k values into a column called paretoK
  mutate( # use map and mutate on the vector of paretoK values for each species
    total = map(paretoK, length), # get the total number of data points for the species
    m3_good = map(paretoK, ~ sum(. <=0.7)),  # find the number of data points for each species that have good pareto k values (< or = 0.7)
    m3_bad = map(paretoK, ~ sum(. >0.7 & .<=1))) %>% # find the number of data points for each species that have bad pareto k values (>0.7)
  select(-paretoK) %>% # drop the vector of pareto k values
  rename(SPEC = name)

# vector of species that are OK
m3_mm_goodPK <- m3_mm_paretok_summary %>% filter(m3_bad == 0) %>% pull(SPEC)

# vector of species that still have high PK values
m3_mm_highPK <- m3_mm_paretok_summary %>% filter(m3_bad > 0) %>% pull(SPEC)

# get elpd and SE for species where moment matching worked
ELPD_m3_mm_goodPK <-  m3_momentmatch %>% 
  mutate(elpd_loo = map(loo_mm_m3, ~ pluck(., 1, 1)), # extract elpd
         SE_elpd = map(loo_mm_m3, ~ pluck(., 1, 4))) %>% # extract SE for elpd
  select(SPEC, elpd_loo, SE_elpd) %>%
  filter(SPEC %in% m3_mm_goodPK) # use vector from previous step to filter results

# finally, use reloo() for species that still have high PK values
# this will compute exact cross-validation for problematic observations where approximate LOO-CV may return incorrect results
m3_to_reloo <- m3_momentmatch %>% 
  select(SPEC, m3, loo_mm_m3) %>%
  filter(SPEC %in% m3_mm_highPK) 

rm(m3_momentmatch) # free up memory before using reloo function

options(future.globals.maxSize = +Inf) # 32GB
m3_reloo <- m3_to_reloo %>%
  mutate(reloo_m3 = map2(m3, loo_mm_m3, ~reloo(.x, .y)))

# save
saveRDS(m3_reloo, here("Models", "m3_reloo.rds"))

options(future.globals.maxSize = 500 *1024^2) # return to default of 500MB

# get elpd and SE for models that needed reloo
ELPD_m3_reloo_goodPK <-  m3_reloo %>% 
  mutate(elpd_loo = map(reloo_m3, ~ pluck(., 1, 1)), # extract elpd
         SE_elpd = map(reloo_m3, ~ pluck(., 1, 4))) %>% # extract SE for elpd
  select(SPEC, elpd_loo, SE_elpd)

# combine elpd and SE for all species for model 2
ELPD_m3_allspp <- bind_rows(ELPD_m3_goodPK, # species that were good with just loo()
                            ELPD_m3_mm_goodPK, # species where moment matching helped
                            ELPD_m3_reloo_goodPK) %>% # species that needed reloo()
  rename(m3_elpd = elpd_loo, m3_SE = SE_elpd)

nrow(ELPD_m3_allspp) # should be 62

# save into Outputs folder
saveRDS(ELPD_m3_allspp, here("Outputs", "ELPD_m3_allspp.rds"))

# remove object from the environment
rm(m3_reloo)

########### combine the 3 summaries ##########

compare_paretoK <- left_join(ELPD_m1_allspp, ELPD_m2_allspp) %>%
  left_join(., ELPD_m3_allspp)

glimpse(compare_paretoK) # some columns are lists
View(compare_paretoK)
# it appears that for most species, m2 is superior to m1, but m3 is often equivalent to m2 and isn't offering additional advantages

saveRDS(compare_paretoK, here("Outputs", "compare_paretoK_m1m2m3.rds")) 

# change list columns to numeric to enable writing to csv
compare_paretoK_int <- compare_paretoK %>%
  mutate(across(m1_elpd:m3_SE, as.numeric))
glimpse(compare_paretoK_int)

write.csv(compare_paretoK_int, here("Outputs", "compare_paretoK_m1m2m3.csv"))

#### Conclusion: Use Model 2 with random slopes and intercepts for Station and random intercepts for Year
