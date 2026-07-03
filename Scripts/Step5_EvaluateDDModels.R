###### MAPS Project: Density Dependence #######
### Script name: Step5_EvaluateDDModels.R
### Author(s): SLJ

########### Objective/Description of Script #####################

#################################################################

### Setup ###

# load packages
library(tidyverse)
library(tidybayes)
library(broom)
library(broom.mixed) 
library(confintr)
library(performance)
library(bayesplot) 
library(here)
library(brms)
library(scales)
library(ggforce)
library(MetBrewer)
library(gridExtra)
library(marginaleffects)
library(loo)

#####################################################################
# set some colors for plotting
clrs <- MetBrewer::met.brewer("Juarez")

#####################################################################

# import models from previous script

m2_allspp_loo <- readRDS(here("Models", "m2_allspp_loo.rds"))

# to examine any of the models individually, use this example:
# this extracts the model for HOWA (hooded warbler)
m2_allspp_loo %>% filter(SPEC == "HOWA") %>% pull(m2) %>% pluck(1) %>% summary()  # get model summary
m2_allspp_loo %>% filter(SPEC == "HOWA") %>% pull(m2) %>% pluck(1) %>% plot() # examine plots

# make an object that only contains the models
m2_onlymodels <- m2_allspp_loo %>% select(SPEC, m2)
rm(m2_allspp_loo)

# get model summary, pull fixed effects, random effects and sigma from each model summary across all species
modelresults_m2 <- m2_onlymodels %>%
  mutate(
    summary = map(m2, ~summary), # add model summary
    fixed = map(m2, ~ { summary(.x)$fixed %>% rownames_to_column(., "Parameter")}), # pull the fixed effects for each model
    random = map(m2, ~ {summary(.x)$random %>% bind_rows(., .id ="Parameter") %>% # pull the random effects for each model
        rownames_to_column(., "Intercept") %>% mutate(Parameter = paste0(Parameter, "_", str_sub(Intercept, end = -5))) %>%
        dplyr::select(-Intercept)}),
    sigma = map(m2, ~ {summary(.x)$spec_pars %>% bind_rows(., .id="Parameter")}) # pull sigma for each model
  ) %>%
  select(-m2)

# extract fixed effects, random effects, and sigma from nested data frame
fixed_m2 <- modelresults_m2 %>% unnest(fixed) %>% select(-summary, -random, -sigma)
random_m2 <- modelresults_m2 %>% unnest(random) %>% select(-summary, -fixed, -sigma)
sigma_m2 <- modelresults_m2 %>% unnest(sigma) %>% select(-summary, -random, -fixed) %>%
  mutate(Parameter = paste0("sigma"))

       
# combine into new data frame and add rounding criteria
all_effects_m2 <- bind_rows(fixed_m2, random_m2, sigma_m2) %>% 
  arrange(SPEC) %>%
  mutate(across(Estimate:Rhat, ~round(. ,4))) %>% # round estimate, SE, 95% CI boundaries and Rhat to 4 decimal places
  mutate(across(Bulk_ESS:Tail_ESS, ~ round(. ,0))) # round bulk and tail ESS to nearest integer

write.csv(all_effects_m2, here("Outputs", "DDModelEffects_m2.csv"))


######################################################################
################ Examine Model Convergence #########################

####### Trace and Density Plots ##########

color_scheme_set("brightblue")
color_scheme_view()

tracedensity_m2 <- m2_onlymodels %>%
  mutate(tracedensity = map(m2, ~ plot(.x, theme=theme_minimal(),
                                       nvariables = 5, # put 5 variables/plots per page
                                       plot = T, newpage = T, ask= F)), .keep="none") %>% # .keep = "none" drops everything but the newly created plots
  deframe() # this converts the plots into a list

# export all plots as a single pdf and label each page with model/species
pdf(here("Figures", "tracedensityplots_m2.pdf"), onefile = TRUE)
for (i in seq(length(tracedensity_m2))) {
  for(j in 1:3) { # there are two pages of plots for each species
    do.call("grid.arrange", c(tracedensity_m2[[i]][[j]], top = names(tracedensity_m2)[[i]]))   
  }
}
dev.off()   

####### Visualize the same thing using Violin Plots ##########
violin_m2 <- m2_onlymodels %>%
  mutate(violin = map(m2, ~ mcmc_plot(.x, type = "violin")), .keep="none") %>% # violin plots
  deframe() # this converts the plots into a list

# export all violin plots as a single pdf and label each page with model/species
pdf(here("Figures", "violinplots_m2.pdf"), onefile = TRUE)
for (i in seq(length(violin_m2))) {
  do.call("grid.arrange", c(violin_m2[i], top = names(violin_m2)[i]))   
}
dev.off()

####### Examine Autocorrelation Plots ##########

autocorr_m2 <- m2_onlymodels %>%
  mutate(autocorr = map(m2, ~ mcmc_plot(.x, type = "acf")), .keep="none") %>% # acf = autocorrelation plots
  deframe() # this converts the plots into a list

# export all autocorrelation plots as a single pdf and label each page with model/species
pdf(here("Figures", "autocorrelationplots_m2.pdf"), onefile = TRUE)
for (i in seq(length(autocorr_m2))) {
  do.call("grid.arrange", c(autocorr_m2[i], top = names(autocorr_m2)[i]))   
}
dev.off()

# remove large files from the environment
rm(tracedensity_m2)
rm(violin_m2)
rm(autocorr_m2)

######################################################################
################ Posterior Predictive Checks #########################

########### Density plots ################

ppc_density_m2 <- m2_onlymodels %>%
  mutate(density = map(m2, ~ pp_check(.x, ndraws = 100) + xlim(-0.5, 20)), .keep="none") %>%
  deframe()

# export all plots as a single pdf and label each page with model/species
pdf(here("Figures", "ppc_density_plots_m2.pdf"), onefile = TRUE)
for (i in seq(length(ppc_density_m2))) {
  do.call("grid.arrange", c(ppc_density_m2[i], top = names(ppc_density_m2)[i]))   
}
dev.off()

############# Histograms #################

ppc_hist_m2 <- m2_onlymodels %>%
  mutate(hist = map(m2, ~ pp_check(.x, type="hist", ndraws = 11) + xlim(-0.5, 20)), .keep="none") %>%
  deframe()

# export all plots as a single pdf and label each page with model/species
pdf(here("Figures", "ppc_histograms_m2.pdf"), onefile = TRUE)
for (i in seq(length(ppc_hist_m2))) {
  do.call("grid.arrange", c(ppc_hist_m2[i], top = names(ppc_hist_m2)[i]))   
}
dev.off()

############### Empirical Cumulative Distribution Function ##################

ppc_ecdf_m2 <- m2_onlymodels %>%
  mutate(ecdf = map(m2, ~ pp_check(.x, type="ecdf_overlay", ndraws = 100)), .keep="none") %>% 
  deframe()

# export all plots as a single pdf and label each page with model/species
pdf(here("Figures", "ppc_ecdfoverlay_m2.pdf"), onefile = TRUE)
for (i in seq(length(ppc_ecdf_m2))) {
  do.call("grid.arrange", c(ppc_ecdf_m2[i], top = names(ppc_ecdf_m2)[i]))   
}
dev.off()

############### Interval Plots ##################

ppc_interval_m2 <- m2_onlymodels %>%
  mutate(interval = map(m2, ~ pp_check(.x, type="intervals", ndraws = 100)), .keep="none") %>% 
  deframe()

# export all plots as a single pdf and label each page with model/species

pdf(here("Figures", "ppc_intervalplots_m2.pdf"), onefile = TRUE)
for (i in seq(length(ppc_interval_m2))) {
  do.call("grid.arrange", c(ppc_interval_m2[i], top = names(ppc_interval_m2)[i]))   
}
dev.off()

# remove large files from the environment
rm(ppc_density_m2)
rm(ppc_ecdf_m2)
rm(ppc_hist_m2)
rm(ppc_interval_m2)

####################################################################
################### Plot Conditional Effects #######################

# to make the plots, we need to get the max number of Adults observed for each species

# to do this, we need to import the data for the models
DD_breed_dat <- readRDS(here("Outputs", "DD_breed_dat.rds"))

# Now, find the max # of adults for each species
max <- DD_breed_dat %>% 
  group_by(SPEC) %>%
  summarize(maxAdult = max(Adult))

epred_m2 <- m2_onlymodels %>%
  left_join(., max) %>%
  mutate(
    epred_draws = map2(m2, maxAdult, # map over both the model and the max number of adults for each species
                       ~epred_draws(newdata = tibble(Adult=seq(0, .y, 1)), object =.x, re_formula = NA)) 
  )

# put results into a data frame to make them easier to manipulate
epred_df_m2 <- epred_m2 %>% pull(epred_draws, name=SPEC) %>% bind_rows(., .id="SPEC")


# create plots of the conditional effects of Adult abundance on productivity
# we want one plot per species
CEplots_m2 <- epred_df_m2 %>% 
  ggplot(aes(x = Adult, y = .epred)) +
  stat_lineribbon(color = clrs[2]) +
  scale_fill_manual(values = colorspace::lighten(clrs[2], c(0.95, 0.75, 0.5))) +
  guides(fill = "none") +
  labs(x = "Adult Abundance", y = "Predicted Productivity") +
  theme_classic() 

# print onto multiple pages 
# two columns and three rows of plots for each page
CEplots_m2 + facet_wrap_paginate(~SPEC, scales="free", nrow = 3, ncol = 2, page = 1) # view page 1

# export plots in a single pdf
pdf(here("Figures", "conditionaleffectsplots_m2.pdf"), onefile = TRUE)
for(i in 1:11){ 
  print(CEplots_m2 + facet_wrap_paginate(~SPEC, scales="free", nrow = 3 , ncol = 2, page = i))
}
dev.off()

# remove from the environment
rm(CEplots_m2)

###################### Calculate Threshold & Intensity ##########################
# get the value of the slope across the range of values for Adults
slopes_m2 <- m2_onlymodels %>%
  left_join(., max) %>%
  mutate(
    slopes = map2(m2, maxAdult, # map over both the model and the max number of adults for each species
                  ~slopes(.x,  # use slopes() in marginaleffects package
                          newdata = datagrid(Adult=seq(0, .y, 1)), # sequence along from zero to max number of adults for each species
                          variables = "Adult", 
                          type = "response",  # this incorporates both hu and mu models and back transforms the values to be in original scale of variables
                          re_formula = NA)),
    post_draws = map(slopes, ~ posterior_draws(.x)) # extract the posterior draws from the object produced by slopes()
  ) %>%
  select(-m2)

# put results into a data frame to make them easier to manipulate
slopedraws_m2_df <- slopes_m2 %>% 
  select(SPEC, post_draws) %>%
  pull(post_draws, name=SPEC) %>% 
  bind_rows(., .id="SPEC")

# classify the slope as either positive, negative or zero (flat) using the confidence interval
slopesign_m2 <- slopedraws_m2_df %>% 
  select(SPEC, estimate, conf.low, conf.high, Adult) %>% 
  distinct() %>%
  mutate(slope_sign = factor(if_else( conf.low > 0 & conf.high >0, "positive",
                                      if_else( conf.low < 0 & conf.high < 0, "negative", "zero")), levels = c("positive", "zero", "negative")))

# For each species, find the range of Adults where slope is positive, negative or zero
rangeadults_m2 <- slopesign_m2 %>% 
  group_by(SPEC, slope_sign) %>%
  summarize(min_adult = min(Adult),
            max_adult = max(Adult)) 

# we want to have a data frame with positive, negative and zero for every species, even if they don't have slopes that meet these criteria
slope_summary_m2 <- rangeadults_m2 %>% 
  ungroup() %>%
  expand(SPEC, slope_sign) %>% # this will use range adults to create every combination of species and slope sign
  left_join(., rangeadults_m2) # join the max and min values from range adults to the expanded data frame
# there will be NAs in rows where a species does not have that type of slope

saveRDS(slope_summary_m2, here("Outputs", "SlopeSummarybySpecies_m2.rds"))
write.csv(slope_summary_m2, here("Outputs", "SlopeSummarybySpecies_m2.csv"))

# keep all the values of adult abundance where the slope is negative
negative_m2 <- slope_summary_m2 %>% 
  filter(slope_sign == "negative")

ave_slope_m2 <-  m2_onlymodels %>%
  left_join(., negative_m2, by = "SPEC") %>%
  mutate(
    ave_slope = pmap(list(m2, min_adult, max_adult),  # use pmap to map over 3 inputs. They need to be put into a list
                     \(mod, min, max) # name the 3 arguments for the function 
                     avg_slopes(mod,  # use the avg_slopes function from marginal effects package. First, give the function the models
                                newdata = datagrid(Adult = seq(min, max, by=1)), # use the min and max values of adults for each species
                                variables = "Adult",
                                type = "response",
                                re_formula = NA)),
    post_draws = map(ave_slope, ~ posterior_draws(.x)) # extract the posterior draws from the object produced by slopes()
  ) 

ave_negative_slope_m2 <- ave_slope_m2 %>% 
  unnest(ave_slope) %>% # unnest the results
  dplyr::select(-m2, -post_draws) # remove columns that are not needed

saveRDS(ave_negative_slope_m2, here("Outputs", "AverageNegativeSlopebySpecies_m2.rds"))
write.csv(ave_negative_slope_m2, here("Outputs", "AverageNegativeSlopebySpecies_m2.csv"))


############# Examine correlation between threshold and intensity #################

confintr::ci_cor(ave_negative_slope_m2$min_adult, ave_negative_slope_m2$estimate, method="spearman", type="bootstrap")

############# Measure phylogenetic signal in threshold and intensity #################

# load additional packages
library(ape)
library(geiger)
library(ggtree)
library(ggtreeExtra)

# import data
DD_breed_dat <- readRDS(here("Outputs", "DD_breed_dat.rds"))
head(DD_breed_dat)

# import bird list used by MAPS to convert 4-letter bird codes to scientific names
birdcodes <- read.csv(here("Data", "IBP-AOS-LIST23.csv"), header=T) 
head(birdcodes)  
# downloaded from IBP website: https://www.birdpop.org/pages/birdSpeciesCodes.php

# import conversion sheet from avonet (Tobias et al. 2022) that aligns bird tree, bird life and eBird taxonomies 
nameconvert <- read.csv(here("Data", "BirdNamesConversion.csv"), header=T)
colnames(nameconvert)

# import phylogenetic tree
tree_out <- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

# get a list of the 62 species with number of data points for each
DD_SPEC <- DD_breed_dat %>%
  group_by(SPEC) %>% # group the data by species
  summarize(n = n()) %>% # get the number of data points for each species
  arrange(desc(n)) # arrange is descending order

# combine names with species codes
names <- DD_SPEC %>% 
  left_join(., birdcodes, by="SPEC") %>% # left join keeps all the rows in DD_breed_dat and adds anything that matches from birdcodes
  dplyr::select(-SP, -CONF, -SPEC6, -CONF6) # remove some unnecessary columns

nrow(names) # check all species are present. Should equal 62

# check all 4 letter codes were paired with a name. Check there are no duplicate names
check <- names %>% 
  dplyr::select(SPEC, COMMONNAME, SCINAME) %>% 
  distinct() # keep only unique rows
print(check, n=Inf) # print all the rows

# there are a few names that need editing
# MAPS differentiates the subspecies of yellow-rumped wabler (Myrtle and Audubon's)
# we have models with AUWA or Audubon's warbler. The trait files do not have this level of specificity
# Similarly, MAPS differentiates subspecies of dark-eyed junco and we have models for Oregon dark-eyed junco
# this next step changes the scientific name for these two to reflect the broader species name, rather than the subspecies name
names$SCINAME[names$SCINAME=="Setophaga coronata auduboni"]<-"Setophaga coronata"
names$SCINAME[names$SCINAME=="Junco hyemalis oreganus"]<-"Junco hyemalis"

# combine with bird names and 4-letter codes from previous steps
DDspp_names <- names %>% 
  rename(Species1_BirdLife=SCINAME) %>% # change the name of the column that was called SCINAME to Species1_BirdLife
  left_join(., nameconvert, by="Species1_BirdLife") %>% # join names to nameconvert keeping all the rows in names and only rows from nameconvert that match
  dplyr::select(-Avibase.ID) %>% # drop column with Avibase.ID
  filter(Species3_BirdTree != "Dendroica aestiva") # remove aestiva group for yellow warbler (it shows up twice and we only want one record per species)

# confirm we still have 62 species/rows
nrow(DDspp_names)

# one is missing. It is BUOR
BUOR_name <- names %>%
  filter(SPEC == "BUOR") %>%
  rename(Species3_BirdTree = SCINAME) %>%
  left_join(., nameconvert) %>%
  dplyr::select(-Avibase.ID) # drop column with Avibase.ID

# join
DDspp_allnames <- bind_rows(DDspp_names, BUOR_name)

# confirm we have one-to-one match for scientific names between the 3 naming schemes
# these will print "TRUE" if we have one-to-one match
length(unique(DDspp_allnames$Species1_BirdLife)) == length(unique(DDspp_allnames$Species2_eBird)) 
length(unique(DDspp_allnames$Species1_BirdLife)) == length(unique(DDspp_allnames$Species3_BirdTree))

# Join the species density dependence measures and the bird scientific names 
DDspp_dat <- left_join(DDspp_names, ave_negative_slope_m2, by = "SPEC") %>%
  mutate(Species3_BirdTree = str_replace(Species3_BirdTree, " ", "_")) %>% # replace blanks in Genus species with underscore
  column_to_rownames(., var="Species3_BirdTree")

head(DDspp_dat)

# link data to phylogenetic tree
et <- treedata(tree_out, DDspp_dat, sort=T)

# Is there phylogenetic signal in the minimum adults to trigger density dependence (threshold)?
# using lambda as a measure of phylogenetic signal
lambda_threshold <- fitContinuous(et$phy, DDspp_dat[7], model = "lambda")
lambda_threshold 


## Is there phylogenetic signal in the average slope (intensity)?
lambda_slope <- fitContinuous(et$phy, DDspp_dat[11], model = "lambda")
lambda_slope

# Visualize estimates of threshold and intensity on tree

hist(DDspp_dat$min_adult)
hist(sqrt(DDspp_dat$min_adult))
hist(log1p(DDspp_dat$min_adult))

phytree <-et$phy # get phylo tree for plotting
circ <- ggtree::ggtree(phytree , layout='circular') # circular phylogeny

threshold <- data.frame(DDspp_dat[7]) # threshold values for each species
threshold_log1p <- data.frame(log1p(DDspp_dat[7])) # log x+1 transformed threshold values for each species

intensity <- data.frame(DDspp_dat[11])
intensity <- data.frame(log(abs(DDspp_dat[11])))


(threshold_plot <- gheatmap(circ, threshold_log1p, offset=.8, width=.2, colnames =F,
               colnames_angle=95, colnames_offset_y = .25) +
  scale_fill_viridis_c(option="A", name="lambda = 0.56\n \nMinimum\nAdult Abundance"))


(intensity_plot <- gheatmap(circ, intensity, offset=.8, width=.2, colnames =F,
                           colnames_angle=95, colnames_offset_y = .25) +
  scale_fill_viridis_c(option="A", name="lambda = 0.73\n \nMean\nSlope"))
