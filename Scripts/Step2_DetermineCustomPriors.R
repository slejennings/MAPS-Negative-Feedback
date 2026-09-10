###### MAPS Project: Negative Feedback #######
### Script name: Step2_DetermineCustomPriors.R
### Author(s): SLJ

########### Objective/Description of Script #####################
# identify potential custom priors to use for species-specific models
# we have many species (n=62) and need to use the same model settings across all of them
# this means we have to find priors that are weakly informative and work well across all models
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
library(gridExtra)
library(ggforce)
library(priorsense)
library(MetBrewer)
library(patchwork)
library(colorspace)
library(cmdstanr)
library(scales)

#####################################################################

# set theme for plots
theme_set(theme_default() +
            theme_minimal())

# set some colors for plotting
clrs <- MetBrewer::met.brewer("Juarez")

#####################################################################
# import data 
NF_breed_dat <- readRDS(here("Outputs", "NF_breed_dat.rds"))
summary(NF_breed_dat)
######################################################################
######### Examine data to decide on model structure #################


NF_spp_examine <- NF_breed_dat %>% 
  group_by(SPEC) %>% # group the data by species
  summarize(n = n(), # get the number of data points
            zeros = sum(FY_to_A==0), # get count of number of rows where productivity = 0
            min_adult= min(Adult), max_adult= max(Adult)) %>% # add min and max number of adults in data for each species (to be used later)
  mutate(percentzero = zeros/n*100) %>% # calculate the % of points that are zeros
  arrange(desc(n)) # arrange is descending order

# look at summary statistics for sample sizes for each species
hist(NF_spp_examine$n)
mean(NF_spp_examine$n) # 841
range(NF_spp_examine$n) # 353 up to 2689

# look at summary statistics for percent zeros for each species
hist(NF_spp_examine$percentzero)
mean(NF_spp_examine$percentzero) # 33%
range(NF_spp_examine$percentzero) # 12 up to 53

# make plots to show range of productivity values for each species
# color everything where productivity = 0 as red
# all other values of productivity in blue
prodplots <-
  NF_breed_dat %>%
  mutate(is_zero = FY_to_A == 0) %>% 
  mutate(FY_to_A = ifelse(is_zero, -0.1, FY_to_A))  %>%  # change all data points that equal zero to -0.1 to help with plotting
  ggplot(aes(x = FY_to_A)) +
  geom_histogram(aes(fill = is_zero), 
                 linewidth = 0.25, boundary = 0, color = "white") +
  geom_vline(xintercept = 0) + # draw a vertical line at zero. Bar to the left are productivity = 0 points
  scale_y_continuous(labels = label_comma()) +
  scale_fill_manual(values = c(clrs[2], clrs[4]), 
                    guide = guide_legend(reverse = TRUE)) +
  labs(x = "productivity", y = "count", fill = "is zero?") +
  theme(legend.position = "bottom") 

# print on multiple pages with one plot per species
prodplots + facet_wrap_paginate(~SPEC, scales="free", nrow=3, ncol=3, page = 1) # view page 1
prodplots + facet_wrap_paginate(~SPEC, scales="free", nrow=3, ncol=3, page = 2) # view page 2

# export plots in a single pdf
pdf(here("Figures", "productivityplots.pdf"), onefile = TRUE)
for(i in 1:7){ 
  print(prodplots + facet_wrap_paginate(~SPEC, scales="free", nrow =3 , ncol = 3, page = i))
}
dev.off()


# log transform all the data where productivity > 0 and plot this for each species
# if lognormal is a suitable distribution, the plots should look normally distributed
logplots <-
  NF_breed_dat %>%
  filter(FY_to_A > 0) %>% 
  mutate(log_prod = log(FY_to_A)) %>%
  ggplot(aes(x = log_prod)) +
  geom_histogram( linewidth = 0.25, boundary = 0, color = "white") +
  labs(x = "log transformed productivity", y = "count")

# print on multiple pages with one plot per species
logplots + facet_wrap_paginate(~SPEC, scales="free", nrow=3, ncol=3, page = 1) # view page 1
logplots + facet_wrap_paginate(~SPEC, scales="free", nrow=3, ncol=3, page = 2) # view page 2

# export plots in a single pdf
pdf(here("Figures", "logtransformplots.pdf"), onefile = TRUE)
for(i in 1:7){ 
  print(logplots + facet_wrap_paginate(~SPEC, scales="free", nrow =3 , ncol = 3, page = i))
}
dev.off()

# conclusion: we will use a hurdle lognormal model
# this uses a two part model: 
# 1) a logistic (hurdle) model for zeros and non zeros
# All the non-zero values pass the hurdle and go into: 
# 2) a lognormal model for all positive values

#####################################################################
# there are 3 possible models that may work well for these data
# all are hurdle lognormal models with different combinations of random effects for stations and years
# we will perform model comparison in a later step
# here, we are using the 2nd model, which is an intermediate level of complexity
# this is the model structure we think is likely best suited to the data

# this model contains random slopes and intercepts for stations, and random intercepts for years
# note: remember that we are ultimately running one model per species

# get the default priors for the model 
# (here we are pulling the default prior using all the data even though we will actually model each species separately)
default_prior(FY_to_A  ~ Adult + (1 + Adult|STA) + (1|year), # this is the lognormal model that models all the positive values
              hu ~ Adult + (1 + Adult|STA) + (1|year),
              data = NF_breed_dat, 
              family = hurdle_lognormal())

# we need priors for:
# the intercept and slope of the hurdle model (fixed effects)
# the intercept and slope of the lognormal model (fixed effects)
# the group-level standard deviation for the random intercepts of Year and Station in both the hurdle and lognormal models
# the correlation between random slopes and intercepts for (1 + Adult|STA) in both the hurdle and lognormal models
# and sigma (residual standard deviation) in the lognormal model

# for each prior, we will generate two options that are similar (both weakly informative)
# in the next script, we'll explore both options using prior sensitivity testing to see which works best across multiple species with varying sample sizes
# the goal is to find a set of weakly informative priors that work well across all species

######################## Hurdle/Logistic Model #############################

# priors for the logistic model are in log odds
# we need to convert back and forth between probability (easier to think about) and log odds (for specifying prior)
# use plogis to convert between log odds and probability
# use qlogis to convert between probability and log odds

###### Prior for the intercept of the hurdle model #######

# we know the data have many zeros, but we also know that positive values of productivity are, on average, more common than not
# the proportion of zeros for each species is likely between ~0.10 (10%) and  ~0.5 (50%), with a lot of variation among species
# this was exploratory analysis we already did to identify that the hurdle lognormal model was appropriate for our data

# examine histogram of productivity across the entire dataset (all species)
NF_breed_dat %>%
  mutate(is_zero = FY_to_A == 0) %>% 
  mutate(FY_to_A = ifelse(is_zero, -0.1, FY_to_A))  %>%  # change all data points that equal zero to -0.1 to help with plotting
  ggplot(aes(x = FY_to_A)) +
  geom_histogram(aes(fill = is_zero), 
                 linewidth = 0.25, boundary = 0, color = "white", binwidth=1) +
  geom_vline(xintercept = 0) + # draw a vertical line at zero. Bar to the left are productivity = 0 points
  scale_y_continuous(labels = label_comma()) +
  scale_fill_manual(values = c(clrs[2], clrs[4]), 
                    guide = guide_legend(reverse = TRUE)) +
  labs(x = "productivity", y = "count", fill = "is zero?") +
  theme(legend.position = "bottom") + xlim (-1,70)
# clearly, there are more positive values of productivity than zero values but zeros are common

# However, this covers the probability that productivity is equal to zero across all values of adult abundance
# the intercept in the hurdle model is the predicted probability that productivity is zero when adult abundance equals zero
# we have less knowledge about this value, so we want a prior that allows for a wide range of possibilities
# we'll try priors around 0.5 probability to reflect our lack of knowledge and keep wide standard deviation around the mean

# If we set the mean for the intercept prior at the probability of zero productivity = 0.5 
# then, we can use qlogis to convert between probability and log odds
qlogis(0.5) #  0

# we'll also explore the probability of zero productivity = 0.4
qlogis(0.4) # -0.4

# we also need a value for the standard deviation around the mean
# if we use 0 log odds as the mean, what does it mean to try 1 log odds as the std dev?
# convert log odds to probability using plogis
plogis(0 - (2*1)) # 0.12
plogis(0 + (2*1)) # 0.88
# this gives a range of 12% to 88% for the intercept around the mean of 50%
# this seems like a good, wide range for our weak prior

# if we use -0.4 log odds as the mean, what does it mean to try 1 log odds as the std dev?
# convert log odds to probability using plogis
plogis(-0.4 - (2*1)) # 0.08
plogis(-0.4 + (2*1)) # 0.83
# this gives a range of 8% to 83% for the intercept around the mean of 40%

# we will use student_t with 3 degrees of freedom to test these two options to leave more area in the tails (for extreme values) compared to a normal distribution

# plot the default and the two custom priors for hurdle model intercept to visualize how they differ
hu_prior_intercept <- tibble(
  x = seq(from = -10, to = 10, by = 0.1),
  default = dlogis(x, 0, 1), # default from brms
  custom1 = dstudent_t(x, 3, 0, 1), # custom prior option 1
  custom2 = dstudent_t(x, 3, -0.4, 1)) %>% # custom prior option 2
  pivot_longer(., cols = -x, names_to = "prior")

ggplot(hu_prior_intercept,
  aes(x = x, y = value, linetype = fct_rev(prior), color = fct_rev(prior))) +
  geom_line() +
  labs(x = "log odds", y = "", linetype = "Prior", color = "Prior") +
  scale_x_continuous(breaks = seq(-10, 10, 5)) +
  labs(title="Hurdle model intercept") +
  theme_classic()

###### Prior for the slope of the hurdle model #######

# the slope is linear on the log odds scale
# for every increase in adult abundance of 1 adult, what will happen to the log odds?
# it can help to think about this on the nonlinear odds scale to make sense of it

# to calculate odds from probability:
# odds = probability/(1-probability)
# if probability = 0.5
0.5/(1-0.5) # odds equal 1
# odds range from 0 up to infinity (whereas probability is restricted to 0 to 1)
# when odds = 1,there is a 50-50 chance of an event
# if odds > 1, then an event's chances are greater than 50-50 (i.e. probability of the event is > 50%)
# if odds < 1, then an event's chances are less than 50-50 (i.e. probability of the event is < 50%)

# to calculate probability from odds:
# probability = odds/(1+odds)

# finally, to convert log odds to odds, use exp(logodds)
# and log odds from odds, use log(odds)

# For our models, we predict that the slope will be negative
# for each increase of one adult, the probability of obtaining productivity = 0 should decrease 
# however, we do not have prior knowledge of how steep the slope is
# so we want a prior for the slope that skews negative but allows for a large range of possible values

# if the mean of the prior for the slope = -0.1, what does this mean in odds?
exp(-0.1) # 0.90 <- the outcome of productivity = 0 is ~ 10 % less likely due to an increase of one adult
# if we use a standard deviation of 1 around -0.1 ?
exp(-0.1 - (2*1)) # 0.12
exp(-0.1 + (2*1)) # 6.69
# mostly negative slopes with a wide range and some possibility of positive responses

# if the mean of the prior for the slope = -0.2, what does this mean in odds?
exp(-0.2) # 0.82 <- the outcome of prod = 0 is ~ 18 % less likely due to an increase of one adult
# this is a slightly steeper slope than above
# if we use a standard deviation of 1 around -0.2 ?
exp(-0.2 - (2*1)) # 0.11
exp(-0.2 + (2*1)) # 6.049
# mostly negative slopes with a wide range and some possibility of positive responses

# we will use student t distribution with 3 degrees of freedom to leave more area in the tails compared to a normal distribution

# plot the default and custom priors for hurdle model slope to visualize how they differ
# note: the default slope is flat
hu_prior_slope <- tibble(
  x = seq(from = -3, to = 3, by = 0.1),
  default = dunif(x, min=-3, max=3), # the default prior is flat or uniform
  custom1 = dstudent_t(x, 3, -0.1, 1), # custom prior option 1
  custom2 = dstudent_t(x, 3, -0.2, 1)) %>% # custom prior option 2
  pivot_longer(. , cols = -x, names_to = "prior")

ggplot(hu_prior_slope,
  aes(x = x, y = value, linetype = fct_rev(prior), , color = fct_rev(prior))) +
  geom_line() +
  labs(x = "log odds", y = "", linetype = "Prior", color="Prior") +
  scale_x_continuous(breaks = seq(-3, 3, 1)) + 
  labs(title="Hurdle model slope") + 
  theme_classic()

######################## Lognormal Model #############################

# the priors for the lognormal model are on the log scale
# to get back to regular scale, use exp()
# to convert "normal" values to log scale, use log()

###### Prior for the intercept of the lognormal model #######

# the lognormal model only has positive values of productivity

# if we return to the same histogram of productivity across the entire dataset that we examined before:
NF_breed_dat %>%
  mutate(is_zero = FY_to_A == 0) %>% 
  mutate(FY_to_A = ifelse(is_zero, -0.1, FY_to_A))  %>%  # change all data points that equal zero to -0.1 to help with plotting
  ggplot(aes(x = FY_to_A)) +
  geom_histogram(aes(fill = is_zero), 
                 linewidth = 0.25, boundary = 0, color = "white", binwidth=0.5) +
  geom_vline(xintercept = 0) + # draw a vertical line at zero. Bar to the left are productivity = 0 points
  scale_y_continuous(labels = label_comma()) +
  scale_fill_manual(values = c(clrs[2], clrs[4]), 
                    guide = guide_legend(reverse = TRUE)) +
  labs(x = "productivity", y = "count", fill = "is zero?") +
  theme(legend.position = "bottom") + xlim (-1,40)
# we can see that the most frequent positive value for productivity is 1 or less 
# the first two blue columns in this version of the histogram shows the count for productivity values > 0 and <= 1 (binwdith = 0.5)
# but we can also see that there are some values >1 that pull the average of the positive productivity values above 1


# the intercept in the lognormal model reflects the predicted value of productivity when adult abundance is zero for observations where productivity is > 0
# In our data, we have no positive values of productivity where adult abundance = 0 because we assumed that any hatch year bird must be associated with at least 1 adult, and thus made adult abundance equal to 1 in these situations
# therefore, the value for the y-intercept in the lognormal model for each species is largely being determined data points with low values for adult abundance and productivity > 0

# we also know that:
# the effects of density dependence is likely to be lower when the population of breeding adults is small (low adult abundance)
# many songbirds have clutches of > 1 egg and possibly multiple clutches in a year, but many offspring will fail to hatch/fledge
# moreover, there are likely adults caught in the mist nets who do not breed; these individuals may balance out hatch year birds
# therefore, it seems reasonable that the y-intercept is around productivity = 1 with a lot of wiggle around it

# if the mean of the intercept prior is productivity = 0.9. What is this on the log scale?
log(0.9) #  ~ 0.1

# if the mean of the intercept prior is productivity = 0.80 (slightly less than 1). What is this on the log scale?
log(0.8) #  ~ -0.2
# this is what brms picks as the default prior mean

# if the mean of -0.1 has a standard deviation of 1 on the log scale. What does this mean for productivity?
exp(-0.1 + (2*1)) # 6.68
exp(-0.1 - (2*1)) # 0.12

 # if the mean of -0.2 has a standard deviation of 1 on the log scale. What does this mean for productivity?
exp(-0.2 + (2*1)) # 6.05
exp(-0.2 - (2*1)) # 0.11

# plot the default and custom priors for lognormal model intercept to visualize how they differ
lognormal_prior_intercept <- tibble(
  x = seq(from = -10, to = 10, by = 0.1),
  default = dstudent_t(x, 3, -0.2, 2.5), # default from brms
  custom1 = dnorm(x, -0.1, 1), # custom prior option 1
  custom2 = dnorm(x,-0.2, 1)) %>% # custom prior option 2
  pivot_longer(., cols = -x, names_to = "prior")

ggplot(lognormal_prior_intercept,
       aes(x = x, y = value, linetype = fct_rev(prior), color = fct_rev(prior))) +
  geom_line() +
  labs(x = "log", y = "", linetype = "Prior", color = "Prior") +
  scale_x_continuous(breaks = seq(-10, 10, 5)) +
  labs(title="Lognormal model intercept") +
  theme_classic()


###### Prior for the slope of the lognormal model #######

# if negative density dependence is occurring, then slope of the lognormal model should be negative
# but we have no prior knowledge of how steep the slope is, so we need a weak prior that allows for a wide range of values

# if we say slope = -0.1 on the log scale
# what does that mean when converted?
1-exp(-0.1) # this is the easiest way to interpret
# each additional adult results in a 9.5% decline in productivity 
# NOTE: in reality the line isn't straight when not on the log scale so this approach is only so helpful

# if we say slope = -0.2 on the log scale. This is a slightly steeper slope
1-exp(-0.2) # this is the easiest way to interpret
# each additional adult results in an 18.1% decline in productivity 
# NOTE: in reality the line isn't straight when not on the log scale so this approach is only so helpful

# with standard deviation of 1 around the mean
(-0.1 + (2*1)) # 1.9
(-0.1 - (2*1)) # -2.1

1 - exp(1.9) 
1 - exp(-2.1) 

# this second prior option is skewed slightly more towards negative slopes than the option above
(-0.2 + (2*1)) # 1.8
(-0.2 - (2*1)) # -2.2
1 - exp(1.8) 
1 - exp(-2.2)

#  we will use a normal distribution for this prior as this already seems fairly weak

# plot the default and custom priors for lognormal model slope to visualize how they differ
lognormal_prior_slope <- tibble(
  x = seq(from = -4, to = 4, by = 0.01),
  default = dunif(x, min =-4, max=4), # the default prior is flat or uniform
  custom1 = dnorm(x, -0.1, 1), # custom prior option 1
  custom2 = dnorm(x, -0.2, 1)) %>% # custom prior option 2
  pivot_longer(. , cols = -x, names_to = "prior")

ggplot(lognormal_prior_slope,
       aes(x = x, y = value, linetype = fct_rev(prior), , color = fct_rev(prior))) +
  geom_line() +
  labs(x = "log", y = "", linetype = "Prior", color="Prior") +
  scale_x_continuous(breaks = seq(-4, 4, 1)) + 
  labs(title="Lognormal model slope") +
  theme_classic()


######################################################################################
# there are 3 other classes of priors to consider: 
# 1) priors for sd around random intercepts
# 2) sigma in the lognormal model
# 3) correlation between random slopes and intercepts for Stations
# we will leave all these priors as using the brms defaults

# plot the default and custom priors for SD of random intercepts
# standard deviation values must be positive, so brms implements a half student t distribution
sd_prior <- tibble(
  x = seq(from = 0, to = 10, by = 0.01),
  default = dstudent_t(x, 3, 0, 2.5)) %>% # the default prior
  pivot_longer(. , cols = -x, names_to = "prior")

ggplot(sd_prior,
       aes(x = x, y = value, linetype = fct_rev(prior), , color = fct_rev(prior))) +
  geom_line() +
  labs(x = "", y = "", linetype = "Prior", color="Prior") +
  scale_x_continuous(breaks = seq(0, 10, 1)) + 
  labs(title="SD for random intercepts of STA and Year") +
  theme_classic()

# plot the default and custom priors for sigma in the lognormal model
# sigma must be positive, so brms implements a half student t distribution
sigma_prior <- tibble(
  x = seq(from = 0, to = 10, by = 0.01),
  default = dstudent_t(x, 3, 0, 2.5)) %>% # the default prior
  pivot_longer(. , cols = -x, names_to = "prior")

ggplot(sigma_prior,
       aes(x = x, y = value, linetype = fct_rev(prior), , color = fct_rev(prior))) +
  geom_line() +
  labs(x = "", y = "", linetype = "Prior", color="Prior") +
  scale_x_continuous(breaks = seq(0, 10, 1)) + 
  labs(title="Sigma") + 
  theme_classic()


###############################################################################
# define two groups of custom priors

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
