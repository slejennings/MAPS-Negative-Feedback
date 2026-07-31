###### MAPS Project: Density Dependence #######
### Script name: Step1_DataPrepforDensityDependenceModels.R
### Author(s): SLJ

########### Objective/Description of Script #####################
# prepare and filter MAPS banding/capture data for use in models
# data used here comes from 1989 through 2018
#################################################################

####### Step 1: Prep MAPS Data for Density Dependence Models ########

### Setup ####

# load packages
library(tidyverse)
library(here)


# Import bird capture data 
capturedat <- read.csv(here("Data", "MAPS_BANDING_capture_data.csv"))
head(capturedat)
summary(capturedat)

################################################

# Examine species that are in the data
unique(sort(capturedat$SPEC))

# look at PSFL, COFL, and WEFL, which should all be classified under WEFL
flycatchers <-capturedat %>% filter(SPEC %in% c("WEFL", "PSFL", "COFL"))
nrow(flycatchers) # there are 16219 rows
flycatchers %>% count(SPEC) # most were identified as PSFL

# make the change to label all as WEFL
dat <- capturedat %>% mutate(SPEC = case_match(SPEC, c("PSFL", "COFL") ~ "WEFL", .default = SPEC))

# do we still see PSFL and COFL in the list of species names?
sort(unique(dat$SPEC)) # no

# confirm that all flycatchers are now listed as WEFL
dat %>% filter(SPEC=="WEFL") %>% nrow() == nrow(flycatchers) # should equal TRUE

################################################

# Examine data for recaptured birds and other codes that indicate birds should not be incorporated into productivity or adult abundance analyses

# examine entries using column "C" which denotes recaptured individuals
recap <- dat %>% filter(C=="R") %>% arrange(BAND) # R is the code for recaptured birds

# look at values for BAND, which denotes the band number
unique(recap$BAND) # some recaptured birds have no band number (coded as "")
recap_noband <- recap %>% filter(BAND=="") # keep only entries that are missing band numbers
# which recaptured species are missing band numbers?
unique(recap_noband$SPEC)
# multiple species of hummingbirds, Northern bobwhite (NOBO), Bewick's wren, Wrentit

# Next, examine data using BAND (band number) as an alternative to using column C
band <- dat %>% 
  filter(BAND!="") %>% # remove entries that are missing band numbers
  group_by(BAND) %>% # group the data by band number
  filter(n()>1) %>% # retain all band numbers that occur more than once (birds that are recaptured)
  arrange(BAND, year) # sort data frame by band number and year

# By examining "band", we can see that recaptured birds can occur in several ways
# First way: an individual is recaptured within a breeding season.
# We need to remove instances where individual birds were recaptured within a particular year/breeding season
# Second way: recaptured across multiple breeding seasons
# As we are looking at changes in the productivity and abundance over time, we should count each individual once for each year it was caught
# So we do not want to remove all instances where column C == R

# Make 3 separate data frames to keep track of different types of captures
# data frame 1: entries with no band number
no_band <- dat %>%
  filter(BAND=="")

# data frame 2: banded individuals that were not recaptured (band number only occurs once)
one_band <- dat %>%
  filter(BAND!="") %>% # remove entries that are missing band numbers
  group_by(BAND) %>% # group the data by band number
  filter(n()==1) # keep all band numbers that occur only once
# confirm there are no band numbers with duplicate entries:
sum(duplicated(one_band[,5])) # BAND is the 5th column. This should equal 0 if there are no duplicates

# data frame 3: banded individuals that were recaptured, but modified to only include one entry for each year an individual was sighted/caught
band_distinct <- dat %>%
  filter(BAND!="") %>% # remove entries that are missing band numbers
  group_by(BAND) %>%
  filter(n()>1) %>% # retain all band numbers that occur more than once (birds that are recaptured)
  distinct(year, .keep_all = T) %>% # keep only one entry per year
  arrange(BAND)

# combine 3 types back into a single data frame
combine <- bind_rows(no_band, one_band, band_distinct)
nrow(dat) - nrow(combine) # the above steps removed 288466 rows

# Finally, use column "N" for a final cleaning step
# this column has the value of "-" for records that can be used in productivity and abundance analyses
unique(combine$N) # "G" for gallinaceous bird, "H" for hummingbird, "U" for unbanded should be removed
dat2 <- combine %>% filter(N=="-") # keep only records coded as "-"
head(dat2)
nrow(combine)-nrow(dat2) # removed 103347 rows

######################################################
#### Classifying age of captured birds and removing certain stations and/or years from data

# Classify birds caught based on their age
dat_age <- dat2 %>%
  filter(AGE != "0") %>% # remove birds with unknown age
  mutate(AgeCat = if_else(AGE== 4 | AGE == 2, "FirstYear", "Adult")) # classify hatch year birds as "FirstYear", all others "Adult"

# Find number of years of operation for each station and reduce the list of stations to those that have 10 or more years
NumYrsOp <- dat_age %>%
  group_by(STA) %>% # group by STA which is a unique identifier for station
  summarize(NumYrsOp = n_distinct(year)) %>% # find the number of years for each station
  filter(NumYrsOp >=10) # reduce the list to stations with 10 or more years of data

head(NumYrsOp)
nrow(NumYrsOp) # 394 stations with 10 or more yrs of data
print(NumYrsOp, n=Inf)

# Identify stations with more than 10 years of data that are non-consecutive in time
NonConsecutive <- dat_age %>%
  group_by(STA) %>% # group by STA which is a unique identifier for station
  arrange(year) %>%
  distinct(year) %>%
  summarize(NumYrsOp = n_distinct(year), # find the number of years for each station
            MinYr = min(year), # the earliest year
            MaxYr = max(year), # the most recent year
            YrList = str_c(year, collapse=" ")) %>%  # print list of the years with data
  filter(NumYrsOp >= 10) %>% # restrict to stations with 10 or more years
  mutate(TimeRange = ((MaxYr-MinYr)+1), # find number of years between earliest (min) and most recent (max)
         Consecutive = ifelse(TimeRange==NumYrsOp, "yes", "no")) %>% # identify stations with non consecutive data 
  filter(Consecutive=="no") %>% # filter to keep only stations with non consecutive data
  mutate(YrsDiff = TimeRange - NumYrsOp) %>% arrange(desc(YrsDiff)) # find the difference between the range of years and the number of years of operation 

print(NonConsecutive, n=20)
nrow(NonConsecutive) # 82 stations out of 394 with broken time series

# How much of a lag (in years) exists for each station?
Lag <- dat_age %>%
  group_by(STA) %>% # group by STA which is a unique identifier for station
  arrange(year) %>%
  distinct(year) %>%
  mutate(yrlag = year - lag(year)) %>% # calculates the lag between the consecutive rows of year
  arrange(STA, year) %>%
  filter(yrlag > 5) # look at everywhere that the lag is more than 5 years
# Note: this is all the stations, not restricted to 10 or more years of data

# combine lag and non consecutive lists to get a list of which stations with 10 or more years of data may be problematic
STA_10yr_Lag <- inner_join(NonConsecutive, Lag) %>%
  select(STA, NumYrsOp, year, yrlag, YrList) %>%
  arrange(yrlag)

print(STA_10yr_Lag, n=Inf)
# 12 rows but 11 unique stations
# We decided what to do with each of these station by manually examining it's time series
# The decisions were as follows:
# 15591 - Drop. Has two 6 year gaps
# 15650 - Drop. Doesn't have 10 yrs after removing the years before gap
# 16658 and 16659 - Keep but remove 2015, 2016
# 16696 - Keep 2007 - 2016 but drop 2001
# 13319 - Drop. Seven year gap in middle of time series
# 11109 and 11110 - Keep. Remove 1992 and 1993
# 15595 - Drop. Gap in middle of time series
# 16609 - Drop. Gap in middle of time series
# 15510 - Keep. Remove 2015-2018 and keep 1992-2001


# Make a data frame of the stations that are being kept but have certain years that need to be removed
# these will all still have 10 or more years of data after dropping these years
remove <- data.frame(STA = as.integer(c(rep("16658", 2), rep("16659", 2),"16696", rep("11109",2), rep("11110",2), rep("15510",4))),
                     year = as.integer(c(rep(c("2015", "2016"), 2), "2001", rep(c("1992","1993"),2), "2015", "2016", "2017", "2018")))
head(remove)

# make final list of stations and the number of years they operated after making the changes for lags
STAYrsOp <- dat_age %>%
  filter(!STA %in% c("15591", "15650", "13319", "15595", "16609")) %>% # remove stations to be dropped
  anti_join(., remove) %>% # drop certain years for stations with large time lags
  group_by(STA) %>% # group by STA which is a unique identifier for station
  summarize(NumYrsOp = n_distinct(year)) %>% # find the number of years for each station
  filter(NumYrsOp >=10) # reduce the list to stations with 10 or more years of data

nrow(STAYrsOp) # 389 stations  
print(STAYrsOp, n=20)  
  
# Reduce data to only include stations in the final list
# Remove problem years for stations with large effort lags
dat_10yr <- dat_age %>%
  inner_join(., STAYrsOp) %>% # join with STAYrsOp to reduce dat_age to so that it only includes the list of stations in STAYrsOp
  anti_join(., remove) # remove years from stations with long time lags

nrow(dat_age) - nrow(dat_10yr) # removed 511866 rows
head(dat_10yr)

saveRDS(dat_10yr, here("Outputs", "dat_10yr.rds")) # save this into outputs as we will use it later to extract morphometric measurements

# Find counts of adults and hatch years for each station, species, and year
allspp <- dat_10yr %>% 
  group_nest(STA, SPEC, AgeCat) %>% # put each station, species and age category into its own tibble
  mutate(
    yrcount = map(data,. %>% group_by(year) %>% summarize(NumObs=n())) # within each tibble, group by the year and find the number observed
  ) %>%
  unnest(yrcount) %>% # unnest to see results in a flat tibble
  arrange(STA, SPEC)

print(allspp, n=100) # check this worked as expected by looking at the first 100 rows

# Reorganize the data frame slightly 
allspp_age <- allspp %>% 
  dplyr::select(-data) %>% # remove the data
  pivot_wider(., names_from=AgeCat, values_from = NumObs) # make wider to put Adult counts and FirstYear counts into their own columns

head(allspp_age)

###############################################################################
# Investigate the Adult and Hatch year counts before proceeding

# examine distributions and values for Adult:
unique(allspp_age$Adult) # does not include zero. But does include NA
sum(is.na(allspp_age$Adult)) # NA occurs in multiple places in the column
sort(unique(allspp_age$Adult)) # look at the same thing as above but in ascending order to more easily examine the range of values. This will drop NA though
hist(allspp_age$Adult) # look at a histogram
boxplot(allspp_age$Adult)# look at a boxplot

# examine distributions and values for FirstYear (hatch year birds):
unique(allspp_age$FirstYear) # does not include zero. But does include NA
sum(is.na(allspp_age$FirstYear)) # NA occurs in multiple places in the column
sort(unique(allspp_age$FirstYear))# look at the same thing as above but in ascending order to more easily examine the range of values. This will drop NA though
hist(allspp_age$FirstYear) # look at a histogram
boxplot(allspp_age$FirstYear) # look at a histogram

# Another check: confirm that within each station there is a row for each species in each of the years the station was in operation
dat.check <- allspp_age %>% 
  group_by(STA, SPEC) %>% # group by the station and species
  summarize(SppYrs = n()) %>% # count the number of years associated with each species
  inner_join(., STAYrsOp, by="STA") %>% # add the number of years of operation for each station
  mutate(check = if_else(SppYrs==NumYrsOp, "Ok", # compare number of years for a species with number of years at that station
                         if_else(SppYrs < NumYrsOp, "Less", "More"))) # identify rows where they are not equal. Mark if SppYrs is either Less than or More than NumOpYrs

unique(dat.check$check) # there are some instances where a species does not have a row for every year (denoted by "Less")
# how many?
dat.check %>% group_by(check) %>% summarize(count=n()) # a lot!
# A species only shows up at a station in a particular year if either 1 or more hatch year birds or 1 or more adults were caught
# if none of either was caught, then the data does not have an entry for that species

# There are a couple of things need to happen at this stage to resolve the Adult and FirstYear counts:
# step 1: we need to convert the rows where the counts for Adult or FirstYear are currently NA to zero
# step 2: we also need to add rows at each station for the years where the species was not observed/caught

# First, change all the NAs to zero:
allspp_age <- allspp_age %>%
  mutate(Adult=replace_na(Adult,0), # replace all NAs in column Adult with 0
         FirstYear = replace_na(FirstYear, 0)) # replace all NAs in column FirstYear with 0

# check there are no NAs left in these two cxfolumns:         
sum(is.na(allspp_age$Adult)) # should equal 0
sum(is.na(allspp_age$FirstYear)) # should equal 0

# Next, expand the data to get every combination of years/species at each station
allspp_expand <- allspp_age %>% 
  group_by(STA) %>% # group by station ID
  tidyr::expand(SPEC, year) # expand the data to make sure all combinations of species and years are present

# did this work as expected?
allspp_expand %>% group_by(STA) %>% summarize(count=n_distinct(SPEC)) # do we have different number of species at each station? yes
# next, confirm that SppYrs==NumYrsOp in allspp_expand
equal.test <- allspp_expand %>%
  group_by(STA, SPEC) %>% # group by the station and species
  summarize(SppYrs = n()) %>% # count the number of years associated with each species
  inner_join(., STAYrsOp, by="STA") %>% # add the number of years of operation for each station
  mutate(check = if_else(SppYrs==NumYrsOp, "Ok", # compare number of years for a species with number of years at that station
                         if_else(SppYrs < NumYrsOp, "Less", "More"))) # identify rows where they are not equal. Mark if SppYrs is either Less than or More than NumOpYrs

unique(equal.test$check) # looks good! All should be labeled "Ok"

# Finally, join adult and hatch year counts with new expanded data and make all rows with missing values have a count of zero
allspp_combine <- allspp_expand %>%
  left_join(., allspp_age) %>% # left_join will keep all rows in allspp_expand and add rows that match in allspp_age
  mutate(Adult=replace_na(Adult,0), # replace all NAs in column Adult with 0
         FirstYear = replace_na(FirstYear, 0)) # replace all NAs in column FirstYear with 0
print(allspp_combine, n=50)

###############################################################################
#### Calculate Productivity ####

# Calculate productivity as the ratio of hatch year (FirstYear) birds to adults
# Apply several rules when obtaining this estimate:
# 1. if Adult = 0 and FirstYear = 0, then productivity = 0
# 2. if Adult > 0 and FirstYear = 0, then productivity = 0
# 3. if Adult = 0 and FirstYear > 0, then productivity = FirstYear/(Adult+1)
# 4. if Adult > 0 and FirstYear > 0, then productivity = FirstYear/Adult
# the column HY_to_A contains productivity estimate
allspp_prod <- allspp_combine %>%
  mutate(FY_to_A = if_else(FirstYear==0 & Adult ==0, 0, # if both adults and hatch years were not observed, make productivity equal to zero
                           if_else(FirstYear == 0 & Adult > 0, 0, # if some adults but no hatch year birds were observed, make productivity equal to zero
                                   if_else(FirstYear > 0 & Adult == 0, ((FirstYear)/(Adult+1)), # if no adults were observed but some hatch years were detected, add one adult to denominator 
                                           (FirstYear/Adult))))) # in all other situations, calculate productivity as the ratio of hatch year to adult birds

###############################################################################
#### Prep Data for Density Dependence Models ####

# In this section we apply some restrictions to remove stations/species combinations where either Adults or Productivity are heavily zero-inflated
# For density dependence, we are not making a model for each station (rather for each species)
# However, we want to use station as a random intercept in these models, so it is still useful to only keep stations where data isn't mostly zeros

# Filter data using the adult abundance values

# Steps:
# 1. Find the total number of adults of each species caught at a station over all the years it has operated. 
# 2. Count number of data points within each species/station combination that are equal to zero
# 3. Filter to retain only species/station combinations where number of adults is >= number of years of operation
# 4. Filter to retain only species/station combinations where 1/3 or more of Adult data points are >0 (filter out highly zero-inflated models)

Adult_summary <- allspp_prod %>% 
  group_by(STA, SPEC) %>% # 
  summarize(TotalAdult=sum(Adult), # find total for Adults of each species at each station across all years the station operated
            NumAdultZeros = sum(Adult==0)) %>% # count number of rows/data points that are equal to zero
  inner_join(.,STAYrsOp, by="STA") %>% # combine with number yrs of operation for station
  mutate(ZeroAdultCheck= NumYrsOp-NumAdultZeros) %>% # make column to identify (stations/species) where they have only a few data points where Adults > 0
  filter(TotalAdult >= NumYrsOp) %>% # keep only rows (stations/species) where number of adults is greater than or equal to number of years of operation
  filter(ZeroAdultCheck >= (round((NumYrsOp/3),0))) # keep only rows (stations/species) where 1/3 of the data points are non-zeros

# Combine summary stats from Adult_summary with data in allspp_prod
allspp_Adult <- allspp_prod %>%
  inner_join(., Adult_summary, by=c("STA", "SPEC")) 

print(allspp_Adult, n=30)

# Next, apply similar criteria for hatch year birds

# Steps:
# 1. Find the total number of hatch year (FirstYear) birds  of each species caught at a station over all the years it has operated
# 2. Count number of data points within each species/station that are equal to zero
# 3. Filter to retain only species/station combinations where 1/3 or more of hatch year data points for the model are >0 (filter out highly zero-inflated models)

FY_summary <- allspp_Adult %>% 
  group_by(STA, SPEC) %>% # 
  summarize(TotalFY = sum(FirstYear),  # find total hatch year birds of each species at each station across all years the station operated
            NumFYZeros = sum(FirstYear==0)) %>% # count number of rows/data points that are equal to zero
  inner_join(.,NumYrsOp, by="STA") %>% # combine with number years operation for stations
  mutate(ZeroFYCheck= NumYrsOp-NumFYZeros) %>% # make column to identify (stations/species) where they have only a few data points where FirstYear > 0
  filter(ZeroFYCheck >= (round((NumYrsOp/3),0))) # keep only rows (stations/species) where 1/3 or more of the data points are non-zeros


# Combine summary stats from FY_summary with data in allspp_Adult
allspp_AdultProd <- allspp_Adult %>%
  inner_join(., FY_summary, by=c("STA", "SPEC", "NumYrsOp")) %>%
  select(STA, SPEC, year, Adult, FirstYear, FY_to_A) # select columns to keep in final data frame

print(allspp_AdultProd, n=30)
# STA is the station ID
# SPEC is the species
# year is the year
# Adult is the number of adults (aka adult abundance)
# FirstYear is the number of hatch year birds (abudance of hatch years)
# FY_to_A is the productivity value, which was calculated as the ratio of hatch (or first) years birds to adults


# examine the structure of this object
str(allspp_AdultProd)
# for the models, we will want SPEC, STA, and year to be stored as factors, which they are currently not

# change certain columns to be stored as factors
allspp_AdultProd <- allspp_AdultProd %>%
  mutate_at(c("STA", "SPEC", "year" ), as.factor) # change multiple columns to factors

# examine sample size for each species
DD_spp <- allspp_AdultProd %>%
  group_by(SPEC) %>%
  count() %>%
  arrange(desc(n))

# Above, we made the decision that if # of hatch year birds in a given year for a species/station combination was >0, then we would calculate productivity using Adult =1
# to be consistent with that decision, we also need to change the raw adult abundance counts (whereas so far we have only modified this number inside the calculation for productivity)
# Look at the number of times we see adult = 0, prod > 1
adult0 <- allspp_AdultProd %>%
  mutate(class = if_else(Adult == 0 & FY_to_A > 0, "Yes", "No")) %>%
  group_by(SPEC) %>%
  count(class) %>%
  pivot_wider(names_from = class, values_from=n) %>%
  mutate(percentage = round((Yes/No)*100, 1)) %>%
  select(-Yes, -No) %>%
  inner_join(., DD_spp) %>%
  arrange(desc(n))

print(adult0, n=Inf)

# fix this. Make Adult = 1 if FY_to_A is > 0
DD_dat <- allspp_AdultProd %>%
  mutate(Adult = if_else(Adult == 0 & FY_to_A > 0, 1, Adult))

# For the final filtering steps, check that all species-station observations being used are from birds classified as breeding at that location
# And reduce to species that have at least 350 data points

# breeding status codes are contained in BRTSTAT column. Look at the unique values
spp_status <- dat_10yr %>% 
  select(STA, SPEC, year, BRSTAT) %>%
  distinct()

brstat <- spp_status %>% 
  group_by(STA, SPEC, BRSTAT) %>% 
  count() %>%
  arrange(STA, SPEC)

# find station/species combinations that were given more than one BRSTAT across the years of operation
multistatus <- brstat %>% 
  group_by(STA, SPEC) %>%
  filter(n() >1) %>%
  arrange(STA, SPEC)
print(multistatus) # all WEFL. Probably coming from splitting and then binning this species earlier in this script

# I examined the classifications for WEFL (and associated species) at these stations
# conclusion - all stations had years where WEFL was classified as some kind of breeder, so changing this to be the status across all years
WEFL_fix <- multistatus %>%
  pivot_wider(., names_from=BRSTAT, values_from=n) %>%
  mutate(BRSTAT = "B") %>%
  select(STA, SPEC, BRSTAT)

# combine fixed WEFL data with other species-stations
brstat_clean <- anti_join(brstat, multistatus) %>%
  select(-n) %>%
  bind_rows(WEFL_fix) %>%
  mutate(STA = factor(STA))

# combine with DD_dat
DD_brstat <- DD_dat %>% 
  left_join(., brstat_clean)

# identify if any are missing breeding status
missing <- DD_brstat %>%
  filter(is.na(BRSTAT))
nrow(missing) # all rows have a breeding status

# look at list of species that have 350 or more data points
DD_spp_examine <- DD_brstat %>% 
  group_by(SPEC) %>% # group the data by species
  summarize(n = n(), # get the number of data points
            zeros = sum(FY_to_A==0), # get count of number of rows where productivity = 0
            min_adult= min(Adult), max_adult= max(Adult)) %>% # add min and max number of adults in data for each species (to be used later)
  mutate(percentzero = zeros/n*100) %>% # calculate the % of points that are zeros
  filter(n > 350) %>% # restrict to species with 350 data points or more
  arrange(desc(n)) # arrange is descending order

SPEC350 <- DD_spp_examine %>% pull(SPEC) # get a vector of species with more than 350 points

# filter DD_brstat to retain only these species
DD_brstat_SPEC <- DD_brstat %>%
  filter(SPEC %in% SPEC350)

# make a vector of non-breeding classification codes
nonbreeder <- c("T", "A", "M")

# pull all the data points for non-breeding species/station combinations
nonbreeding <- DD_brstat_SPEC %>%
  filter(BRSTAT %in% nonbreeder) 
nrow(nonbreeding) # 1078

# get clean dataset with final species (n = 62 species) and only breeding birds
DD_breed_dat <- DD_brstat %>%
  filter(SPEC %in% SPEC350) %>% # reduce to species with > 350 points
  anti_join(., nonbreeding) # remove records from stations where species are not breeding

length(unique(DD_breed_dat$SPEC)) # 62 species

# export data
saveRDS(DD_breed_dat, here("Outputs", "DD_breed_dat.rds"))

