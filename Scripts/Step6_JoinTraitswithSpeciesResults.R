###### MAPS Project: Density Dependence #######
### Script name: Step6_JoinTraitswithSpeciesResults.R
### Author(s): SLJ, JML

########### Objective/Description of Script #####################
# join together the species specific measures of density dependence with species specific functional traits
# we will examine traits that fall into 4 groups: morphometric traits, trophic/diet traits, life history traits, and sexual selection traits
# species values for these traits come from a variety of previously published databases. References for each can be found in the readme file
#################################################################

### Setup ###

# load packages
library(tidyverse)
library(here)
library(ape)
library(geiger)
library(confintr)

#################### Combine Bird Codes with Bird Scientific Names ##################


# import data 
DD_breed_dat <- readRDS(here("Outputs", "DD_breed_dat.rds"))

# find the number of data points for each species
DD_SPEC <- DD_breed_dat %>%
  group_by(SPEC) %>% # group the data by species
  summarize(n = n()) %>% # get the number of data points for each species
  arrange(desc(n)) # arrange is descending order
  

# import bird list used by MAPS to convert 4-letter bird codes to scientific names
birdcodes <- read.csv(here("Data", "IBP-AOS-LIST23.csv"), header=T) 
head(birdcodes)  
# downloaded from IBP website: https://www.birdpop.org/pages/birdSpeciesCodes.php

# combine names with species codes
names <- DD_SPEC %>% 
  left_join(., birdcodes, by="SPEC") %>% # left join keeps all the rows in DD_breed_dat and adds anything that matches from birdcodes
  dplyr::select(-SP, -CONF, -SPEC6, -CONF6) # remove some unnecessary columns

head(names)
nrow(names) # check all species are still present. Should equal 62

# check all 4 letter codes were paired with a name. Check there are no duplicate names
check <- names %>% 
  dplyr::select(SPEC, COMMONNAME, SCINAME) %>% 
  distinct() # keep only unique rows
print(check, n=Inf) # print all the rows
# looks good

# there are a few names that need editing
# MAPS differentiates the subspecies of yellow-rumped wabler (Myrtle and Audubon's)
# we have models with AUWA or Audubon's warbler. The trait files do not have this level of specificity
# Similarly, MAPS differentiates subspecies of dark-eyed junco and we have models for Oregon dark-eyed junco
# this next step changes the scientific name for these two to reflect the broader species name, rather than the subspecies name
names$SCINAME[names$SCINAME=="Setophaga coronata auduboni"]<-"Setophaga coronata"
names$SCINAME[names$SCINAME=="Junco hyemalis oreganus"]<-"Junco hyemalis"

# import conversion sheet from avonet (Tobias et al. 2022) that aligns bird tree, bird life and eBird taxonomies 
nameconvert <- read.csv(here("Data", "BirdNamesConversion.csv"), header=T)
colnames(nameconvert)

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

#################### Import DD Model Stats ##################

DDmod_stats <- readRDS(here("Outputs", "AverageNegativeSlopebySpecies_m2.rds"))
head(DDmod_stats)
# our two response variables for each species are min_adult and estimate (average of negative slopes)

#################### Join model statistics and species names ##################

DDnames_stats <- left_join(DDmod_stats, DDspp_allnames, by = "SPEC") %>%
  select(SPEC, min_adult, estimate:Species3_BirdTree)
head(DDnames_stats)

############## Import trait files ################################################

#### diet/trophic and morphometric trait files #####
avonet <- read.csv(here("Data", "avonet.csv"), header=T) # Tobias et al. 2022 (avonet)
elton <- read.csv(here("Data", "elton.csv"), header=T) # Wilman et al. 2014
pigot <- read.csv(here("Data", "pigot.csv"), header=T) # Pigot et al. 2020

#### life history trait files ####  
longevity <- read.csv(here("Data", "longevity.csv"), header=T) # Bird et al. 2020 (annual adult survival trait)

amniote <- read.csv(here("Data", "amniote.csv"), header=T) # Myhrvold et al. 2015 (egg/clutch traits)
clutch <- amniote %>% 
  filter(class=="Aves") %>% # filter data to retain info for birds
  mutate(scientific_name = paste(genus, species, sep=" ")) %>% # there is no column with full sci name (genus and species), so creating one
  mutate(across(where(is.numeric), ~na_if(., -999))) # this spreadsheet uses -999 if there is missing data. Changing missing values to NA for now
colnames(clutch)
head(clutch) # scientific_name is the sci name with genus and species separated by a space
# note that there are species with missing values in this database

### sexual selection trait files ###
delhey <- read.csv(here("Data", "delhey.csv"), header=T) # Delhey et al. 2023
# Dehley et al. 2023 has scores that reflect the intensity of sexual selection on males and females 

# plumage scores from males and females that we will use to calculate sexual dichromatism
plumage <- read.csv(here("Data", "plumage.csv"), header=T) # Dale et al. 2015

##########################################################################################
############################### Join traits with species results #########################
##########################################################################################

##########################################################################################
#### trophic/diet traits #####

# elton traits from Wilman et al. 2014
# this contains diet percentages by different food groups and foraging strata use for each species
head(elton)
mapsDD_elton_trophic <- elton %>%
  rename(Species3_BirdTree = Scientific) %>% # rename column as Species3_BirdTree to facilitate join
  select(Species3_BirdTree, Diet.Inv:ForStrat.SpecLevel) %>% # keep only traits of interest
  left_join(DDnames_stats, ., by="Species3_BirdTree") # join modified elton traits data frame and DDnames_stats data frame

# most species have diet information with certainty = A, meaning high confidence/quality data
dietcert_A <- mapsDD_elton_trophic %>% filter(Diet.Certainty == "A")
print(dietcert_A)

# however, 3 species have diet information that does not have a high confidence associated with it
dietcert_C <- mapsDD_elton_trophic %>% filter(Diet.Certainty == "C")
print(dietcert_C)

# we updated these 3 species using current Birds of the World species pages
# import updated values
diet_update <- read.csv(here("Data", "eltondiet_update.csv"), header=T)

diet_update_spp <- dietcert_C %>% select(SPEC:Species3_BirdTree, ForStrat.watbelowsurf:ForStrat.SpecLevel) %>%
  left_join(., diet_update) %>%
  select(SPEC:Species3_BirdTree, Diet.Inv:Diet.EnteredBy, ForStrat.watbelowsurf:ForStrat.SpecLevel)

# combine update species with other species to get complete list
mapsDD_elton_update <- bind_rows(dietcert_A, diet_update_spp)

### Trophic trait 1: create a continuous trophic level score
# to do this, combine all diet categories that come from animals to reflect percentage of diet in each species that is animal-based
# plant diet categories include: Diet.Fruit, Diet.Nect, Diet.Seed, Diet.PlantO
# animal diet categories include: Diet.Inv, Diet.Vend, Diet.Vect, Diet.Vfish, Diet.Vunk, Diet.Scav
# use percentage of animals in diet to create a continuous score for each species that reflect its trophic level
# this score ranges between 0 and 1. Scores of 0 are 100% plant-based diets (herbivores)
# Scores of 1 are 100% animal-based diets (carnivores)
# Scores of 0.5 are 50/50 plants and animals (omnivores)
mapsDD_elton_trophiclevel <- mapsDD_elton_update %>%
  rowwise() %>%
  mutate(AnimalDietPerc = (Diet.Inv + Diet.Vend + Diet.Vect + Diet.Vfish + Diet.Vunk + Diet.Scav),
         TrophicLevel = (AnimalDietPerc/100))

range(mapsDD_elton_trophiclevel$TrophicLevel) # examine range of scores for our species
boxplot(mapsDD_elton_trophiclevel$TrophicLevel) # examine boxplot

### Trophic traits 2 and 3: diet diversity and foraging strata diversity scores

# use information in Diet columns and Strata columns to calculate scores
# Higher scores indicate more generalized diets and/or more diverse use of different foraging space
# Lower scores reflect more specialized diets/foraging space use
# Diet.score and Strata.score only count a category if it makes up 20% or more of the diet/foraging strata
# This decision was made because Wilman et al. 2014 state that a value of 10% for either diet or strata reflects very occasional or rare use
# However, if a species had multiple categories that are listed as 10%, they were combined
# For example, a species with 10% for fruit and 10% for nectar would get another point added to their score

mapsDD_elton_score <- mapsDD_elton_trophiclevel %>%
  rowwise() %>%
  mutate(Diet.score.1 = sum(c_across(Diet.Inv:Diet.PlantO) > 10), # get number of categories with 20% or higher
         Diet.score.2 = sum(c_across(Diet.Inv:Diet.PlantO) == 10), # get number of categories with 10%
         Diet.score = Diet.score.1 + floor(Diet.score.2/2), # Divide Diet.score.2 by two and round down to get total number of categories used by each species
         Strata.score.1 = sum(c_across(ForStrat.watbelowsurf:PelagicSpecialist) > 10), # get number of categories with 20% or higher
         Strata.score.2 = sum(c_across(ForStrat.watbelowsurf:PelagicSpecialist) == 10), # get number of categories with 10%
         Strata.score = Strata.score.1 + floor(Strata.score.2/2)) %>% # Divide Strata.score.2 by two and round down
  select(-Diet.score.1, -Diet.score.2,
         - Strata.score.1, - Strata.score.2) # remove these columns as they were intermediate steps to generate scores

# examine these scores (how are the species distributed across scores?)
mapsDD_elton_score %>% group_by(Diet.score) %>% count()
mapsDD_elton_score %>% group_by(Strata.score) %>% count()

# get final trophic traits data frame by simplifying mapsDD_elton_score
mapsDD_trophic <- mapsDD_elton_score %>% 
  select(SPEC:Species3_BirdTree, TrophicLevel, Diet.score, Strata.score)

##########################################################################################
#### morphometric traits #####

### Morphometric trait 1: hand wing index
# we will use hand wing index from avonet
# avonet uses BirdLife species names

head(avonet) # scientific names are stored in a column called Species1 with genus and spp separated by a space
head(DDnames_stats)
colnames(avonet)
# head() pulls the first 6 rows. colnames() prints the column names

# we need to change the column name in avonet to Species1_BirdLife to enable join
mapsDD_HWI <- avonet %>%
  rename(Species1_BirdLife = Species1) %>% # rename the column
  select(Species1_BirdLife, Hand.Wing.Index) %>% # keep only the trait of interest and drop the other columns
  left_join(DDnames_stats, ., by="Species1_BirdLife") # keep all rows in DDnames_stats (our 62 species) and add only the data from avonet for those species

# are any species missing values for HWI?
mapsDD_HWI %>% filter(is.na(Hand.Wing.Index)) %>% nrow() # keep any rows where HWI = NA (aka is missing)
# this should print 0 if all species have values

# examine the distribution of values for HWI
range(mapsDD_HWI$Hand.Wing.Index) # min = 11.4, max = 33.5
hist(mapsDD_HWI$Hand.Wing.Index) # look at a histogram of the values
boxplot(mapsDD_HWI$Hand.Wing.Index)


### Morphometric trait 2: body mass
# pull body mass from MAPS data
dat_10yr <- readRDS(here("Outputs", "dat_10yr.rds")) # import capture records with morphometrics

# reduce to the species in our analysis
DDspp_allnames # we will use DDspp_names to reduce data

mapsDDspp_morphs <- left_join(DDspp_allnames, dat_10yr) %>%
  filter(AgeCat == "Adult") %>%
  filter(BRSTAT %in% c("B", "U", "O")) %>% # keep only breeders 
  filter(STATUS %in% c("300", "301", "325")) %>% # keep only healthy banded birds
  filter(N == "-") %>% # keep only records to be used in productivity and survivorship analyses (may have already filtered for this)
  select(SPEC:STATION, SEX, WNG, WEIGHT, STATUS, N, BRSTAT, year, AgeCat)

# get average mass for each species and combine with HWI values from previous step
mapsDD_massHWI <- mapsDDspp_morphs %>%
  group_by(SPEC) %>%
  summarize(Mass=mean(WEIGHT, na.rm=T)) %>%
  right_join(mapsDD_HWI, .)


### Morphometric traits 3 and 4: beak PC1 and beak PC2 
# these come from Pigot et al. 2020
# PC1 for beak is relative beak size
# PC2 for beak is relative beak shape

mapsDD_pigot <- pigot %>%
  mutate(Species3_BirdTree = str_replace(Binomial, "_", " ")) %>% # rename column as Species3_BirdTree to facilitate join
  select(Species3_BirdTree, Beak_PC1, Beak_PC2) %>% # keep only traits of interest
  left_join(DDnames_stats, ., by="Species3_BirdTree") # join modified pigot traits data frame and DDnames_stats data frame

# combine with mass and HWI to get complete data frame of morphs
mapsDD_morphometrics <- left_join(mapsDD_massHWI, mapsDD_pigot)

##########################################################################################
#### life history traits #####

### life history trait 1: maximum longevity 
# Bird et al. 2020 uses BirdLife taxonomy
# we were originally interested in Maximum.longevity and Adult.survival from this data set, but they were strongly correlated
# proceeding with only longevity to avoid multicollinearity in models

mapsDD_longevity <- longevity %>% 
  rename(Species1_BirdLife = Scientific.name) %>% # Rename the species name column
  select(Species1_BirdLife, Maximum.longevity) %>% # Select max longevity
  left_join(DDnames_stats, ., by="Species1_BirdLife") # Join with DD stats by species name column 

# Check if there are missing values for max. longevity 
mapsDD_longevity %>% filter(is.na(Maximum.longevity)) %>% nrow() #There are no missing values 

# Examine max longevity
range(mapsDD_longevity$Maximum.longevity) # The range is 6.10 - 14.49
mean(mapsDD_longevity$Maximum.longevity) # The mean is 9.04 
hist(mapsDD_longevity$Maximum.longevity)

### Life history trait 2: clutch size
# from Myhrvold et al. 2015
head(clutch) 

# joining first using BirdLife scientific names, second using BirdTree names, and then joining both together to maximize matches
mapsDD_clutch_BL <- clutch %>% rename(Species1_BirdLife = scientific_name) %>% # rename columns to match MAPS
  left_join(DDnames_stats, ., by="Species1_BirdLife") %>% filter(!is.na(class))

mapsDD_clutch_BT <- clutch %>% rename(Species3_BirdTree = scientific_name) %>% 
  left_join(DDnames_stats, ., by="Species3_BirdTree") %>% filter(!is.na(class))

mapsDD_clutch <- full_join(mapsDD_clutch_BL, mapsDD_clutch_BT) %>%
  select(SPEC, COMMONNAME, Species1_BirdLife, Species2_eBird, Species3_BirdTree,
         litter_or_clutch_size_n)
nrow(mapsDD_clutch) # all 62 species

# combine longevity and adult survival with clutch traits
mapsDD_lifehistory <-
  left_join(mapsDD_longevity, mapsDD_clutch) 

##########################################################################################
### sexual selection traits #####

### SS trait 1: sexual size dimorphism
# use mapsDDspp_morphs from above
# these are measurements collected from breeding birds at MAPS stations
# separate by Sex and find degree of sexual dimorphism using wing length
mapsDD_sexdimorphism <- mapsDDspp_morphs %>%
  filter(SEX %in% c("M", "F")) %>% # keep only records with known sex
  group_by(SPEC, SEX) %>%
  summarize(Wing = mean(WNG, na.rm = T)) %>% # find mean wing length for each species/sex combination
  pivot_wider(., names_from = SEX, values_from = Wing) %>%
  rename(Wing_F = F, Wing_M = M) %>%
  rowwise() %>%
  mutate(Wing_DM = log(Wing_M) - log(Wing_F)) %>% # calculate dimorphism as log(male) minus log(female)
  left_join(DDnames_stats, .)

### SS trait 2: plumage sexual dichromatism
head(plumage) # uses BirdTree taxonomy

# join plumage scores and calculate sexual dichromatism
mapsDD_plumage <- plumage %>%
  rename(Species3_BirdTree = Scientific_name) %>% # rename column as Species3_BirdTree to facilitate join
  select(Species3_BirdTree, Female_plumage_score, Male_plumage_score) %>% # keep only traits of interest
  left_join(DDnames_stats, ., by="Species3_BirdTree") %>% # join modified plumage traits data frame and DDnames_stats data frame
  rowwise() %>%
  mutate(Plumage_DC = Male_plumage_score - Female_plumage_score) # calculate dichromatism

# are any missing?
mapsDD_plumage %>% filter(is.na(Plumage_DC))
# Yes. Two species: RBSA and DOWO

# combine with size dimorphism
mapsDD_sizeplumage <- left_join(mapsDD_sexdimorphism, mapsDD_plumage) %>%
  select(-Male_plumage_score, -Female_plumage_score, - Wing_F, -Wing_M)

### SS trait 3: mating system

# we are interested in sex.sel.m and sex.sel.f 
head(delhey) 
# phylo is the column that contains the species names. The genus and spp is separated by an underscore

mapsDD_delhey <- delhey %>%
  select(phylo,sex.sel.m, sex.sel.f) %>% # keep two columns only
  distinct() %>% # keep only one row per species. drop extra row associated with male vs female
  mutate(Species3_BirdTree = str_replace( phylo, "_", " ")) %>% # make a new column called Species3_BirdTree using the phylo column. Replace underscore with a space
  left_join(DDnames_stats, ., by="Species3_BirdTree") %>%
  select(-min_adult, -estimate, -conf.low, -conf.high, -n)


# export sexual selection intensity scores for Males and Females from Delhey et al. 2023
write.csv(mapsDD_delhey, here("Data", "mapsDD_delhey.csv")) # save the updated delhey data 


# look up and add missing species on Birds of the World
# Additionally, we looked up all species labeled as SSM = 0 and SSF = 0 on Birds of the World, which are monogamous species
# read section in the BOW species' account Behavior > Sexual Behavior to confirm mating system
# if there was evidence of non-monogamy described on BOW, we updated the mating system score to reflect the provided information

# import updated file
mapsDD_matingsystem <- read.csv(here("Data", "mapsDD_delhey_update.csv")) # read in updated delhey data 
# added 3 species to updated file (BUOR, SPTO and WOTH)
# 1 species updated - NOCA
#View(mapsDD_matingsystem)

# how many species per category of sex.sel.m?
mapsDD_matingsystem %>% group_by(sex.sel.m) %>% count()
# most species fall into sex.sel.m = 0 or 1

# how many species per category of sex.sel.f?
mapsDD_matingsystem %>% group_by(sex.sel.f) %>% count()
# most species fall into sex.sel.f = 0 

# are any species missing values?
mapsDD_matingsystem %>% filter(is.na(sex.sel.m))
mapsDD_matingsystem %>% filter(is.na(sex.sel.f))
# should not have any missing values because we updated the NAs 

# combine mating system with plumage dichromatism and size dimorphism scores
mapsDD_sexualselection <- left_join(mapsDD_sizeplumage, mapsDD_matingsystem)

View(mapsDD_sexualselection) 

# check degree of correlation between some of the variables to make sure they can be used in the same model
mapsDD_SScheck <- mapsDD_sexualselection %>% filter(!is.na(Plumage_DC))

# plumage dichromatism and size dimorphism
set.seed(285)
confintr::ci_cor(mapsDD_SScheck$Plumage_DC, mapsDD_SScheck$Wing_DM, method="spearman", type="bootstrap") # low

# plumage dichromatism and degree of polygyny in males
set.seed(597)
confintr::ci_cor(mapsDD_SScheck$Plumage_DC, mapsDD_SScheck$sex.sel.m, method="spearman", type="bootstrap") # low

# size dimorphism and degree of polygyny in males
set.seed(463)
confintr::ci_cor(mapsDD_SScheck$Wing_DM, mapsDD_SScheck$sex.sel.m, method="spearman", type="bootstrap") # low

# lastly, look at Beak PC2 and trophic level
# Both ended up showing the same trend in their respective trait models and it looks like Beak PC2 (in our 62 species) reflects trophic level
# use a correlation test to further explore this relationship
mapsDD_combo <- left_join(mapsDD_morphometrics, mapsDD_trophic)
set.seed(175)
confintr::ci_cor(mapsDD_combo$Beak_PC2, mapsDD_combo$TrophicLevel, method="spearman", type="bootstrap")
plot(mapsDD_combo$Beak_PC2, mapsDD_combo$TrophicLevel, xlab="Beak PC2", ylab = "Trophic Level")
# higher Beak PC2 scores (long beaks relative to their depth and width) occupy higher trophic levels (have a diet with a greater % of animal based sources)

##########################################################################################
####### Export objects as rds files to be used in trait pgls models ###########

saveRDS(mapsDD_morphometrics, here("Outputs", "mapsDD_morphometrics_m2.rds"))
saveRDS(mapsDD_lifehistory, here("Outputs", "mapsDD_lifehistory_m2.rds"))
saveRDS(mapsDD_trophic, here("Outputs", "mapsDD_trophic_m2.rds"))
saveRDS(mapsDD_sexualselection, here("Outputs", "mapsDD_sexualselection_m2.rds"))
