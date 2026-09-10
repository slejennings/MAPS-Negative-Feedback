# MAPS Negative Feedback README

############################################################################
### OVERVIEW
############################################################################

This repository contains data and code for the analyses in the manuscript "Traits shape demographic regulation in North American birds". The analyses are organized into an R Project, which will reproduce the results, tables, and figures presented in the manuscript. This Read Me file describes the required software and the organization of the R Project and the associated files.

Important: certain files are too large to sync with this GitHub repository and need to be downloaded from the accompanying Google Drive folder. This applies to any file marked with an asterisk (*) at the end of the filename. Link to folder: https://shorturl.at/VuObK

############################################################################
### CORRESPONDENCE
############################################################################

Please direct questions about the data, analysis, and results to:

Sarah L. Jennings, Email: sjenni02@calpoly.edu

Clinton D. Francis, Email: cdfranci@calpoly.edu

############################################################################
### SOFTWARE
############################################################################

Rstudio 2025.05.0+496 R version 4.4.3 (2025-02-28)

R Packages and versions: ape (5.8-1), bayesplot (1.15.0), brms (2.22.0), broom (1.0.7), broom.mixed (0.2.9.6), cmdstanr (0.8.1.9000), colorspace (2.1-1), confintr (1.0.2), geiger (2.0.11), ggdist (3.3.3.9000), ggeffects (2.2.0), ggExtra (0.11.0), ggforce (0.4.2), ggnewscale (2.2.0), ggspatial (1.1.9), ggtext (0.1.2), ggtree (4.3.0), ggtreeExtra (1.23.1), gridExtra (2.3), here (1.0.1), loo (2.8.0), marginaleffects (0.25.0), MetBrewer(0.2.0), patchwork (1.3.2), performance (0.15.1), Polychrome (1.5.4), posterior (1.6.1), priorsense (1.0.4), rstan(2.32.6), scales (1.4.0), sf (1.0-19), tidybayes (3.0.7), tidyverse (2.0.0)

############################################################################
### R PROJECT FOLDER STRUCTURE
############################################################################

The R Project contains folders that are organized using the following organizational structure:

Data

Figures

Models
- Trait Models

Outputs

Results
- Manuscript Tables

Scripts

############################################################################
### DESCRIPTIONS OF R PROJECT FOLDERS AND THEIR CONTENTS
############################################################################

### DATA

**Folder Contents:** 19 files, file format(s): .csv, .tre, .pdf

**Folder Description:** contains the data files that were used in the analysis, mostly stored as .csv. 

**Note:** files with * at end of name are too large to sync with the GitHub repository and need to be downloaded from the accompanying Google drive folder

**List of Files:**

***File name:*** aminote.csv

***File Description:*** Clutch traits from the Amniote life-history database on birds, mammals and reptiles (Myhrvold et al. 2015).

***Columns:***
The amniote life-history dataset along with descriptions of the various columns can be found in the following FigShare folder:
https://doi.org/10.6084/m9.figshare.c.3308127

-------------------------------------------------------------------------------------------------------------------

***File name:*** avonet.net

***File Description:*** AVONET database that contains species morphological, ecological, and geographical data for the World's birds (Tobias et al. 2020). We used this to obtain Hand-Wing Index values for species.

***Columns:***

The full AVONET dataset along with descriptions of the various columns can be found in the following FigShare folder:
https://figshare.com/s/b990722d72a26b5bfead

-------------------------------------------------------------------------------------------------------------------

***File name:*** BirdNamesConversion.csv

***File Description:*** 

Conversion between bird scientific names using different taxonomies. Adapted from AVONET (Tobias et al. 2019).

***Columns:***

Avibase.ID: alphanumeric code assigned by the Avibase database to each species

Species1_BirdLife: scientific name using the BirdLife International taxonomy

Species2_eBird: scientific name using the eBird/Clements taxonomy

Species3_BirdTree: scientific name using the Bird Tree (Jetz et al. 2012) taxonomy

-------------------------------------------------------------------------------------------------------------------

***File name:*** delhey.csv

***File Description:*** Sexual selection scores (degree of polygyny in males and degree of polyandry in females) for species from Delhey et al. 2023

***Columns:***

See the Supporting Information and specifically Dataset S01 that accompany Delhey et al. 2023 for descriptions of this dataset:
https://www.pnas.org/doi/10.1073/pnas.2217692120

-------------------------------------------------------------------------------------------------------------------

***File name:*** ebird_taxonomy_v2022.csv

***File Description:*** The 2022 version of the eBird/Clements Checklist. The eBird/Clements Checklist is an integrated global taxonomy for the birds of the world

***Columns:***

See information at this site for how to interpret the columns in this file: https://science.ebird.org/en/use-ebird-data/the-ebird-taxonomy/2022-ebird-taxonomy-update

-------------------------------------------------------------------------------------------------------------------
***File name:*** elton.csv

***File Description:*** Species foraging and diet traits from EltonTraits 1.0 (Wilman et al. 2014)

***Columns:***

See the following FigShare folder for descriptions of the EltonTraits 1.0 dataset:

Wilman, Hamish; Belmaker, Jonathan; Simpson, Jennifer; de la Rosa, Carolina; Rivadeneira, Marcelo M.; Jetz, Walter (2016). EltonTraits 1.0: Species-level foraging attributes of the world's birds and mammals. Wiley. Collection. https://doi.org/10.6084/m9.figshare.c.3306933.v1

-------------------------------------------------------------------------------------------------------------------

***File name:*** eltondiet_update.csv

***File Description:*** 

Contains updated diet information for 3 of the species dataset. These species had low confidence associated with the diet information provided in Elton Traits 1.0 database (Wilman et al. 2014) noted by a score of "C" in the Diet-Certainty column. We researched their diets using the Birds of the World Database (Billerman et al. 2026) and updated their scores for the various diet categories. This sheet contains the updated information

***Columns:***

SPEC: species 4-letter code 

COMMONNAME: English common name

Diet-Inv: Percentage of diet made up of invertebrates. In increments of 10%

Diet-Vend: Percentage of diet made up of endothermic vertebrates. In increments of 10%

Diet-Vfish:  Percentage of diet made up of fish. In increments of 10%

Diet-Vunk: Percentage of diet made up of unknown vertebrates. In increments of 10%

Diet-Scav: Percentage of diet made up of scavenged items (e.g., carrion, offal, carcasses, garbage). In increments of 10%

Diet-Fruit: Percentage of diet made up of fruit. In increments of 10%

Diet-Nect: Percentage of diet made up of nectar. In increments of 10%

Diet-Seed: Percentage of diet made up of seeds. In increments of 10%

Diet-PlantO: Percentage of diet made up of other plant material (not accounted for in other categories). In increments of 10%

Diet-5Cat: Assigns each species to one of 5 diet groups (PlantSeed, FruiNect, Invertebrate, VertFishScav or Omnivore)

Diet-Source: Specifies the source(s) used to assign the species' diet information

Diet-Certainty: Contains the confidence code to reflect the reliability of species assignments

Diet-EnteredBy: Initials of person who entered the data

-------------------------------------------------------------------------------------------------------------------
***File name:*** IBP-AOS-LIST23.csv

***File Description:*** Institute for Bird Populations (IBP) standardized alpha bird species code list updated for 2023

***Columns:***

SPEC: 4 letter species code

COMMONNAME: English common name

SCINAME: Scientific name

SPEC6: 6 letter species code

-------------------------------------------------------------------------------------------------------------------

***File name:*** Jetz_ConsensusPhy.tre

***File Description:*** consensus phylogeny based on 10,000 for 9,993 bird species based on Jetz et al. 2012

-------------------------------------------------------------------------------------------------------------------

***File name:*** longevity.csv

***File Description:*** simplified version of species life history traits from Bird et al. 2020

***Columns:***

Scientific.name: species scientific name

Adult.survival: the estimated adult survival rate for each species

Age.at.first.reproduction: the age of first reproductive event for each species in years

Maximum.longevity: the maximum lifespan for each species in years

-------------------------------------------------------------------------------------------------------------------

***File name:*** MAPS-database_codes-and-structures.pdf

***File Description:*** Explains the MAPS program and the contents and structure of all of the subsequent files in this folder that have MAPS in their filename

-------------------------------------------------------------------------------------------------------------------

***File name:*** MAPS_BANDING_capture_data.csv*

***File Description:*** Contains all the banding records for captured birds

***Columns:*** 

Refer to “MAPS-database_codes-and-structures.pdf” for a complete description of the columns and data within this file

-------------------------------------------------------------------------------------------------------------------

***File name:*** MAPS_EFFORT_net_open_and_close_times.csv*

***File Description:*** Contains information about each station's operation, including the time each mist net was opened and closed

***Columns:*** 

Refer to “MAPS-database_codes-and-structures.pdf” for a complete description of the columns and data within this file

-------------------------------------------------------------------------------------------------------------------

***File name:*** MAPS_missing_effort.csv

***File Description:*** Certain stations/years were missing effort information in MAPS_EFFORT_net_open_and_close_times.csv. This file was provided to us by MAPS staff to fill in the missing values. Its structure is similar to MAPS_EFFORT_net_open_and_close_times.csv

***Columns:*** 

Refer to “MAPS-database_codes-and-structures.pdf” for a complete description of the columns and data within this file

-------------------------------------------------------------------------------------------------------------------

***File name:*** MAPS_STATION_location_and_operations.csv

***File Description:*** Contains information about each MAPS station, including geographic coordinates and years of operation

***Columns:***

Refer to “MAPS-database_codes-and-structures.pdf” for a complete description of the columns and data within this file

-------------------------------------------------------------------------------------------------------------------

***File name:*** mapsNF_delhey_update.csv

***File Description:*** 

We exported the sexual selection scores obtained from Delhey et al. 2023 for the 62 species in our analysis, which are contained in mapsNF_delhey.csv. If a species had scores of 0 for sex.sel.m and 0 for sex.sel.f, we investigated it using the species accounts on Birds of the World (Billerman et al. 2026). These are species marked as monogamous and yet many species exhibit some degree of non-monogamy when carefully studied. If the Birds of the World species account contain more recent information to suggest an alternative mating system, we updated their score. This sheet contains the information for any updated species

***Columns:***

SPEC: species 4 letter code

COMMONNAME: English common name

Species1_BirdLife: Scientific name using BirdLife International taxonomy

Species2_eBird: Scientific name using eBird/Clements taxonomy

Species3_BirdTree: Scientific name using Bird Tree taxonomy (Jetz et al. 2012)

sex.sel.m: Sexual selection score that reflects the degree of polygyny in males. Ranges from 0 to 5

sex.sel.f: Sexual selection score that reflects the and degree of polyandry in females. Ranges from 0 to 5

-------------------------------------------------------------------------------------------------------------------

***File name:*** mapsNF_delhey.csv

***File Description:*** 

We exported the sexual selection scores obtained from Delhey et al. 2023 for the 62 species in our analysis. This sheet contains that information. If a species had scores of 0 for sex.sel.m and 0 for sex.sel.f, we investigated it using the species accounts on Birds of the World (Billerman et al. 2026). These are species marked as monogamous and yet many species exhibit some degree of non-monogamy when carefully studied. If the Birds of the World species account contain more recent information to suggest an alternative mating system, we updated their score (see mapsNF_delhey_update.csv)


***Columns:***

SPEC: species 4 letter code

COMMONNAME: English common name

Species1_BirdLife: Scientific name using BirdLife International taxonomy

Species2_eBird: Scientific name using eBird/Clements taxonomy

Species3_BirdTree: Scientific name using Bird Tree taxonomy (Jetz et al. 2012)

sex.sel.m: Sexual selection score that reflects the degree of polygyny in males. Ranges from 0 to 5

sex.sel.f: Sexual selection score that reflects the and degree of polyandry in females. Ranges from 0 to 5


-------------------------------------------------------------------------------------------------------------------

***File name:*** pigot.csv

***File Description:*** Species morphological PC scores from Pigot et al. 2020

***Columns:***

See Supplementary Information (and specifically Dataset 1) available at the following link for descriptions of this dataset: https://www.nature.com/articles/s41559-019-1070-4

-------------------------------------------------------------------------------------------------------------------

***File name:*** plumage.csv

***File Description:*** species male and female plumage coloration scores from Dale et al. 2015

***Columns:***

Please see the following Dryad folder where the data associated with Dale et al. (2015) are archived for how to interpret the columns in this .csv file: https://datadryad.org/dataset/doi:10.5061/dryad.1rp0s

############################################################################

### SCRIPTS

**Folder Contents:** 7 files, file format(s): .R

**Folder Description:** contains R scripts to analyze the data and to generate the findings in the manuscript. The scripts are sequential, should be run in numerical order, and are labeled accordingly. At the beginning of each R script is a section titled “Objectives/Description of the Script” which details the main tasks accomplished by each file. We refer you to this section in lieu of providing an additional description of each file here.

**List of Files:**

Step1_DataPrepforNegativeFeedbackModels.R

Step2_DetermineCustomPriors.R

Step3_PriorSensitivityTesting.R

Step4_ModelComparison.R

Step5_EvaluateNFModels.R

Step6_JoinTraitswithSpeciesResults.R

Step7_TraitAnalyses.R

############################################################################

### OUTPUTS

**Folder Contents:** X files, file format(s): .rds, .csv

**Folder Description:** Contains various files stored as .rds or .csv that contain organized data that is being moved between R scripts and used in subsequent analysis steps. All of these objects are generated using the scripts, so we have not described them individually here.

**List of Files:**

AverageNegativeSlopesbySpecies_m2.rds 

compare_paretoK_m1m2m3.rds

dat_10yr.rds

dat_LH.rds

dat_morph.rds

dat_ss.rds

dat_trophic.rds

NF_breed_dat.rds

NFModelEffects_m2.csv

ELPD_m1_allspp.rds

ELPD_m2_allspp.rds

ELPD_m3_allspp.rds

epred_df_m2.rds

mapsNF_lifehistory_m2.rds

mapsNF_morphometrics_m2.rds

mapsNF_sexualselection_m2.rds

mapsNF_trophic_m2.rds

priors1_priorsensitivity.rds

priors2_priorsensitivity.rds

priors3_priorsensitivity.rds

SlopeSummarybySpecies_m2.rds

SppNames_STA.rds

STA_RE_m2.rds


############################################################################

### MODELS

Models Folder Structure (contains one subfolders):

Models
  - Trait Models


**Folder:** Models

**Folder Contents:** 12 files, files format(s): .rds

**Folder Description:** contains files related to the species-specific models that are produced by the R Scripts. Models in this folder were run using different priors, different random effects structures, and also include steps to obtain accurate ELPD values from models to be used in model comparison. As all files are generated by the accompanying code, we have not provided individual descriptions. Files are exported and saved as .rds

**List of Files:**

NOTE: all these files are too large to sync with GitHub and need to be downloaded form the Google Drive folder linked above

m_allspp_priors1.rds*

m_allspp_priors2.rds*

m_allspp_priors3.rds*

m1_allspp_loo.rds*

m1_momentmatch.rds*

m1_reloo.rds*

m2_allspp_loo.rds*

m2_momentmatch.rds*

m2_reloo.rds*

m3_allspp_loo.rds*

m3_momentmatch.rds*

m3_reloo.rds*


**Folder:** Trait Models

**Folder Contents:** subfolder within Models, 8 files, files format(s): .rds

**Folder Description:** contains the trait-focused phylogenetic models that are produced by the R scripts. As all files are generated by the accompanying code, we have not provided individual descriptions. Files are exported and saved as .rds.

**List of Files:**

minadult_LH_mod_nb.rds

minadult_Morph_mod_nb.rds

minadult_SS_mod_nb.rds

minadult_Trophic_mod_nb.rds

slope_LH_mod_logn.rds

slope_Morph_mod_logn.rds

slope_SS_mod_logn.rds

slope_Trophic_mod_logn.rds

############################################################################

### FIGURES

**Folder Contents:** 26 files, file format(s): .pdf, .R

**Folder Description:** contains copies of the figures produced and included in the accompanying publication, as well as plots used to examine convergence and validate of species-specific models. The R scripts are named to reflect the figure(s) they generate (e.g., Figure1.R generates Figure 1 in the manuscript). Each script has an accompanying .pdf file of the same name, which is the exported version of the final figure(s). For example, Figure1.R generates Figure1.pdf. Additionally, at the beginning of each R script is a section titled “Objectives/Description of the Script” which details the main tasks accomplished by each file. We refer you to this section in lieu of providing an additional description of each R script here. All remaining .pdf files are generated by the code in the Scripts folder (see Step2_DetermineCustomPriors.R and Step5_EvaluateNFModels.R), so we have not provided individual descriptions

**List of Files:**

autocorrelationplots_m2.pdf

conditionaleffectplots_m2.pdf

Figure1.pdf

Figure1.R

Figure2.pdf

Figure2.R

Figure3.pdf

Figure3.R

Figure4.pdf

Figure4&S2.R

FigureS2.pdf

FigureS1.R

FigureS1_Type1_Page1.pdf

FigureS1_Type2_Page1.pdf

FigureS1_Type2_Page2.pdf

FigureS1_Type2_Page3.pdf

FigureS1_Type2_Page4.pdf

FigureS1_Type3_Page1.pdf

logtransformplots.pdf

ppc_density_plots_m2.pdf

ppc_ecdfoverlay_m2.pdf

ppc_histograms_m2.pdf

ppc_intervalplots_m2.pdf

productivityplots.pdf

tracedensityplots_m2.pdf

violinplots_m2.pdf

############################################################################

### RESULTS

Results Folder Structure (contains one subfolders):

Results
  - Manuscript Tables


**Folder:** Results

**Folder Contents:** 10 files, files format(s): .csv

**Folder Description:** contains files related to results of the trait-focused phylogenetic models that are produced by the R Scripts. As all files are generated by the accompanying code, we have not provided individual descriptions. Files are exported and saved as .csv

**List of Files:**

all_traitmodels_summary.csv

NFCurveClassificationbySpecies.csv

MinAdultLH.csv

MinAdultMorph.csv

MinAdultSS.csv

MinAdultTrophic.csv

SlopeLH.csv

SlopeMorph.csv

SlopeSS.csv

SlopeTrophic.csv


**Folder:** Manuscript Tables

**Folder Contents:** subfolder within Results, 4 files, files format(s): .xlxs

**Folder Description:** Contain copies of formatted tables that are including in the main text and/or supplementary files of the associated manuscript. As all files are presented and fully explained in the the accompanying manuscript, we have not provided individual descriptions. Files are saved as .xlsx

**List of Files:**

Table1_TraitPredictions.xlsx

TableS1_Priors.xlsx

TableS2_ELPDModelComparisions.xlsx

TableS3_TraitModelResults.xlsx

############################################################################
### REFERENCES
############################################################################

Billerman, SM, BK Keeney, GM Kirwan, F Medrano, ND Sly, and MG Smith, Editors (2026). Birds of the World. Cornell Laboratory of Ornithology, Ithaca, NY, USA. https://doi.org/10.2173/bow

Bird, JP, R Martin, HR Akçakaya, J Gilroy, IJ Burfield, ST Garnett, A Symes, J Taylor, ÇH Şekercioğlu, and SHM Butchart. 2020. Generation Lengths of the World’s Birds and Their Implications for Extinction Risk. Conservation Biology 34 (5): 1252–61. https://doi.org/10.1111/COBI.13486

Dale, D, CJ Dey, K Delhey, B Kempenaers, and M Valcu. 2015. The effects of life history and sexual selection on male and female plumage coloration. Nature 527(7578): 367-70. https://doi.org/10.1038/nature15509

Delhey K, M Valcu, C Muck, B Kempenaers. 2023. Evolutionary predictors of the specific colors of birds. Proceedings of the National Academy of Sciences 120 (34): e2217692120. https://doi.org/10.1073/pnas.2217692120

Jetz, W, GH Thomas, JB Joy, K Hartmann, and AO Mooers. 2012. The Global Diversity of Birds in Space and Time. Nature 491 (7424): 444–48. https://doi.org/10.1038/nature11631

Myhrvold, NP, E Baldridge, B Chan, D Sivam, DL Freeman, and SK Morgan Ernest. 2015. An Amniote Life-History Database to Perform Comparative Analyses with Birds, Mammals, and Reptiles. Ecology 96 (11): 3109–3109. https://doi.org/10.1890/15-0846R.1

Pigot, AL, C Sheard, ET Miller, TP Bregman, BG Freeman, U Roll, N Seddon, CH Trisos, BC Weeks, and JA Tobias. 2020. Macroevolutionary convergence connects morphological form to ecological function in birds. Nature Ecology & Evolution 4:230–239. https://doi.org/10.1038/s41559-019-1070-4

Tobias, JA, C Sheard, AL Pigot, AJM Devenish, J Yang, F Sayol, MHC Neate-Clegg, et al. 2022. AVONET: Morphological, Ecological and Geographical Data for All Birds. Ecology Letters 25 (3): 581–97. https://doi.org/10.1111/ele.13898

Wilman, H, J Belmaker, J Simpson, C de la Rosa, MM Rivadeneira, and W Jetz. 2014. EltonTraits 1.0: Species-Level Foraging Attributes of the World’s Birds and Mammals. Ecology 95 (7): 2027–2027. https://doi.org/10.1890/13-1917.1

