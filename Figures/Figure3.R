###### MAPS Project: Density Dependence #######
### Script name: Figure3.R
### Author(s): SLJ

########### Objective/Description of Script #####################
# create Figure 3, which contains two panels
# Panel 1: distribution for threshold and intensity and the correlation between them
# Panel 2: phylogeny of 62 species with heatmap showing species' values for threshold and intensity
#################################################################

### Setup ###

library(here)
library(ape)
library(geiger)
library(ggtree)
library(ggtreeExtra)
library(ggExtra)
library(ggnewscale)
library(colorspace)
library(Polychrome)
library(tidyverse)
library(patchwork)

#########################################################################
### Load required files ###

# import data
DD_breed_dat <- readRDS(here("Outputs", "DD_breed_dat.rds"))

ave_negative_slope_m2 <- readRDS(here("Outputs", "AverageNegativeSlopebySpecies_m2.rds"))

# import list of species with names and classifications of curve type and the MAPS stations for each species
SppNames_STA <- readRDS(here("Outputs", "SppNames_STA.rds"))

# import bird list used by MAPS to convert 4-letter bird codes to scientific names
birdcodes <- read.csv(here("Data", "IBP-AOS-LIST23.csv"), header=T) 
head(birdcodes)  
# downloaded from IBP website: https://www.birdpop.org/pages/birdSpeciesCodes.php

# import conversion sheet from avonet (Tobias et al. 2022) that aligns bird tree, bird life and eBird taxonomies 
nameconvert <- read.csv(here("Data", "BirdNamesConversion.csv"), header=T)
colnames(nameconvert)

# import eBird taxonomy to add Family info for each species
ebird_tax <- read.csv(here("Data","ebird_taxonomy_v2022.csv"), header=T) 
colnames(ebird_tax)

# import phylogenetic tree
tree_out <- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

#########################################################################

### Create Phylogeny Figure ###

# get a list of the 62 species 
DD_SPEC <- DD_breed_dat %>%
  ungroup() %>%
  select(SPEC) %>%
  distinct()

# combine names with species codes in DD_SPEC
names <- DD_SPEC %>% 
  left_join(., birdcodes, by="SPEC") %>% # left join keeps all the rows in DD_breed_dat and adds anything that matches from birdcodes
  dplyr::select(-SP, -CONF, -SPEC6, -CONF6) # remove some unnecessary columns

nrow(names) # check all species are present. Should equal 62

# check all 4 letter codes were paired with a name and that there are no duplicate names
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

# confirm we have 62 species/rows
nrow(DDspp_names)

# one is missing -> BUOR
BUOR_name <- names %>%
  filter(SPEC == "BUOR") %>%
  rename(Species3_BirdTree = SCINAME) %>%
  left_join(., nameconvert) %>%
  dplyr::select(-Avibase.ID) # drop column with Avibase.ID

# join to get 62 species
DDspp_allnames <- bind_rows(DDspp_names, BUOR_name)

# confirm we have one-to-one match for scientific names between the 3 naming schemes
# these will print "TRUE" if we have one-to-one match
length(unique(DDspp_allnames$Species1_BirdLife)) == length(unique(DDspp_allnames$Species2_eBird)) 
length(unique(DDspp_allnames$Species1_BirdLife)) == length(unique(DDspp_allnames$Species3_BirdTree))

# add curve type for each species
DDspp_curve <- SppNames_STA %>% select(SPEC, Curve_Type) %>% distinct() %>%
  left_join(DDspp_allnames, .)

# add family information to bird species list using eBird taxonomy
DDspp_tax <- ebird_tax %>% 
  rename(Species2_eBird = SCI_NAME) %>% # use Species2_eBird for joining data as these are using eBird names
  select(Species2_eBird, FAMILY) %>%
  left_join(DDspp_curve, .) %>%
  mutate(Family = word(FAMILY, 1), # only keep first word in FAMILY columnn (scientific family name)
         Sci_Name = Species3_BirdTree,
         Spp_Name= str_replace(Sci_Name, " ", "_")) %>%
  select(SPEC, COMMONNAME, Family, Sci_Name, Spp_Name, Curve_Type)


# prep a list of species, their families, and curve types to be joined to phylogeny
families <- DDspp_tax %>% 
  select(Spp_Name, Family, Curve_Type) %>%
  rename(label = Spp_Name, family = Family, curve = Curve_Type)

# get vector of scientific names to prune the phylogeny
sppnames <- families %>% pull(label)

# prune phylogenetic tree so it only contains the 62 relevant species
tree_trim <- tidytree::keep.tip(tree_out, sppnames)

#  convert the phylogeny into tibble
phy_tibble <- as_tibble(tree_trim)

# add family information to the phylogeny
phy_join <- left_join(phy_tibble, families, by="label")

# convert back into class phylo
new_tree <- as.phylo(phy_join) 

# create a list of family names that are associated with each tree tip label (species name) 
family_info <- split(phy_join$label, phy_join$family)

# create a list of family names that are associated with each tree tip label (species name) 
curve_info <- split(phy_join$label, phy_join$curve)

# update the tree with the family as the grouping info using list created in previous step
family_tree <- groupOTU(new_tree, family_info, group_name="Family") 
# this will allow us to label sections of the tree based on family

# update the tree with the curve as the grouping info using list created in previous step
curvefamily_tree <- groupOTU(family_tree, curve_info, group_name="Curve")
# this will allow us to label tippoints using curve type

# join the species density dependence measures with the bird scientific names 
DDspp_dat <- left_join(DDspp_tax, ave_negative_slope_m2, by = "SPEC") %>%
  column_to_rownames(., var="Spp_Name")

# link data to phylogenetic tree
tree_dat <- treedata(curvefamily_tree, DDspp_dat, sort=T)

phytree <-tree_dat$phy # get phylo tree for plotting

# we don't want names with underscores in the plot, so create a df for converting them
newlabs <- DDspp_tax %>% select(Spp_Name, Sci_Name)
# rename the tips of the tree to scientific names without underscores
phytree_name <- treeio::rename_taxa(phytree, data=newlabs, key=Spp_Name, value=Sci_Name)
# note: ended up not plotting species names on the phylogeny but retaining the above steps for future reference

# we have 18 families represented by the 62 species. Create a 18-color distinct palette
# using package Polychrome to make this palette
set.seed(8502) # for reproducibility
palette18 <- createPalette(18, c("#00ffff", "#ff00ff", "#ffff00"), M=5000)
swatch(palette18) # view the palette
palette18_hex <- as.character(palette18) # save hex codes for colors in a vector of characters


# plot circular phylogeny with Family labels
# we need to see nodes on the tree to figure out where to apply family labels
ggtree::ggtree(phytree_name, layout="circular") +
    geom_tippoint(aes(color=Family), size=2, shape=19) +
    scale_color_manual(values = palette18_hex) +
    geom_text(aes(label=node)) # label nodes

# make a data frame of all families and their nodes that we want to add using geom_cladelabel
dt <- data.frame(node =c(123, 65, 70, 74, 84, 87, 92, 95, 102, 109),
                 familyname =c("Picidae", "Tyrannidae", "Vireonidae", "Paridae", "Troglodytidae", 
                               "Turdidae", "Fringillidae", "Cardinalidae", "Passerellidae", "Parilidae"),
                 vjust = c(1, 1.75, 2.5, 2, 2.25,
                           1, 0.75, -1.5, -1.5, -1),
                 hjust = c(0.5, 0.3, 0.1, 0.5, 0.9,
                           0.5, 0.6, 0.7, 0.8, 0.2))


# now plot curve type as the tippoint and add family labels using dt 
(circ <- ggtree::ggtree(phytree_name, layout="circular") +
    geom_tippoint(aes(color=Curve), size=2.5, shape=19) +
    scale_color_manual(values = c("#DA4167", "#29335C", "#F5AF00")) + # CHANGE COLORS
    geom_cladelab(data = dt, 
                  mapping=aes(node=node, label=familyname, vjust=vjust, hjust=hjust), offset = 30, offset.text = 30, fontsize=4.5) +
    guides(color = guide_legend( # to style this legend individually to have different settings that Intensity and Threshold legends, use guides()
    title = "Curve",
    position ="right",
    theme(legend.title.position = "top", # put legend title at top
          legend.title = element_text(size=12, margin=margin(b=8)), # change text for legend title
          legend.text = element_text(size=11),
          legend.key.size = unit(0.4, "cm")
          )))
)
          

# get data for threshold and intensity
# this needs to have species names formatted without the underscore to match the updated labels for tree tips

intensity_dat <- left_join(DDspp_tax, ave_negative_slope_m2, by = "SPEC") %>%
  mutate(intensity = abs(estimate)) %>%
  select(Sci_Name, intensity) %>%
  column_to_rownames(., var="Sci_Name")
  
threshold_dat <- left_join(DDspp_tax, ave_negative_slope_m2, by = "SPEC") %>%
  mutate(threshold = min_adult) %>%
  select(Sci_Name, threshold) %>%
  column_to_rownames(., var="Sci_Name")


# need to add heatmap for the two response variables using multiple steps
# we are adding to the circular phylogeny plot created above

# because we are adding two continous fill scales in this next set of steps, we need to use ggnewscale:new_scale_fill()
# to format the legends associated with each step, we need to put the guide arguments/details inside of scale_fill() otherwise we get an error message

# first add values for threshold


(thresholdplot <- gheatmap(circ, threshold_dat, offset=-2, width=.15, colnames =F) +
    scale_fill_continuous_sequential(palette = "Purples", name="Threshold", na.value="white", 
                                     guide= guide_colorbar(title = "Threshold", 
                                                           theme = theme(
                                                           legend.title = element_text(size=12, margin=margin(b=8)),
                                                           legend.text = element_text(size=11),
                                                           legend.key.size = unit(0.4, "cm")))))
                                                         
part1 <- thresholdplot + ggnewscale::new_scale_fill()

# now add values for intensity
(threshold_intensity_plot <- gheatmap(part1, intensity_dat, offset=10, width=.15, colnames = F) +
    scale_fill_continuous_sequential(palette = "DarkMint", name="Intensity", na.value="white", breaks = c(0.025, 0.075, 0.125),
    guide = guide_colorbar(title = "Intensity",
                           theme = theme(legend.title = element_text(size=12, margin=margin(b=8)),
                                         legend.text = element_text(size=11),
                                         legend.key.size = unit(0.4, "cm")))))
                                                     

#########################################################################

### Create Other Figure Panel ###

# pick colors from same palettes used for phylogeny

hcl_palettes(palette="DarkMint", n=9, plot=T) # look at palette used above for Intensity
sequential_hcl(9, "DarkMint") # get hex codes
# try "#3F8489". Darkest is "#0E3F5C" for outline


hcl_palettes(palette="Purples", n=9, plot=T) # look at palette used above for Threshold
sequential_hcl(9, "Purples") # get hex codes
# try "#8B7EBB", Darkest is "#3D1778" for outline


# combine threshold and intensity values into a single data frame
vars <- threshold_dat %>%
  rownames_to_column("Species") %>%
  left_join(intensity_dat %>% rownames_to_column("Species"), by="Species")
  
# make a scatterplot to show correlation between the two variables
(scatter <- ggplot(vars, aes(x=intensity, y=threshold)) +
  geom_point(size=2, color = "#484554") +
  theme_classic() +
  labs(x="Intensity", y="Threshold") +
  xlim(0, 0.15) + ylim(0,20) +
  theme(legend.position="none",
        axis.title = element_text(size=14),
        axis.text = element_text(size=12))
  )

# add histograms to sides of scatterplot with ggMarginal from ggExtra package
(panelA <- ggMarginal(scatter, type = "histogram",
           size = 2, # size of center scatterplot relative to histograms
           yparams = list(fill = "#8B7EBB", # set specific parameters for threshold histogram
                          col="#3D1778", # col is outline for histogram bars
                          bins=25), # number of bins in the histogram
           xparams = list(fill = "#3F8489",  # set specific parameters for intensity histogram
                          col= "#0E3F5C", 
                          bins=25))) 

#########################################################################

### Combine Figure Panels & Save ###

wrap_elements(plot_spacer() + panelA + plot_spacer() + plot_layout(widths=c(0.08, 0.55, 0.08))) / 
  wrap_elements(threshold_intensity_plot) +
  plot_layout(heights=c(7,10)) & plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size=14, family = "Arial", face="bold"))


