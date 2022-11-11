
# analyse the Vile et al. (2006) data

# load relevant libraries
library(readr)
library(dplyr)
library(tidyr)
library(here)
library(rethinking)

# load the species RGR values from Vile et al. 2006
vile_rgr <- read_csv(here("data/Vile_et_al_2006_Ecology_Letters_RGR.csv"))
head(vile_rgr)

# fix the names
vile_rgr <- 
  vile_rgr %>%
  mutate(Species = gsub(pattern = " ", replacement = "_", x = Species))

# load the relative abundance data
vile_ra <- read_csv(here("data/Vile_et_al_2006_Ecology_Letters_RA.csv"))
head(vile_ra)

# get the first 13 columns
vile_ra <- vile_ra[,1:13]

# process the relative abundance data into a usable format i.e. site x sp
vile_ra <- 
  vile_ra %>%
  mutate(Species = gsub(pattern = " ", replacement = "_", x = Species)) %>%
  pivot_longer(cols = names(vile_ra)[-1],
               names_to = "field_age",
               values_to = "live_biomass") %>%
  arrange(field_age, Species) %>%
  mutate(field_age = substr(field_age, start = 1, stop = 2)) %>%
  mutate(field_age = gsub(pattern = "\\.", replacement = "", x = field_age)) %>%
  mutate(field_age = as.integer((field_age)))
  
# add a plot id column
vile_ra <- bind_cols(tibble(id = rep(1:12, each = length(unique(vile_rgr$Species))) ), 
                     vile_ra)

# get the two most dominant species per plot id
vile_ra <- 
  vile_ra %>%
  group_by(id) %>%
  mutate(q12 = live_biomass[rank(live_biomass) == max(rank(live_biomass))],
         q11 = live_biomass[rank(live_biomass) == (max(rank(live_biomass))-1) ]) %>%
  filter(live_biomass == q12 | live_biomass == q11) %>%
  ungroup() %>%
  select(-q12, -q11) %>%
  arrange(field_age, Species)

# how many species are there?
unique(vile_ra$Species)

# load the SANPP data
vile_SANPP <- read_csv(here("data/Vile_et_al_2006_Ecology_Letters_SANPP.csv"))
head(vile_SANPP)

# add a plot id column
vile_SANPP <- bind_cols(tibble(id = 1:nrow(vile_SANPP)), vile_SANPP)
head(vile_SANPP)

# reduce the number of species to the two most dominant per field






