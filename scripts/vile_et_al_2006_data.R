
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
un_sp <- unique(vile_ra$Species)

# fill in missing species
vile_ra <- 
  lapply(split(vile_ra, vile_ra$id), function(y) {
  
  df <- 
    tibble(id = y$id[1],
           Species = un_sp[!(un_sp %in% y$Species)],
           field_age = y$field_age[1],
           live_biomass = 0)
  
  df <- 
    bind_rows(df, y) %>%
    arrange(Species)
  
  return(df)
  
} )

vile_ra <- bind_rows(vile_ra)

# convert to relative biomass
vile_ra <- 
  vile_ra %>%
  group_by(id) %>%
  mutate(sum_biomass = sum(live_biomass)) %>%
  ungroup() %>%
  mutate(live_biomass = live_biomass/sum_biomass) %>%
  select(-sum_biomass)

# load the SANPP data
vile_SANPP <- read_csv(here("data/Vile_et_al_2006_Ecology_Letters_SANPP.csv"))
head(vile_SANPP)

# add a plot id column
vile_SANPP <- bind_cols(tibble(id = 1:nrow(vile_SANPP)), vile_SANPP)
head(vile_SANPP)

# join these data
vile_dat <- 
  left_join(vile_ra, 
             select(vile_rgr, Species, RGRmax), 
             by = "Species"
             )

vile_dat <- full_join(vile_dat, vile_SANPP, by = "id")

# can we get realistic SANPP data from the RGR values?
vile_dat %>%
  group_by(id) %>%
  summarise(SANPP_est = log10(sum( live_biomass*exp(RGRmax*6) ))/90,
            SANPP = median(SANPP)/1000)



# make a list to fit the productivity model
d <- 
  list(S = as.integer(as.factor(vile_dat$Species)),
       RGR = vile_dat$RGRmax,
       pi = vile_dat$live_biomass,
       NPP = vile_dat$SANPP/1000)

# write this as a proper stan model...

# useful looking thread: 
# https://discourse.mc-stan.org/t/finite-mixture-model-where-the-sum-of-a-groups-characteristics-is-known/20146
# https://discourse.mc-stan.org/t/sum-vector-by-groups-month-year/5588/5


