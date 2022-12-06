
# analyse the TRY data for the species in Vile et al. (2006, Ecology Letters)

# load relevant libraries
library(readr)
library(dplyr)
library(tidyr)
library(here)

# load the species RGR values from Vile et al. 2006
try_dat <- read_tsv(here("data/TRY_data_23382.txt"))

# check basic data structure
head(try_dat)
names(try_dat)
str(try_dat)

# remove the empty column
try_dat <- 
  try_dat %>%
  select(-...29)

try_dat %>%
  filter(!is.na(TraitID)) %>%
  View()

summary(try_dat)
nrow(try_dat)

# standard trait filtering

# 1. only take observations with standardised trait values
try_dat <- 
  try_dat %>%
  filter(!is.na(StdValue))

# 2. extract all the rows with trait measurements (i.e. !is.na(TraitID))
try_trait <- 
  try_dat %>%
  filter(!is.na(TraitID))
summary(try_dat)

# 3. DataID = 327 (experimental treatment), 
# 413 (mature and juvenile plants), 1961 (unhealthy plants), 1444 (leaves)
try_trait %>%
  filter(DataID %in% c(327, 413, 1961, 1444)) %>%
  nrow()

# 4. assess error risk: Remove species with ErrorRisk greater than five
try_trait <- 
  try_trait %>%
  filter(ErrorRisk < 5) 

# access metadata of the geographical coordinates
try_meta <- 
  try_dat %>%
  filter(ObservationID %in% unique(try_trait$ObservationID)) %>%
  filter(is.na(TraitID)) %>%
  filter(DataName %in% c("Latitude", "Longitude", "Altitude")) %>%
  select(DatasetID, AccSpeciesID, AccSpeciesName, ObservationID,
         DataName, OrigValueStr, OrigUnitStr, StdValue, UnitName)
  
# put latitude, longitude and altitude into separate variables
try_meta <- 
  try_meta %>%
  pivot_wider(id_cols = c("DatasetID","AccSpeciesID","AccSpeciesName", "ObservationID"),
              names_from = c("DataName"),
              values_from = c("StdValue"))

# check the try_meta data
summary(try_meta)

# how many species are left?
unique(try_trait$AccSpeciesName)

# get the units for the different traits
try_trait %>%
  select(TraitID, UnitName) %>%
  distinct()

# select the relevant columns
try_trait <- 
  try_trait %>%
  select(DatasetID, AccSpeciesID, AccSpeciesName, ObservationID,
         TraitID, DataID, StdValue, UnitName)

# check if there is only one unique value per data entry
try_trait %>%
  group_by(DatasetID, AccSpeciesID, AccSpeciesName, ObservationID, TraitID, DataID ) %>%
  summarise(n = n()) %>%
  pull(n) %>%
  max()

# convert the traits into the wide format
try_trait <- 
  try_trait %>%
  pivot_wider(id_cols = c("DatasetID","AccSpeciesID","AccSpeciesName", "ObservationID", "DataID"),
              names_from = "TraitID",
              values_from = "StdValue")

# rename the traits from their IDs and specify their units

# 3115: Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole excluded
# unit: mm2 mg-1
# name: SLA1

# 109: Leaf area per plant dry mass (leaf area ratio; LAR)
# unit: mm2 mg-1
# name: LAR

# 3116:	Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole included
# unit: mm2 mg-1
# unit: SLA2

# 12: Leaf lifespan (longevity)
# unit: month
# name: L_life

# 3117: Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded
# unit: mm2 mg-1
# name: SLA3

# 144: Leaf length
# unit: mm
# name: L_len

# 660: Leaf nitrogen (N) content organic per leaf dry mass
# unit: mg/g
# name: L_N

# 3086: Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA) petiole, rhachis and midrib excluded
# unit: mm2 mg-1
# name: SLA4

try_trait <- 
  try_trait %>%
  rename(SLA1 = `3115`,
         SLA2 = `3116`,
         SLA3 = `3117`,
         SLA4 = `3086`,
         LAR = `109`,
         L_life = `12`,
         L_len = `144`,
         L_N = `660`) %>%
  select(DatasetID, AccSpeciesID, AccSpeciesName, ObservationID, DataID,
         SLA1, SLA2, SLA3, SLA4, LAR, L_life, L_len, L_N)
  
# join the latitude-longitude information
try_trait <- right_join(select(try_meta, ObservationID, Latitude, Longitude), 
                        try_trait, 
                        by = "ObservationID"
                        )

# filter by North America (lat-lon)?
# S - 26
# N - 49
# W - -125
# E - -68

# check how many species
unique(try_trait$AccSpeciesName)

# calculate the differences in missing value number between the different SLAs
lapply(try_trait[, grepl("SLA", names(try_trait))], function(x){
  length(x[!is.na(x)])
} )

# fill in SLA with that priority i.e. (abundance)
try_trait <- 
  try_trait %>%
  mutate(SLA = ifelse(!is.na(SLA3), SLA3, 
                      ifelse( !is.na(SLA2), SLA2, 
                              ifelse( !is.na(SLA1), SLA1,
                                      ifelse(is.na(SLA4), SLA4, NA) ) )  )  )

# collapse trait information by observation id
try_trait <- 
  try_trait %>%
  group_by(Latitude, Longitude, ObservationID, AccSpeciesID, AccSpeciesName) %>%
  summarise(across(.cols = c("SLA", "LAR", "L_life", "L_len", "L_N"), mean, na.rm = TRUE), .groups = "drop")

# get a list of observationIDs for the references
obs_id <- unique(try_trait$ObservationID)

# summarise by species
sp_traits <- 
  try_trait %>%
  pivot_longer(cols = c("SLA", "LAR", "L_life", "L_len", "L_N"),
               names_to = "Trait",
               values_to = "Value") %>%
  group_by(AccSpeciesID, AccSpeciesName, Trait) %>%
  summarise(Trait_m = mean(Value, na.rm = TRUE),
            Trait_sd = sd(Value, na.rm = TRUE),
            Trait_min = min(Value, na.rm = TRUE),
            Trait_max = max(Value, na.rm = TRUE),
            Trait_n = length(Value[!is.na(Value)]), .groups = "drop") %>%
  arrange(AccSpeciesID, AccSpeciesName, Trait) %>%
  filter(Trait_n > 0)

# warnings occur when there are no values for the traits in question but these values are removed  
# check the data
View(sp_traits)

# export a table of traits
saveRDS(sp_traits, file = here("data/TRY_species_traits.rds"))

### END
