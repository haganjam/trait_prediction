#'
#' @title: Theoretical limits to the correlation between plant traits and productivity
#' 
#' @description: Here, we match the species names from Vile et al. (2006, Ecology Letters)
#' to names in the TRY database so that we can efficiently request the data.
#' 

# load relevant libraries
library(readr)
library(dplyr)
library(here)

# load the species names from Vile et al. 2006
vile_sp <- read_csv(here("data/Vile_et_al_2006_Ecology_Letters_RGR.csv"))
unique(vile_sp)

# load the try species list
try_sp <- read_tsv(here("data/TRY_species_list.txt"))
head(try_sp)

# only get the relevant columns from try
try_sp <- 
  try_sp %>%
  select(AccSpeciesID, AccSpeciesName) %>%
  rename(Species = AccSpeciesName)

# join the species names and the try database
vile_sp <- left_join(vile_sp, try_sp, by = "Species")

# pull the species names separated by a comma
sp_id <- unique(vile_sp$AccSpeciesID)
sp_id <- sp_id[!is.na(sp_id)]
paste(sp_id, collapse = ", ")

