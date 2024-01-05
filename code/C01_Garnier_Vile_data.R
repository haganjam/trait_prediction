
# analyse the Vile et al. (2006) data

# load relevant libraries
library(readr)
library(dplyr)
library(tidyr)

# load the species RGR values from Vile et al. 2006
vile_rgr <- readr::read_csv("data/Vile_et_al_2006_Ecology_Letters_RGR.csv")
head(vile_rgr)

# fix the names
vile_rgr <- 
  vile_rgr |>
  dplyr::mutate(Species = gsub(pattern = " ", replacement = "_", x = Species))

# load the relative abundance data
vile_ra <- readr::read_csv("data/Vile_et_al_2006_Ecology_Letters_RA.csv")
head(vile_ra)

# get the first 13 columns
vile_ra <- vile_ra[,1:13]

# process the relative abundance data into a usable format i.e. site x sp
vile_ra <- 
  vile_ra |>
  dplyr::mutate(Species = gsub(pattern = " ", replacement = "_", x = Species)) |>
  tidyr::pivot_longer(cols = names(vile_ra)[-1],
                      names_to = "field_age",
                      values_to = "live_biomass") |>
  arrange(field_age, Species) |>
  dplyr::mutate(field_age = substr(field_age, start = 1, stop = 2)) |>
  dplyr::mutate(field_age = gsub(pattern = "\\.", replacement = "", x = field_age)) |>
  dplyr::mutate(field_age = as.integer((field_age)))
  
# add a plot id column
vile_ra <- dplyr::bind_cols(dplyr::tibble(id = rep(1:12, each = length(unique(vile_rgr$Species))) ), vile_ra)

# check if these values are relative abundance
vile_ra |>
  dplyr::group_by(field_age, id) |>
  dplyr::summarise(sum_live_biomass = sum(live_biomass))

# join these data
vile_dat <- 
  dplyr::full_join(vile_ra, 
                   dplyr::select(vile_rgr, Species, RGRmax), 
                   by = "Species")

# load the trait data
trait_dat <- readRDS("data/TRY_species_traits.rds")

# get SLA from this data
trait_dat <- 
  trait_dat |>
  dplyr::filter(Trait == "SLA") |>
  dplyr::select(AccSpeciesName, Trait_m) |>
  dplyr::rename(Species = AccSpeciesName, SLA = Trait_m) |>
  dplyr::mutate(Species = gsub(pattern = " ", "_", Species))
head(trait_dat)
range(trait_dat$SLA)

# check if all the species are present
vile_dat$Species[which( !(unique(vile_dat$Species) %in% trait_dat$Species) )]

# add the SLA data to the RGR and RA data
vile_dat <- dplyr::left_join(vile_dat, trait_dat, by = "Species")

# export these data
saveRDS(vile_dat, file = "data/vile-rgr-sla-data.rds")

# get the correct labels for the fields
ids <- 
  vile_dat |>
  dplyr::select(id, field_age) |>
  dplyr::distinct()

# load the SANPP data
vile_SANPP <- readr::read_csv("data/Garnier_2004_Ecology_ecosystem_properties.csv")
head(vile_SANPP)

# convert field age to a factor
vile_SANPP$Field_age <- factor(vile_SANPP$Field_age, levels = unique(ids$field_age) )
vile_SANPP <-
  vile_SANPP |>
  dplyr::arrange(Field_age)

# add the field ids
vile_SANPP$id <- ids$id

# reorder the columns
vile_SANPP <-
  vile_SANPP |>
  dplyr::select(id, field_age = Field_age, BIOmax, ANPP, SANPP) |>
  dplyr::mutate(field_age = as.integer(field_age))

# export as a .rds file
saveRDS(vile_SANPP, file = "data/vile-sanpp-data.rds")

### END
