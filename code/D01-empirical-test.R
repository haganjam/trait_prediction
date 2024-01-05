
# perform an empirical test

# load relevant scripts
source("code/helper-plotting-theme.R")

# load the trait data
trait_dat <- readRDS("data/vile-rgr-sla-data.rds")
head(trait_dat)

# plot the relationship between RGRmax and SLA across species
SLA_RGR <- 
  trait_dat |>
  dplyr::select(Species, RGRmax, SLA) |>
  dplyr::distinct()

# run a correlation analysis
cor_test <- cor.test(SLA_RGR$RGRmax, SLA_RGR$SLA)
(cor_test$estimate^2)

# make the plot
ggplot(data = SLA_RGR,
       mapping = aes(x = SLA, y = RGRmax)) +
  geom_point(size = 2.5, alpha = 0.75) +
  ylab("RGR max (g g-1 day-1)") +
  xlab("Specific leaf area (SLA, cm g-1)") +
  annotate(geom = "text", 
           label = paste0("r = ", round(cor_test$estimate, 2), " ; ",
                          "P = ", round(cor_test$p.value, 2)),
           x = Inf, y = -Inf, vjust = -1.25, hjust = 1.1, size = 5) +
  theme_meta()

# calculate relative abundance
trait_dat <- 
  trait_dat |>
  dplyr::group_by(id) |>
  dplyr::mutate(RA = live_biomass/sum(live_biomass)) |>
  dplyr::ungroup()

# calculate the aggregated RGR
trait_sum <- 
  trait_dat |>
  dplyr::group_by(id, field_age) |>
  dplyr::summarise(RGRmax_CWM = ( sum(RA*RGRmax) ),
                   RGRmax_FD = sum(RA*( (RGRmax-mean(RGRmax))^2 ) ),
                   SLA_CWM = ( sum(RA*SLA, na.rm = TRUE) ),
                   SLA_FD = sum(RA*( (SLA-mean(SLA, na.rm = TRUE))^2), na.rm = TRUE),
                   .groups = "drop")

# load the measured SANPP data
sanpp_data <- readRDS("data/vile-sanpp-data.rds")

# join the community weighted means and function diversity
sanpp_data <- dplyr::left_join(trait_sum, sanpp_data, by = c("id") )
summary(sanpp_data)

# remove field 11 years old
sanpp_data <- 
  sanpp_data |> 
  filter( !(id %in% c(1, 2)) )

# plot the bivariate correlation
plot(SANPP ~ RGRmax_CWM, sanpp_data)
plot(SANPP ~ SLA_CWM, sanpp_data)

# fit a linear model to these data
lm1 <- lm(SANPP ~ RGRmax_CWM, data = sanpp_data)
plot(lm1)

# check the model summary
summary(lm1)
AIC(lm1)

# fit a linear model to these data
lm2 <- lm(SANPP ~ SLA_CWM, data = sanpp_data)
plot(lm2)

# check the model summary
summary(lm2)
AIC(lm2)


