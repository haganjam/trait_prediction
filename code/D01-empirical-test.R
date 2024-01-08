
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
print(cor_test)
(cor_test$estimate^2)

# make the plot
p1 <- 
  ggplot(data = SLA_RGR,
       mapping = aes(x = SLA, y = RGRmax)) +
  geom_point(size = 2.5, alpha = 0.75) +
  ylab("RGR max (g g-1 day-1)") +
  xlab("Specific leaf area (SLA, cm g-1)") +
  annotate(geom = "text", 
           label = paste0("r = ", round(cor_test$estimate, 2), " ; ",
                          "P = ", round(cor_test$p.value, 2)),
           x = Inf, y = -Inf, vjust = -1, hjust = 1.1, size = 4) +
  scale_y_continuous(limits = c(0.16, 0.34)) +
  theme_meta() +
  theme_transparent()
plot(p1)

# export the figure for further modification
ggsave(filename = "figures-tables/fig6.pdf", p1,
       unit = "cm", width = 10, height = 8, bg = "transparent")

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
lm1_s <- summary(lm1)
print(lm1_s)
AIC(lm1)

# fit a linear model to these data
lm2 <- lm(SANPP ~ SLA_CWM, data = sanpp_data)
plot(lm2)

# check the model summary
lm2_s <- summary(lm2)
print(lm2_s)
AIC(lm2)

# these two models as a comparison
p2 <- 
  ggplot(data = sanpp_data,
         mapping = aes(x = RGRmax_CWM, y = SANPP)) +
  geom_smooth(method = "lm", size = 0.5, alpha = 0.2, fill = "darkred", colour = "darkred") +
  geom_point(shape = 16, alpha = 0.6, size = 2.5, colour = "darkred") + 
  scale_colour_viridis_d() +
  xlab("CWM RGR (g g-1 day-1)") +
  ylab("SANPP (g kg-1 day-1)") +
  annotate(geom = "text", 
           label = lm_eqn(lm1),
           x = Inf, y = -Inf, vjust = -1.25, hjust = 1.25, size = 4.5, parse = TRUE) +
  annotate(geom = "text", 
           label = paste0("AIC = ", round(AIC(lm1), 1)),
           x = Inf, y = -Inf, vjust = -4, hjust = 1.23, size = 4.5) +
  theme_meta() +
  theme(legend.position = "none") +
  theme_transparent()
plot(p2)

# export the figure for further modification
ggsave(filename = "figures-tables/fig7.pdf", p2,
       unit = "cm", width = 12, height = 10, bg = "transparent")

# these two models as a comparison
p3 <- 
  ggplot(data = sanpp_data,
         mapping = aes(x = SLA_CWM, y = SANPP)) +
  geom_smooth(method = "lm", size = 0.5, alpha = 0.2, fill = "darkblue", colour = "darkblue") +
  geom_point(shape = 16, alpha = 0.6, size = 2.5, colour = "darkblue") + 
  scale_colour_viridis_d() +
  xlab("CWM SLA (cm g-1)") +
  ylab("SANPP (g kg-1 day-1)") +
  annotate(geom = "text", 
           label = lm_eqn(lm2),
           x = Inf, y = -Inf, vjust = -1.25, hjust = 1.25, size = 4.5, parse = TRUE) +
  annotate(geom = "text", 
           label = paste0("AIC = ", round(AIC(lm2), 1)),
           x = Inf, y = -Inf, vjust = -4, hjust = 1.23, size = 4.5) +
  theme_meta() +
  theme(legend.position = "none") +
  theme_transparent()
plot(p3)

# export the figure for further modification
ggsave(filename = "figures-tables/fig8.pdf", p3,
       unit = "cm", width = 12, height = 10, bg = "transparent")



