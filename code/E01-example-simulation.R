
# Example replicate

# load relevant functions
source("code/A01_sim_abun_func.R")
source("code/A02_sim_traits_func.R")
source("code/helper-plotting-theme.R")

# load relevant libraries
library(dplyr)
library(ggplot2)

# get traits from the sim_traits() function
df.traits <- 
  sim_traits(sp = 5, com = 100,
             t1 = "SLA", t2 = "RGR", 
             mu_t1 = 100, mu_t2 = 5, sd_t1 = 5, sd_t2 = 2.5, r = 0.5,
             Eind_t1 = 2.5, Eind_t2 = 0.5,
             t2_scale = 10)

# bind into a data.frame
df.traits <- dplyr::bind_rows(df.traits)

# plot the traits
p1 <- 
  ggplot(data = df.traits) +
  geom_point(mapping = aes(x = RGR, y = SLA, colour = sp), alpha = 0.75) +
  geom_density_2d(mapping = aes(x = RGR, y = SLA, colour = sp), alpha = 0.5, n = 100) +
  scale_colour_viridis_d() +
  ylab("Specific Leaf Area (SLA, cm g-1") +
  xlab("Relative growth rate (g g-1 day-1") +
  scale_x_continuous(limits = c(-0.1, 1.05)) +
  scale_y_continuous(limits = c(87, 115)) +
  theme_meta() +
  theme(legend.position = "none")
plot(p1)

p2 <- 
  ggplot(data = df.traits) +
  geom_point(mapping = aes(x = RGR, y = SLA, colour = sp), alpha = 0.75) +
  geom_density_2d(mapping = aes(x = RGR, y = SLA, colour = sp), alpha = 0.5, n = 100) +
  stat_ellipse(mapping = aes(x = RGR, y = SLA), linetype = "dashed") +
  scale_colour_viridis_d() +
  annotate(geom = "text", 
           label = "r = 0.5",
           x = -Inf, y = Inf, vjust = 2, hjust = -0.6, size = 5) +
  ylab("Specific Leaf Area (SLA, cm g-1)") +
  xlab("Relative growth rate (g g-1 day-1)") +
  scale_x_continuous(limits = c(-0.1, 1.05)) +
  scale_y_continuous(limits = c(87, 115)) +
  theme_meta() +
  theme(legend.position = "none")
plot(p2)

# simulate different communities from a log-normal distribution
df_even <- sim_abund(sp = 5, com = 1, N = 100, cv_abund = 0.5)

# bind into a data.frame
df_even <- dplyr::bind_rows(df_even)

# change the species names
df_even$sp <- as.character(1:5)

# plot the relative abundance distribution
p3 <- 
  ggplot(data = df_even,
       mapping = aes(x = sp, y = abund, colour = sp, fill = sp)) +
  geom_col(width = 0.5) +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  ylab("Relative abundance") +
  xlab("Species") +
  theme_meta() +
  theme(legend.position = "none")
plot(p3)

# simulate different communities from a log-normal distribution
df_uneven <- sim_abund(sp = 5, com = 1, N = 100, cv_abund = 20)

# bind into a data.frame
df_uneven <- dplyr::bind_rows(df_uneven)

# change the species names
df_uneven$sp <- as.character(1:5)

# plot the relative abundance distribution
p4 <- 
  ggplot(data = df_uneven,
       mapping = aes(x = sp, y = abund, colour = sp, fill = sp)) +
  geom_col(width = 0.5) +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  ylab("Relative abundance") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  xlab("Species") +
  theme_meta() +
  theme(legend.position = "none")
plot(p4)

# phenology





