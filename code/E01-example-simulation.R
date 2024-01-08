
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
  geom_point(mapping = aes(x = RGR, y = SLA, colour = sp), alpha = 0.5) +
  geom_density_2d(mapping = aes(x = RGR, y = SLA, colour = sp), alpha = 0.5, n = 100) +
  scale_colour_viridis_d() +
  ylab("Specific Leaf Area (SLA, cm g-1)") +
  xlab("Relative growth rate (g g-1 day-1)") +
  scale_x_continuous(limits = c(-0.1, 1.05)) +
  scale_y_continuous(limits = c(87, 115)) +
  theme_meta() +
  theme(legend.position = "none") +
  theme_transparent()
plot(p1)

# export the figure for further modification
ggsave(filename = "figures-tables/fig1.pdf", p1,
       unit = "cm", width = 10, height = 8, bg = "transparent")

p2 <- 
  ggplot(data = df.traits) +
  geom_point(mapping = aes(x = RGR, y = SLA, colour = sp), alpha = 0.5) +
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
  theme(legend.position = "none") +
  theme_transparent()
plot(p2)

# export the figure for further modification
ggsave(filename = "figures-tables/fig2.pdf", p2,
       unit = "cm", width = 10, height = 8, bg = "transparent")

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
  theme(legend.position = "none") +
  theme_transparent()
plot(p3)

# export the figure for further modification
ggsave(filename = "figures-tables/fig3.pdf", p3,
       unit = "cm", width = 10, height = 8, bg = "transparent")

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
  theme(legend.position = "none") +
  theme_transparent()
plot(p4)

# export the figure for further modification
ggsave(filename = "figures-tables/fig4.pdf", p4,
       unit = "cm", width = 10, height = 8, bg = "transparent")

# phenology: not necessary


# SLA and SANPP

# get traits from the sim_traits() function
df.traits <- 
  sim_traits(sp = sp, com = com,
             t1 = "SLA", t2 = "RGR", 
             mu_t1 = 100, mu_t2 = 3, sd_t1 = 30, sd_t2 = 2.5, r = r,
             Eind_t1 = 1, Eind_t2 = 0.15,
             t2_scale = 100)

# simulate different communities from a log-normal distribution
df.1 <- sim_abund(sp = sp, com = com, N = N, cv_abund = 0.5)

# bind into a data.frame
df.1 <- dplyr::bind_rows(df.1)

# add the evenness value
df.1$cv_abund <- 0.5

# arrange this data.frame
df1 <- 
  df1 |>
  dplyr::arrange(cv_abund, com, sp)

# bind the trait data to df2
df2 <- dplyr::full_join(df1, dplyr::bind_rows(df.traits), by = c("com", "sp"))

# add the phenology to this data.frame
df2$t0f <- t0f
df2$pheno <- "fixed"

# calculate the parameters that we will use later
df2$SANPP <- with(df2,
                  (abund*(exp(RGR*(t0f)))) )

# summarise to the community-level
df.sum <- 
  df2 |>
  dplyr::group_by(cv_abund, pheno, com) |>
  dplyr::summarise(CWM_SLA = sum(abund*SLA),
                   FD_SLA = sum(abund*( (SLA-mean(SLA))^2 ) ),
                   SANPP = (log(sum(SANPP))/dt), .groups = "drop")
head(df.sum)

# plot
p5 <- 
  ggplot(data = df.sum) +
  geom_point(mapping = aes(x = CWM_SLA, y = SANPP), alpha = 0.75, size = 2.5) +
  scale_colour_viridis_d() +
  xlab("CWM SLA (cm g-1)") +
  ylab("SANPP (g kg-1 day-1)") +
  theme_meta() +
  theme(legend.position = "none") +
  theme_transparent()
plot(p5)

# export the figure for further modification
ggsave(filename = "figures-tables/fig5.pdf", p5,
       unit = "cm", width = 10, height = 8, bg = "transparent")





