
# sim_ANPP2()

# load relevant libraries
library(here)

# load relevant functions
source(here("scripts/A01_sim_abun_func.R"))
source(here("scripts/A02_sim_pheno_func.R"))
source(here("scripts/A03_sim_traits_func.R"))

# load relevant libraries
library(dplyr)
library(ggplot2)

# set up number of species
sp <- 10

# set up the number of communities
com <- 20

# get traits from the sim_traits() function
df.traits1 <- 
  sim_traits(sp = sp, com = com,
             t1 = "SLA", t2 = "RGR", 
             mu_t1 = 100, mu_t2 = 10, sd_t1 = 30, sd_t2 = 2, r = 0.75,
             Eind_t1 = 5, Eind_t2 = 0.5)

# plot this trait distributions
ggplot(data = df.traits1,
       mapping = aes(x = SLA, fill = sp, colour = sp)) +
  geom_density(alpha = 0.75) +
  scale_fill_viridis_d(option = "D") +
  scale_colour_viridis_d(option = "D") +
  theme_classic()

ggplot(data = df.traits1,
       mapping = aes(x = RGR, fill = sp, colour = sp)) +
  geom_density(alpha = 0.75) +
  scale_fill_viridis_d(option = "D") +
  scale_colour_viridis_d(option = "D") +
  theme_classic()

ggplot(data = df.traits1,
       mapping = aes(x = SLA, y = RGR, colour = sp)) +
  geom_point() +
  scale_colour_viridis_d(option = "D") +
  theme_classic()

# calculate the average trait values of the different species across communities
df.traits1

# simulate a set of communities with some relative abundances
sad_list <- 
  
  lapply(1:com, function(x) {
  
  # sp_pool size determined by J and theta together (can we calculate this?)
  sad.x <- mobsim::sim_sad(s_pool = NA, n_sim = 100, sad_type = c("mzsm"),
                       sad_coef = list("J" = 10000, "theta" = 1))
  
  # remove the class attributes
  attr(sad.x, "class") <- NULL
  names(sad.x) <- NULL
  
  # initialise an empty vector with 10 species
  sp_vec <- rep(0, sp)
  sp_vec[sample(1:sp, size = length(sad.x))] <- sad.x
  
  # wrap into a data.frame
  sp_df <- tibble(com = as.character(x),
                  sp = paste0("sp_", 1:sp),
                  abund = sp_vec)
  
  return(sp_df)
  
} )

# bind into a data.frame
sad_df <- bind_rows(sad_list)

# bind the sad_df data.frame to the trait data
df.test <- full_join(sad_df, df.traits1, by = c("com", "sp"))
print(df.test)

# check the correlation between RGR and abundance
plot(df.test$RGR, df.test$abund)
ggplot(data = df.test, 
       mapping = aes(x = RGR, y = abund, colour = com)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_viridis_d() +
  theme_classic() +
  theme(legend.position = "none")

# reorder to make the rank correlation perfect and re-plot
com1 <- df.test[df.test$com == 4,]
rnk <- rank(com1$RGR, ties.method = "random")
round(sort(com1$abund), 1)[rnk]

plot(com1$RGR, round(sort(com1$abund), 1)[rnk])


# run a set of simulations
sim.df <- expand.grid(rep = 1:2,
                      sp = c(5, 10),
                      com = c(20),
                      even_par = c(0.1, 1),
                      lambda = c(100),
                      lambda_equal = c(FALSE),
                      rt1 = c(0.1, 0.5, 0.9),
                      rt2 = c(0.1, 0.5, 0.9),
                      trait_pheno_sd = c(0, 10),
                      dt = 90 )

# add an ID column
sim.df <- bind_cols(data.frame(id = 1:nrow(sim.df)), sim.df)
dim(sim.df)

# simulate a relative abundance distribution
sim.out <- vector("list", length = nrow(sim.df)) 
for(i in 1:nrow(sim.df)) {

  df.x <- sim.df[i,]
  
  df.abun <- 
    sim_abund(sp = df.x[["sp"]], com = df.x[["com"]], com_mat = FALSE,
              even_par = df.x[["even_par"]], even_mix = FALSE, 
              lambda = df.x[["lambda"]], lambda_equal = df.x[["lambda_equal"]])
  
  df.traits1 <- 
    sim_traits(sp = df.x[["sp"]], com = df.x[["com"]],
               t1 = "SLA", t2 = "RGR", 
               mu_t1 = 100, mu_t2 = 0.05, sd_t1 = 30, sd_t2 = 0.02, r = df.x[["rt1"]],
               Eind_t1 = 0, Eind_t2 = 0)
  
  # simulate height and individual level biomass
  df.traits2 <- 
    sim_traits(sp = df.x[["sp"]], com = df.x[["com"]],
               t1 = "FD", t2 = "t0f", 
               mu_t1 = 10, mu_t2 = 60, sd_t1 = 3, sd_t2 = df.x[["pheno_sd"]], r = df.x[["rt2"]],
               Eind_t1 = 0, Eind_t2 = 0)
  
  # calculate ANPP
  
  # join the traits
  df.traits12 <- full_join(df.traits1, df.traits2, by = c("com", "sp") )
  
  # join these data.frames
  df.full <- 
    full_join(df.abun, 
              df.traits12,
              by = c("com", "sp")
    )
  
  # calculate ANPP
  df.full$ANPP<- with(df.full,
                      (N*(exp(RGR*(t0f)) - 1)) )
  
  # summarise to the community-level
  df.sum <- 
    df.full %>%
    group_by(com) %>%
    summarise(sum_N = sum(N),
              CWM_pheno = sum( pi*(t0f) ),
              CWM_SLA = sum(pi*SLA),
              CWM_RGR = sum(pi*RGR),
              CWM_FD = sum(pi*FD),
              ANPP = (sum(ANPP)/df.x[["dt"]]), .groups = "drop" )
  
  # fit a linear model to get the RGR r2 value
  lm.x <- lm(log(ANPP) ~ CWM_RGR, data = df.sum)
  lm.x <- summary(lm.x)
  r2_ln_RGR <- lm.x$r.squared
  
  # fit a linear model to get the SLA r2 value on log(ANPP)
  lm.x <- lm(log(ANPP) ~ CWM_SLA, data = df.sum)
  lm.x <- summary(lm.x)
  r2_ln_SLA <- lm.x$r.squared
  
  # fit a linear model to get the H r2 value
  lm.x <- lm(log(ANPP) ~ CWM_FD, data = df.sum)
  lm.x <- summary(lm.x)
  r2_ln_FD <- lm.x$r.squared
  
  # fit a linear model to get the H r2 value
  lm.x <- lm(log(ANPP) ~ CWM_SLA + CWM_FD + CWM_SLA:CWM_FD, data = df.sum)
  lm.x <- summary(lm.x)
  r2_ln_SLA_FD <- lm.x$r.squared
  
  # pull this into a data.frame
  df.out <- data.frame(ANPP_mean = mean((df.sum$ANPP) ),
                       ANPP_sd = sd((df.sum$ANPP)),
                       pheno_mean = mean(df.sum$CWM_pheno),
                       pheno_sd = sd(df.sum$CWM_pheno),
                       N_negative = any(df.full$N < 0),
                       N_mean = mean(df.sum$sum_N),
                       N_sd = sd(df.sum$sum_N),
                       r2_ln_SLA = round(r2_ln_SLA, 3),
                       r2_ln_RGR = round(r2_ln_RGR, 3),
                       r2_ln_FD = round(r2_ln_FD, 3),
                       r2_ln_SLA_FD = round(r2_ln_SLA_FD, 3))
  
  # write into a list
  sim.out[[i]] <- df.out
  
}

# bind into a data.frame
sim.out <- bind_rows(sim.out)
nrow(sim.out)

# check the summary statistics
summary(sim.out)

# bind these data to simulation parameters
sim.out <- bind_cols(sim.df, sim.out)

# let's see what it looks like to examine r2_ln_SLA and r2_ln_H
library(ggplot2)
sim.out %>%
  mutate(rt1 = as.character(rt1)) %>%
  filter(trait_pheno_sd == 0) %>%
  ggplot(data = .,
         mapping = aes(x = N_sd, y = r2_ln_SLA, colour = rt1)) +
  geom_point(size = 1.8) +
  geom_smooth(method = "lm") +
  theme_classic()








