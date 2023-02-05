
# sim_ANPP2()

# load relevant libraries
library(here)

# load relevant functions
source(here("scripts/A01_sim_abun_func.R"))
source(here("scripts/A03_sim_traits_func.R"))

# load relevant libraries
library(dplyr)
library(purrr)
library(ggplot2)

# set the r between RGR and SLA
r <- 0.75

# set the phenology
t0f <- 20

# set the timeframe
dt <- 30

# set up number of species
sp <- 10

# set up the number of communities
com <- 20

# total individuals
N <- 100

# evenness parameters
even <- c(1, 20)

# get traits from the sim_traits() function
df.traits <- 
  sim_traits(sp = sp, com = com,
             t1 = "SLA", t2 = "RGR", 
             mu_t1 = 100, mu_t2 = 10, sd_t1 = 30, sd_t2 = 2, r = r,
             Eind_t1 = 5, Eind_t2 = 0.5,
             t2_scale = 100)

# loop over the two different evenness parameters
df3 <- vector("list", length = length(even))
for(i in 1:length(even)) {
  
  # simulate different communities from a log-normal distribution
  df.abund <- sim_abund(sp = sp, com = com, N = N, cv_abund = even[i])
  
  # bind df.abund to the traits and then reorder or not
  
  # abundance and RGR are not correlated
  df1 <- 
    mapply(function(x, y){
      
      df <- full_join(x, y, by = c("com", "sp"))
      df$trait_abund_cor <- FALSE
      return(df)
      
    }, df.abund, df.traits, SIMPLIFY = FALSE)
  
  # abundance and RGR are correlated
  df2 <- 
    mapply(function(x, y){
      
      df <- full_join(x, y, by = c("com", "sp"))
      rnk <- rank(df$RGR, ties.method = "random")
      df$abund <- round(sort(df$abund), 1)[rnk]
      df$trait_abund_cor <- TRUE
      return(df)
      
    }, df.abund, df.traits, SIMPLIFY = FALSE)
  
  # bind this together
  df12 <- bind_rows( bind_rows(df1), bind_rows(df2))
  
  # add the evenness value
  df12$cv_abund <- even[i]
  
  # add the df12 to the df.even data.frame
  df3[[i]] <- df12
  
}

# bind the df.even together
df3 <- bind_rows(df3)

# arrange this data.frame
df3 <- 
  df3 %>%
  arrange(cv_abund, trait_abund_cor, com, sp)

# add an abundance variation column to df3
df3a <- 
  df3 %>%
  mutate(abund_var = FALSE)

# add variation in overall abundance between communities
abund_var <- as_tibble(data.frame(com = c(1:com),
                                  abund_factor = seq(1, 10, length.out = com)))

df3b <- 
  full_join(df3, abund_var, by = "com") %>%
  mutate(abund = round(abund*abund_factor, 0 )) %>%
  select(-abund_factor) %>%
  mutate(abund_var = TRUE)

# join these data.frame together
df4 <- bind_rows(df3a, df3b)

# reorder the columns and arrange
df4 <- 
  df4 %>%
  select(abund_var, cv_abund, trait_abund_cor, com, sp, abund, SLA, RGR) %>%
  arrange(abund_var, cv_abund, trait_abund_cor, com, sp)

# calculate the parameters that we will use later
df4$ANPP <- with(df4,
                 (abund*(exp(RGR*(t0f)) - 1)) )

# summarise to the community-level
df.sum <- 
  df4 %>%
  group_by(abund_var, cv_abund, trait_abund_cor, com) %>%
  mutate(pi = abund/sum(abund)) %>%
  summarise(abund_tot = sum(abund),
            CWM_SLA = sum(pi*SLA),
            CWM_RGR = sum(pi*RGR),
            ANPP = (sum(ANPP)/dt), .groups = "drop" )

# set-up a vector of groups to split by
groups <- with(df.sum, as.integer(factor(paste0(abund_var, cv_abund, trait_abund_cor))) )

# split by group and then analyse using the linear models
df.sum <- split(df.sum, groups)

# lapply over the different samples and calculate r2 values
df.out <- 
  
  lapply(df.sum, function(df.x) {
  
  # fit a linear model to get the SLA r2 value on log(ANPP)
  lm.x <- lm(log(ANPP) ~ CWM_SLA, data = df.x)
  lm.x <- summary(lm.x)
  r2_ln_ANPP <- lm.x$r.squared
  
  # fit a linear model to get the SLA r2 value on log(ANPP)
  lm.x <- lm((ANPP) ~ CWM_SLA, data = df.x)
  lm.x <- summary(lm.x)
  r2_ANPP <- lm.x$r.squared
  
  df.out <- tibble(abund_var = df.x$abund_var[1],
                   cv_abund = df.x$cv_abund[1], 
                   trait_abund_cor = df.x$trait_abund_cor[1],
                   r_RGR_SLA = r,
                   r2_ANPP = r2_ANPP,
                   r2_ln_ANPP = r2_ln_ANPP)
  
} )

# bind this into a data.frame
df.out <- bind_rows(df.out)
print(df.out)

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








