
# sim_ANPP2()

# load relevant libraries
library(here)

# load relevant functions
source(here("scripts/A01_sim_abun_func.R"))
source(here("scripts/A03_sim_traits_func.R"))

# load relevant libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# fixed parameters

# abundance variation factor
abund_factor <- 3

# set the phenology
t0f <- 30

# set the timeframe
dt <- 50

# set up number of species
sp <- 10

# set up the number of communities
com <- 20

# total individuals
N <- 100

# evenness parameters
even <- c(5, 20)

# expand grid to get 10 replicates of each correlation value
r_vec <- c(0.01, seq(0.05, 0.95, 0.05), 0.99)
par.grid <- expand.grid(rep = 1:10, r = r_vec)

sim.df <- vector("list", length = nrow(par.grid))
for(j in 1:nrow(par.grid)) {
  
  # extract parameters from the par.grid
  r <- par.grid[j,"r"]
  rep_N <- par.grid[j, "rep"]
  
  # get traits from the sim_traits() function
  df.traits <- 
    sim_traits(sp = sp, com = com,
               t1 = "SLA", t2 = "RGR", 
               mu_t1 = 100, mu_t2 = 10, sd_t1 = 30, sd_t2 = 2.5, r = r,
               Eind_t1 = 1, Eind_t2 = 0.1,
               t2_scale = 100)
  
  # loop over the two different evenness parameters
  df1 <- vector("list", length = length(even))
  for(i in 1:length(even)) {
    
    # simulate different communities from a log-normal distribution
    df.abund <- sim_abund(sp = sp, com = com, N = N, cv_abund = even[i])
    
    # bind into a data.frame
    df.abund <- bind_rows(df.abund)
    
    # add the evenness value
    df.abund$cv_abund <- even[i]
    
    # add the df12 to the df.even data.frame
    df1[[i]] <- df.abund
    
  }
  
  # bind the df.even together
  df1 <- bind_rows(df1)
  
  # arrange this data.frame
  df1 <- 
    df1 %>%
    arrange(cv_abund, com, sp)
  
  # add an abundance variation column to df3
  df1a <- 
    df1 %>%
    mutate(abund_var = FALSE)
  
  # add variation in overall abundance between communities
  abund_var <- as_tibble(data.frame(com = c(1:com),
                                    abund_factor = seq(1, abund_factor, length.out = com)))
  
  df1b <- 
    full_join(df1, abund_var, by = "com") %>%
    mutate(abund = round(abund*abund_factor, 0 )) %>%
    select(-abund_factor) %>%
    mutate(abund_var = TRUE)
  
  # join these data.frame together
  df2 <- bind_rows(df1a, df1b)
  
  # bind the trait data to df2
  df3 <- full_join(df2, bind_rows(df.traits), by = c("com", "sp"))
  
  # reorder the columns and arrange
  df3 <- 
    df3 %>%
    select(abund_var, cv_abund, com, sp, abund, SLA, RGR) %>%
    arrange(abund_var, cv_abund, com, sp)
  
  # calculate the parameters that we will use later
  df3$ANPP <- with(df3,
                   (abund*(exp(RGR*(t0f)) - 1)) )
  
  # summarise to the community-level
  df.sum <- 
    df3 %>%
    group_by(abund_var, cv_abund, com) %>%
    mutate(pi = abund/sum(abund)) %>%
    summarise(abund_tot = sum(abund),
              CWM_SLA = sum(pi*SLA),
              FD_SLA = sum(pi*( (SLA-mean(SLA))^2 ) ),
              ANPP = (sum(ANPP)/dt), .groups = "drop" )
  
  # set-up a vector of groups to split by
  groups <- with(df.sum, as.integer(factor(paste0(abund_var, cv_abund))) )
  
  # split by group and then analyse using the linear models
  df.sum.list <- split(df.sum, groups)
  
  # lapply over the different samples and calculate r2 values
  df.out <- 
    
    lapply(df.sum.list, function(df.x) {
      
      # fit a linear model to get the abun r2 on ANPP
      lm.x <- lm(log(ANPP) ~ abund_tot, data = df.x)
      lm.x <- summary(lm.x)
      r2_abund <- lm.x$r.squared
      
      # fit a linear model to get the SLA r2 value on log(ANPP)
      lm.x <- lm(log(ANPP) ~ CWM_SLA, data = df.x)
      lm.x <- summary(lm.x)
      r2_CWM_SLA <- lm.x$r.squared
      
      # fit a linear model to get the SLA r2 value on log(ANPP)
      lm.x <- lm(log(ANPP) ~ FD_SLA, data = df.x)
      lm.x <- summary(lm.x)
      r2_FD_SLA <- lm.x$r.squared
      
      # fit a linear model to get the SLA r2 value on log(ANPP)
      lm.x <- lm(log(ANPP) ~ CWM_SLA*FD_SLA, data = df.x)
      lm.x <- summary(lm.x)
      r2_CWM_FD_SLA <- lm.x$r.squared
      
      df.out <- tibble(abund_var = df.x$abund_var[1],
                       cv_abund = df.x$cv_abund[1],
                       r_RGR_SLA = r,
                       r2_abund,
                       r2_CWM_SLA = r2_CWM_SLA,
                       r2_FD_SLA = r2_FD_SLA,
                       r2_CWM_FD_SLA = r2_CWM_FD_SLA)
      
    } )
  
  # bind this into a data.frame
  df.out <- bind_rows(df.out)
  
  # add metadata from the par.grid
  df.out$rep <- rep_N
  
  sim.df[[j]] <- df.out
  
}

# bind into a data.frame
sim.df <- bind_rows(sim.df)
head(sim.df)

# what is the highest correlation we find?
max_r2 <- 
  sim.df %>%
  filter(r_RGR_SLA == 0.99) %>%
  group_by(abund_var, cv_abund, r_RGR_SLA) %>%
  summarise(r2 = mean(r2_CWM_FD_SLA), .groups = "drop") %>%
  filter(r2 == max(r2)) %>%
  pull(r2)

# plot these results
ggplot(data = sim.df, 
       mapping = aes(x = r_RGR_SLA,
                     y = r2_CWM_FD_SLA,
                     colour = abund_var)) +
  geom_jitter(width = 0.01, shape = 1, alpha = 0.6) +
  geom_smooth(se = TRUE, alpha = 0.3, size = 0.5, method = "lm",
              formula = y~poly(x, 2), show.legend = FALSE) +
  facet_wrap(~cv_abund) +
  scale_colour_viridis_d(option = "C", end = 0.9) +
  ylab("r2 SLA ~ ANPP") +
  xlab("r RGR ~ SLA") +
  geom_hline(yintercept = max_r2, linetype = "dashed", colour = "red") +
  theme_classic() +
  theme(legend.position = "top")
  
# when there is abundance variation which is better? SLA or just abundance?
sim.df %>%
  filter(abund_var == TRUE) %>%
  select(cv_abund, rep, r_RGR_SLA, r2_abund, r2_CWM_FD_SLA) %>%
  pivot_longer(cols = starts_with("r2"),
               names_to = "abund_SLA",
               values_to = "r2") %>%
  ggplot(data = ., 
         mapping = aes(x = r_RGR_SLA,
                       y = r2,
                       colour = abund_SLA)) +
  geom_jitter(width = 0.01, shape = 1, alpha = 0.6) +
  geom_smooth(se = TRUE, alpha = 0.3, size = 0.5, method = "lm",
              formula = y~poly(x, 2), show.legend = FALSE) +
  facet_wrap(~cv_abund) +
  scale_colour_viridis_d(option = "C", end = 0.9) +
  ylab("r2") +
  xlab("r RGR ~ SLA") +
  theme_classic() +
  theme(legend.position = "top")

  









