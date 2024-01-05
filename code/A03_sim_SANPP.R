
# sim_ANPP2()

# load relevant functions
source("code/A01_sim_abun_func.R")
source("code/A02_sim_traits_func.R")

# load relevant libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# fixed parameters

# set the phenology
t0f <- 30

# set the timeframe
dt <- 50

# set up number of species
sp <- 5

# set up the number of communities
com <- 10

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
               mu_t1 = 100, mu_t2 = 3, sd_t1 = 30, sd_t2 = 2.5, r = r,
               Eind_t1 = 1, Eind_t2 = 0.15,
               t2_scale = 100)
  
  # loop over the two different evenness parameters
  df1 <- vector("list", length = length(even))
  for(i in 1:length(even)) {
    
    # simulate different communities from a log-normal distribution
    df.abund <- sim_abund(sp = sp, com = com, N = N, cv_abund = even[i])
    
    # bind into a data.frame
    df.abund <- dplyr::bind_rows(df.abund)
    
    # add the evenness value
    df.abund$cv_abund <- even[i]
    
    # add the df12 to the df.even data.frame
    df1[[i]] <- df.abund
    
  }
  
  # bind the df.even together
  df1 <- dplyr::bind_rows(df1)
  
  # arrange this data.frame
  df1 <- 
    df1 |>
    dplyr::arrange(cv_abund, com, sp)
  
  # bind the trait data to df2
  df2 <- dplyr::full_join(df1, dplyr::bind_rows(df.traits), by = c("com", "sp"))
  
  # reorder the columns and arrange
  df2 <- 
    df2 |>
    dplyr::select(cv_abund, com, sp, abund, SLA, RGR) |>
    dplyr::arrange(cv_abund, com, sp)
  
  # calculate the parameters that we will use later
  df2$SANPP <- with(df2,
                   (abund*(exp(RGR*(t0f)))) )
  
  # summarise to the community-level
  df.sum <- 
    df2 |>
    dplyr::group_by(cv_abund, com) |>
    dplyr::summarise(CWM_SLA = sum(abund*SLA),
                     FD_SLA = sum(abund*( (SLA-mean(SLA))^2 ) ),
                     SANPP = (log(sum(SANPP))/dt), .groups = "drop")
  
  # set-up a vector of groups to split by
  groups <- with(df.sum, as.integer(factor(paste0(cv_abund))) )
  
  # split by group and then analyse using the linear models
  df.sum.list <- split(df.sum, groups)
  
  # lapply over the different samples and calculate r2 values
  df.out <- 
    
    lapply(df.sum.list, function(df.x) {
      
      # fit a linear model to get the SLA CWM r2 value on SANPP
      lm.x <- lm( (SANPP) ~ CWM_SLA, data = df.x)
      lm.x <- summary(lm.x)
      r2_CWM_SLA <- lm.x$r.squared
      
      # fit a linear model to get the SLA FD r2 value on SANPP
      lm.x <- lm( (SANPP) ~ FD_SLA, data = df.x)
      lm.x <- summary(lm.x)
      r2_FD_SLA <- lm.x$r.squared
      
      # fit a linear model to get the CWM and SD of SLA r2 on SANPP
      lm.x <- lm( (SANPP) ~ CWM_SLA*FD_SLA, data = df.x)
      lm.x <- summary(lm.x)
      r2_CWM_FD_SLA <- lm.x$r.squared
      
      df.out <- dplyr::tibble(cv_abund = df.x$cv_abund[1],
                              r_RGR_SLA = r,
                              r2_CWM_SLA = r2_CWM_SLA,
                              r2_FD_SLA = r2_FD_SLA,
                              r2_CWM_FD_SLA = r2_CWM_FD_SLA)
      
    } )
  
  # bind this into a data.frame
  df.out <- dplyr::bind_rows(df.out)
  
  # add metadata from the par.grid
  df.out$rep <- rep_N
  
  sim.df[[j]] <- df.out
  
}

# bind into a data.frame
sim.df <- dplyr::bind_rows(sim.df)
head(sim.df)

# plot these results

# no variation in phenology among species
ggplot(data = sim.df, 
       mapping = aes(x = r_RGR_SLA,
                     y = r2_CWM_FD_SLA)) +
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
  
# variation in phenology among species (simulate this)
  








