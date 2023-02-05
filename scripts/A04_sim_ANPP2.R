
# sim_ANPP2()

# load relevant libraries
library(here)

# load relevant functions
source(here("scripts/A01_sim_abun_func.R"))
source(here("scripts/A03_sim_traits_func.R"))

# load relevant libraries
library(dplyr)
library(ggplot2)

# set the r between RGR and SLA
r <- 0.99

# set the phenology
t0f <- 10

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
                                  abund_factor = seq(1, 5, length.out = com)))

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




