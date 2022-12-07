
# load relevant libraries
library(dplyr)
library(ggplot2)
library(truncnorm)

# set the number of species and communities (sp = n communities)
sp_com <- 3

# generate three communities where a different species dominates each one
abun_df <- 
  
  lapply(1:sp_com, function(x) {
    
    y <- rep(0, sp_com)
    y[x] <- 1
    
    df <- data.frame(scenario = 1,
                     com = x,
                     sp = paste0("sp_", 1:sp_com),
                     RA = y)
    
    return(df)
  
}) 
abun_df <- bind_rows(abun_df)
print(abun_df)

# generate a data.frame of simulation parameters
sim.df <- expand.grid(rep = 1:10,
                      r = c(0.99, seq(0.1, 0.9, 0.1)),
                      SD_pheno_par = seq(0, 20, 2),
                      SD_M0_par = seq(0, 400, 20)
                      )
dim(sim.df)

# add an id column
sim.df <- bind_cols(data.frame(id = 1:nrow(sim.df)), sim.df)

sim.out <- vector("list", length = nrow(sim.df)) 
for(i in 1:nrow(sim.df)) {
  
  df.x <- sim.df[i,]
  
  # draw RGR and SLA values and trait values
  t12_mu <- 
    faux::rnorm_multi(n = sp_com, 
                      mu = c(100, 0.05),
                      sd = c(20, 0.02),
                      r = df.x[["r"]], 
                      varnames = c("SLA", "RGR"),
                      empirical = TRUE
    ) %>%
    arrange(RGR)
  
  # add a species column
  t12_mu <- bind_cols(data.frame(sp = paste0("sp_", 1:sp_com), t12_mu))
  
  # get the phenology
  pheno <- rtruncnorm(n = sp_com, a = 0, b = Inf, mean = 60, sd = df.x[["SD_pheno_par"]])
  
  # bind the phenology into a data.frame
  pheno_df <- data.frame(sp = paste0("sp_", 1:sp), t0f = pheno)
  
  # bind the phenology to the trait data
  trait_df <- full_join(t12_mu, pheno_df, by = "sp")
  
  # get the total standing biomass of each community
  M0 <- rtruncnorm(n = sp_com, a = 0, b = Inf, mean = 1000, sd = df.x[["SD_M0_par"]])

  # bind the standing biomass into a data.frame
  M0_df <- data.frame(com = 1:sp, M0 = M0)
  
  # bind these data.frames into one large data.frame
  scen.x <- 
    full_join(full_join(abun_df, M0_df, by = "com"),
              trait_df, 
              by = "sp") %>%
    mutate(M0 = M0*RA)
  
  # calculate ANPP
  scen.x$ANPP<- with(scen.x, (M0*(exp(RGR*(t0f)) - 1)) )
  
  # summarise to the community-level
  scen.x <- 
    scen.x %>%
    group_by(com) %>%
    summarise(CWM_SLA = sum(RA*SLA),
              CWM_RGR = sum(RA*RGR),
              ANPP = (sum(ANPP)/90), .groups = "drop" )
  
  # fit a linear model to get the RGR r2 value
  lm.x1 <- lm(log10(ANPP) ~ CWM_SLA, data = scen.x)
  lm.x1 <- summary(lm.x1)
  
  # fit a linear model to get the RGR r2 value
  lm.x2 <- lm((ANPP) ~ CWM_SLA, data = scen.x)
  lm.x2 <- summary(lm.x2)
  
  # pull this into a data.frame
  df.out <- data.frame(ANPP_mean = mean((scen.x$ANPP) ),
                       ANPP_sd = sd((scen.x$ANPP)),
                       pheno_mean = mean(pheno),
                       pheno_range = diff(range(pheno)),
                       M0_mean = mean(M0),
                       M0_range = diff(range(M0)),
                       r2_SLA_log10ANPP = lm.x1$r.squared,
                       r2_SLA_ANPP = lm.x2$r.squared)
  
  # add results to a list
  sim.out[[i]] <- df.out
  
}

# bind these simulations results into a data.frame
sim.out <- bind_rows(sim.out)

# add simulation parameters
sim.out <- bind_cols(sim.df, sim.out)

# generate the three proof of concept graphs

# perfect case
sim.out %>%
  filter(r == 0.99, SD_pheno_par == 0, SD_M0_par == 0) %>%
  ggplot(data = .,
         mapping = aes(x = r2_SLA_ANPP)) +
  geom_density() +
  theme_classic()

# r = 0.99, equal phenology but variation in standing biomass
sim.out %>%
  filter(r == 0.99, SD_pheno_par == 0) %>%
  ggplot(data = .,
         mapping = aes(x = M0_range, y = r2_SLA_ANPP)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE) +
  theme_classic()

# r = 0.99, equal phenology but variation in standing biomass
sim.out %>%
  filter(r == 0.99, SD_M0_par == 0) %>%
  ggplot(data = .,
         mapping = aes(x = pheno_range, y = r2_SLA_ANPP)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE) + 
  theme_classic()
    
# equal M0, equal phenology but variation in r
sim.out %>%
  filter(r != 0.99, SD_M0_par == 0, SD_pheno_par == 0) %>%
  ggplot(data = .,
         mapping = aes(x = r, y = r2_SLA_ANPP)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE) + 
  theme_classic()

### END

