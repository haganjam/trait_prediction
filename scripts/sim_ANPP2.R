
# sim_ANPP2()

# load relevant libraries
library(here)

# load relevant functions
source(here("scripts/sim_abun_func.R"))
source(here("scripts/sim_traits_func.R"))
source(here("scripts/sim_pheno_func.R"))

# run a set of simulations
expand.grid(rep = 1:5,
            sp = c(4, 8, 12),
            com = c(12),
            even_par = c(0.1, 0.5, 1),
            lambda = c(50, 100),
            rt1 = seq(0.15, 0.95, 0.15),
            rt2 = seq(0.15, 0.95, 0.15)
            )

# simulate a relative abundance distribution
df.abun <- 
  sim_abund(sp = 10, com = 20, com_mat = FALSE,
            even_par = 0.1, even_mix = FALSE, 
            lambda = 50, lambda_equal = TRUE)

df.traits1 <- 
  sim_traits(sp = 10, com = 20,
             t1 = "SLA", t2 = "RGR", 
             mu_t1 = 100, mu_t2 = 0.05, sd_t1 = 20, sd_t2 = 0.02, r = 0.5,
             Eind_t1 = 5, Eind_t2 = 0.02)

# simulate height and individual level biomass
df.traits2 <- 
  sim_traits(sp = 10, com = 20,
             t1 = "H", t2 = "BN", 
             mu_t1 = 150, mu_t2 = 50, sd_t1 = 50, sd_t2 = 25, r = 0.5,
             Eind_t1 = 1.5, Eind_t2  = )

# simulate phenology
df.pheno <- 
  sim_pheno(com = 20, sp = 10, pheno_mu = 60, pheno_sd = 0,
            Eind_pheno = 0)


# calculate ANPP

# join the traits
df.traits12 <- full_join(df.traits1, df.traits2, by = c("com", "sp") )

# join these data.frames
df.full <- 
  full_join(df.abun, 
            full_join(df.traits12, df.pheno, by = c("com", "sp") ),
            by = c("com", "sp")
  )

# add a column for the standing biomass
df.full$M0 <- df.full$N*df.full$BN

# calculate ANPP
df.full$ANPP<- with(df.full,
                    (M0*(exp(RGR*(t0f)) - 1)) )

# summarise to the community-level
df.sum <- 
  df.full %>%
  group_by(com) %>%
  summarise(sum_M0 = sum(M0),
            CWM_pheno = sum( pi*(t0f) ),
            CWM_SLA = sum(pi*SLA),
            CWM_RGR = sum(pi*RGR),
            CWM_H = sum(pi*H),
            ANPP = (sum(ANPP)/dt), .groups = "drop" )

# fit a linear model to get the RGR r2 value
lm.x <- lm(log(ANPP) ~ CWM_RGR, data = df.sum)
lm.x <- summary(lm.x)
r2_ln_RGR <- lm.x$r.squared

# fit a linear model to get the SLA r2 value on log(ANPP)
lm.x <- lm(log(ANPP) ~ CWM_SLA, data = df.sum)
lm.x <- summary(lm.x)
r2_ln_SLA <- lm.x$r.squared

# fit a linear model to get the H r2 value
lm.x <- lm(log(ANPP) ~ CWM_H, data = df.sum)
lm.x <- summary(lm.x)
r2_ln_H <- lm.x$r.squared

# fit a linear model to get the H r2 value
lm.x <- lm(log(ANPP) ~ CWM_SLA + CWM_H + CWM_SLA:CWM_H, data = df.sum)
lm.x <- summary(lm.x)
r2_ln_SLA_H <- lm.x$r.squared

# pull this into a data.frame
df.out <- data.frame(ANPP_mean = mean(df.sum$ANPP),
                     ANPP_sd = sd(df.sum$ANPP),
                     pheno_mean = mean(df.sum$CWM_pheno),
                     pheno_sd = sd(df.sum$CWM_pheno),
                     M0_negative = any(df.full$M0 < 0),
                     M0_mean = mean(df.sum$sum_M0),
                     M0_sd = sd(df.sum$sum_M0),
                     r2_ln_SLA = round(r2_ln_SLA, 3),
                     r2_ln_RGR = round(r2_ln_RGR, 3),
                     r2_ln_H = round(r2_ln_H, 3),
                     r2_ln_SLA_H = round(r2_ln_SLA_H, 3))






