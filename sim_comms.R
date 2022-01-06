
# Simulating predictive value of traits

# load relevant libraries
library(mobsim)
library(dplyr)

# how many species in the species pool?
sp_n <- 10

# get sp_n mean values for the trait distributions of each species
sp_u <- rnorm(n = sp_n, 100, 30)
sp_u

# get sp_n sd values for the trait distributions of each species
sp_sd <- runif(sp_n, 0, 10)
sp_sd

# simulate a local abundance distrbution
reps <- 2
sim_divs <- rep(2:sp_n, each = reps)

x <- vector("list", length = length(sim_divs))
for (i in 1:length(sim_divs) ) {
  
  y <- sim_sad(s_pool = sim_divs[i],
               n_sim = round(rnorm(n = 1, mean = 100, sd = 10), 0) ,
               sad_type = "lnorm",
               fix_s_sim = FALSE,
               drop_zeros = TRUE)
  
  attr(y, "class") <- NULL
  
  x[[i]] <- y

}

y <- 
  lapply(x, function(a) {
    
    n <- 1:length(a)
    sp_names <- names(a)
    
    tvals <- vector("list", length = length(n))
    for (i in n) {
      
      v <- rnorm(n = a[i], mean = sp_u[i], sd = sp_sd[i])
      u <- rep(sp_u[i], a[i] )
      df <- data.frame(species = sp_names[i],
                       trait_val = v,
                       trait_u = u)
      tvals[[i]] <- df
      
    }
    
    bind_rows(tvals)
    
  } )

# input list
y

# imperfect abundance detection function
imperfect_abundance <- function(df_in, prop_sample = 0.5) {
  
  s_probs <- runif(nrow(df_in), 0, 1)
  s_probs <- s_probs/sum(s_probs)
  
  s_rows <- sample(1:nrow(df_in), size = round(prop_sample*nrow(df_in), 0), prob = s_probs)
  
  df_in[s_rows, ]
  
}

# does the function work?
imperfect_abundance(df_in = y[[18]])


# CWM summary function
CWM_summary <- function(df_in, trait = "intra") {
  
  if (trait == "intra" ) {
    
    df_in %>%
      group_by(species) %>%
      summarise("intra" = mean(trait_val, na.rm = TRUE))
    
  } else if (trait == "true_inter" ) {
    
    df_in %>%
      group_by(species) %>%
      summarise("true_inter" = mean(trait_u, na.rm = TRUE))
    
  } else if (trait == "other_inter") {
    
    df_in %>%
      group_by(species) %>%
      summarise("other_inter" = rnorm(1, mean(trait_u, na.rm = TRUE), 10) )
    
  }
  
}

# test the CWM summary function
CWM_summary(df_in = y[[18]], trait = "intra")

df <- 
  lapply(y, function(x) {
  
 imperfect_abundance(df_in = x, prop_sample = 0.5)
  
} )


# calculate true community weighted mean and functional diversity
# assuming perfect information
df.sim <- bind_rows(y, .id = "community")

# model biomass as an additive linear combination of these quantities
b1 = 1.3

df.bio <- 
  df.sim %>%
  group_by(community) %>%
  summarise(CWM = mean(trait_val, na.rm = TRUE) ) %>%
  mutate(biomass = b1*CWM)

# we see that there is a perfect correlation
plot(df.bio$CWM, df.bio$biomass)
cor(df.bio$CWM, df.bio$biomass)


# assume different levels of error

# 1. imperfect detection of species abundances and trait means
df.im <- 
  bind_rows(df, .id = "community") %>%
  group_by(community) %>%
  summarise(CWM = mean(trait_u, na.rm = TRUE) )

plot(df.im$CWM, df.bio$biomass)  
cor(df.im$CWM, df.bio$biomass) 

# single mean trait estimate

# error in biomass estimation



  







