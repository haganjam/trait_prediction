
# Simulating predictive value of traits

# How to simulate random communities?

# What is the theoretical limit to the correlation between traits and ecosystem processes

# load relevant libraries
library(mobsim)
library(vegan)
library(betapart)

# function to simulate a regional species pool: Sim_regional

# args

#' @param nspec - number of species in the regional species pool
#' @param nind - number of individuals in the regional species pool (should be a very large number like 100000) 
#' @param sad - type of species abundance distribution (default lognormal but see mobsim for details)
#' @param plot - T or F whether to plot the histogram of abundance distributions

Sim_regional <- function(nspec = 10, nind = 100000, plot = TRUE) {
  
  # simulate a species pool
  reg_vector <- sim_sad(s_pool = nspec, n_sim = nind, sad_type = "lnorm")
  attr(reg_vector, "class") <- NULL
  names(reg_vector) <- NULL
  
  # add species names to these species
  specnames <- paste("sp", 1:length(reg_vector), sep = "")
  names(reg_vector) <- specnames
  
  # plot a histogram of the regional species pool abundance distribution
  if (plot) {
    hist(reg_vector)
  }
  
  # add a species number attribute to the reg_vector
  attr(reg_vector, "specnumber") <- nspec
  
  return(reg_vector)
  
}

# function to sample from a regional pool vector from Sim_regional: Reg_sampler

# args

#' @param reg_vector - regional species pool vector from Reg_sampler 
#' @param Kmean - average carrying capacity i.e. number of individuals that a patch can support (varies with a Poisson distribution)

Reg_sampler <- function(reg_vector,  Kmean = 100) {
  
  # get the specnames variable from the regional pool
  specnames <- attr(reg_vector, "specnumber")
  
  # make a list of specnames to serve as the base
  specnames <- paste("sp", 1:specnames, sep = "")
  
  # draw a sample from the distribution
  x <- sample(rep(names(reg_vector), reg_vector), size = rpois(n = 1, Kmean), replace = TRUE )
  x <- table(x)
  
  # extract the missing names and get an ordered list of missing names and present names
  missingnames <- specnames[!specnames %in% names(x)]
  fullnames <- c(names(x), missingnames)
  
  # add zeros for the missing names and add the missing names
  names(x) <- NULL
  x <- c(x, rep(0, length(missingnames) ))
  names(x) <- fullnames
  x <- x[sort(names(x))]
  
  # return the sample
  return(x)
  
}

# function to simulate n independent random samples from a regional species pool: Sim_comms

# args

#' @param reg_vector - regional species pool vector from Reg_sampler 
#' @param Kmean - average carrying capacity i.e. number of individuals that a patch can support (varies with a Poisson distribution)
#' @param ncomm - number of independent communities to sample 
#' @param dom_variation - TRUE of FALSE: highest abundance species identity varies across communities

Sim_comms <- function(reg_vector, Kmean = 100, ncomm = 5, dom_variation = FALSE) {
  
  y <- 
    lapply(1:ncomm, function(vec) {
      x <- Reg_sampler(reg_vector = reg_vector, Kmean = Kmean)
      return(x)
    })
  
  if (dom_variation) {
    y <- 
      lapply(y, function(z) {
        maxin <- max(z)
        maxid <- which(z == maxin)
        spid <- 1:length(z)
        sid <- sample(spid[!(spid %in% maxid)], 1)
        z[maxid] <- z[sid]
        z[sid] <- maxin
        return(z)
      })
  }
  
  return(y)
  
}

# set the number of species to simulate
nspec <- 30
ncomm <- 1

# simulate a set of species traits
traits.m <- round(rnorm(n = nspec, mean = 10, 3), 1)
names(traits) <- paste("sp", 1:nspec, sep = "")

# simulate a regional pool of species
reg_vector <- Sim_regional(nspec = nspec, nind = 100000, plot = TRUE)

# test the sim_comms function
Sim_comms(reg_vector = reg_vector, Kmean = 100, ncomm = ncomm , dom_variation = TRUE)






# maybe we need I need a bit more reading

# 1. simple model showing that the way we sample can strongly affect our quantification of CWM and FD
# 2. more complex model of using traits to estimate, e.g. photosynthesis, and then biomass

# how to do this... 

# maybe we can also check how intraspecific variation affects different aspects of the trait distribution
# e.g.: CWM, FD, kurtosis and skew

# different types of sampling?

# intraspecific traits may not be a random sample from the regional pool?

# are intraspecific traits within communities a random sample?
# are they selected to maximise variation within or between species?

# load relevant libraries
library(mobsim)
library(dplyr)

# how many species in the species pool?
sp_n <- 10

# get sp_n mean values for the trait distributions of each species
sp_u <- rnorm(n = sp_n, 100, 30)

# get sp_n sd values for the trait distributions of each species
sp_sd <- runif(sp_n, 0, 10)

# simulate a local abundance distrbution
reps <- 10
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

x


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



  







