#'
#' @title: Theoretical limits to the correlation between plant traits and productivity
#' 
#' @description: Here, we use the model of primary productivity based on relative growth
#' rates of individual species and their biomass proportions in the communities to test the
#' theoretical limits of the predictability of primary productivity from plant functional traits.
#' 
#' @details: Garnier et al. (2004, Ecology) proposed a formal, causal model linking species'
#' standing biomass relative growth rates and phenology to primary productivity:
#' 
#' ANPP (g g-1 day-1) = sum(M0*(exp(RGRi*(tf-t0)i) - 1) )/dt
#' 
#' where: 
#' M0 is standing, living biomass in a community of species i
#' RGR is the maximum relative growth rate (g g-1 day-1) determined experimentally
#' tf-t0 is the period of growth for species i
#' dt is the period over which the measurements too place
#' 
#' Thus, we can use this equation to simulate ANPP based on known RGRs, M0 and
#' relative growth rates. 
#' 
#' We then assume as per Garnier et al. (2004) and Reich et al. (1997) that certain
#' easy-to-measure plant traits like SLA, leaf nitrogen etc. can predict RGR and, potentially,
#' M0 (e.g. plant height)
#' 
#' We then can ask if there is uncertainty in that prediction along with phenological
#' uncertainty and uncertainty in the actual SANPP measurements, what the theoretical
#' limits to predicting ANPP from plant traits is.
#' 
#' @function - sim_SP()
#' @param com - number of communities to simulate
#' @param sp - number of species to simulate
#' @param even_par - evenness parameter corresponding to the alpha parameter in the Dirichlet distribution
#' @param even_mix - evenness parameter is the same for all communities (FALSE) or it differs (TRUE)
#' @param com_mat - custom community matrix with sites as rows and species as columns
#' @param mu - vector of mean values for SLA and RGR
#' @param sd - vector of sd values for SLA and RGR
#' @param r -  correlation coefficient (Pearson) between SLA and RGR
#' @param bio_m - average standing biomass of a normal distribution
#' @param bio_sd - standard deviation of a normal distribution describing average standing biomass
#' @param dt - growing season length (days)
#' @param pheno - either picks random phenology or keeps all phenologies equal for all species
#'  1. pheno = "random"
#'  2. pheno = "equal"
#' 

sim_SP <- function(com = 20, sp = 10, even_par = 0.5, even_mix = FALSE, com_mat = NA, 
                   mu = c(100, 0.5), sd = c(40, 0.2), r = 0.5,
                   bio_m = 200, bio_sd = 50,
                   dt = 60, pheno = "random" ) {
  
  # load relevant libraries
  library(dplyr)
  library(betapart)
  library(faux)
  
  # simulate a relative abundance distribution
  if (is.matrix(com_mat)) {
    
    RA <- com_mat
    
  } else {
    
    # get the alpha values randomly if even_mix == TRUE
    if (even_mix) {
      
      alpha <- runif(n = sp, min = 0.5, 5)
      
    } else {
      
      alpha <- rep(even_par, sp)
      
    }
    
    # simualate the relative abundance from the Dirichlet distribution
    RA <- gtools::rdirichlet(n = com, alpha = alpha)
    
    # if species is less than 0.01 then set to zero and round to two decimal places
    RA[RA < 0.01] <- 0
    RA <- round(RA, 2)
    
    # convert the relative abundance distribution into overall biomass
    BIO <- rnorm(n = nrow(RA), mean = bio_m, sd = bio_sd)
    for(i in 1:nrow(RA)) {
      RA[i,] <- RA[i,]*BIO[i]
    }
    
  }
  
  # convert RA into a data.frame
  RA <- 
    apply(RA, 1, function(x) {
      
      df <- data.frame(sp = paste0("sp_", 1:length(x)))
      df$M0 <- x
      
      return(df)
      
    })
  RA <- dplyr::bind_rows(RA, .id = "com")
  
  # calculate the relative abundance as well
  RA <- 
    RA %>%
    group_by(com) %>%
    mutate(pi = M0/sum(M0)) %>%
    ungroup()
  
  # draw RGR and trait values
  traits <- faux::rnorm_multi(n = sp, 
                              mu = mu,
                              sd = sd,
                              r = r, 
                              varnames = c("SLA", "RGR"),
                              empirical = TRUE)
  
  # convert RGR to a data.frame
  traits <- data.frame(sp = paste0("sp_", 1:nrow(traits)),
                       RGR = (traits$RGR/10),
                       SLA = traits$SLA)                    
  
  # draw tf and t0 values
  t0f <- 
    
    lapply(1:sp, function(x) {
      
      if(pheno == "random") {
        
        ymin <- runif(n = 1, min = 0, max = dt/5 )
        ymax <- runif(n = 1, min = (dt/5)*4, max = dt )
        
      } else if (pheno == "equal") {
        
        ymin <- 0
        ymax <- dt
        
      }
      
      # round and sort the min and max values to an integer
      z <- sort(round(c(ymin, ymax), 0))
      
      df <- data.frame(sp = paste0("sp_", x))
      df$t0 <- z[1]
      df$tf <- z[2]
      
      return(df)
      
    } )
  t0f <- dplyr::bind_rows(t0f)
  
  # join these data.frames
  prod.dat <- left_join(left_join(RA, t0f, by = "sp"), 
                        traits, 
                        by = "sp")
  
  return(prod.dat)
  
}

# test the sim_SP() function
sim_SP(com = 20, sp = 10, even_par = 0.25, even_mix = FALSE, com_mat = NA, 
       mu = c(100, 0.5), sd = c(40, 0.2), r = 0.5,
       bio_m = 200, bio_sd = 50,
       dt = 60, pheno = "random")

# test the sim_SP() when using a custom community matrix

# generate a community matrix with only one species per community
RA01 <- 
  sapply(1:20, function(x) {
    y <- rep(0, 10)
    y[sample(1:length(y), 1)] <- 100
    return(y)
  } ) %>%
  t()

sim_SP(com = 20, sp = 10, even_par = 0.25, even_mix = FALSE, com_mat = RA01, 
       mu = c(100, 0.5), sd = c(40, 0.2), r = 0.5,
       bio_m = 200, bio_sd = 50,
       dt = 60, pheno = "random")


# run a set of simulations

# load relevant libraries
library(dplyr)
library(betapart)
library(faux)

# set-up a data.frame for simulations
sim.df <- 
  expand.grid(rep = 1:5, 
              pheno = c("equal", "random"),
              dt = c(100),
              even_par = round(seq(0.1, 20, length.out = 5), 1),
              r = seq(0.05, 0.95, 0.10),
              bio_m = c(1000),
              bio_sd = seq(0, 400, by = 20) ) %>%
  arrange(rep) %>%
  as_tibble()
dim(sim.df)

# set-up input lists
data.full <- vector("list", length = nrow(sim.df))
data.sum <- vector("list", length = nrow(sim.df))
data.out <- vector("list", length = nrow(sim.df))

# loop over each parameter combination
for(i in 1:nrow(sim.df)) {
  
  # get ith row of the data
  df.x <- sim.df[i,]
  
  prod.dat <- 
    sim_SP(com = 20, sp = 10, even_par = df.x[["even_par"]], even_mix = FALSE, com_mat = NA, 
           mu = c(100, 0.5), sd = c(40, 0.2), r = df.x[["r"]],
           bio_m = df.x[["bio_m"]], bio_sd = df.x[["bio_sd"]],
           dt = 100, pheno = df.x[["pheno"]])
  
  # write to a list
  data.full[[i]] <- prod.dat
  
  # calculate SANPP
  prod.dat$ANPP<- with(prod.dat,
                       (M0*(exp(RGR*(tf-t0)) - 1)) )
  
  # summarise to the community-level
  prod.sum <- 
    prod.dat %>%
    group_by(com) %>%
    summarise(sum_M0 = sum(M0),
              cor_dom = cor(RGR, M0),
              CWM_pheno = mean( pi*(tf-t0) ),
              CWM_SLA = sum(pi*SLA),
              CWM_RGR = sum(pi*RGR),
              ANPP = (sum(ANPP)/df.x[["dt"]]) )
  
  # write to a list
  data.sum[[i]] <- prod.sum
  
  # fit a linear model to get the RGR r2 value
  lm.x <- lm(ANPP ~ CWM_RGR, data = prod.sum)
  lm.x <- summary(lm.x)
  r2_RGR <- lm.x$r.squared
  
  # fit a linear model to get the SLA r2 value
  lm.y <- lm((ANPP) ~ CWM_SLA, data = prod.sum)
  lm.y <- summary(lm.y)
  r2_SLA <- lm.y$r.squared
  
  # fit a linear model to get the SLA r2 value on log(ANPP)
  lm.z <- lm(log(ANPP) ~ CWM_SLA, data = prod.sum)
  lm.z <- summary(lm.z)
  r2_ln_SLA <- lm.z$r.squared
  
  # pull this into a data.frame
  df.out <- data.frame(M0_negative = any(prod.dat$M0 < 0),
                       cor_dom = median(prod.sum$cor_dom),
                       r2_SLA = round(r2_SLA, 3),
                       r2_ln_SLA = round(r2_ln_SLA, 3),
                       r2_RGR = round(r2_RGR, 3))
  
  # write to a list
  data.out[[i]] <- df.out
  
}

# convert full output into a data.frame
sim.out.df <- bind_cols(sim.df, bind_rows(data.out))
head(sim.out.df)

# filter rows with M0 that was negative
sim.out.df <- 
  sim.out.df %>%
  filter(M0_negative != TRUE)


# simulate a set of parameters with perfect conditions
sim.df2 <- 
  expand.grid(rep = 1:5, 
              pheno = c("equal"),
              dt = 60,
              even_par = round(seq(0.1, 20, length.out = 5), 1),
              r = 0.99,
              bio_m = c(100) ) %>%
  arrange(rep) %>%
  as_tibble()
dim(sim.df2)

# set-up input lists
data.perfect <- vector("list", length = nrow(sim.df2))

# loop over each parameter combination
for(i in 1:nrow(sim.df2)) {
  
  # get ith row of the data
  df.x <- sim.df2[i,]
  
  RA01 <- 
    sapply(1:20, function(x) {
      y <- rep(0, 10)
      y[sample(1:length(y), 1)] <- df.x[["bio_m"]]
      return(y)
    } ) %>%
    t()
  
  prod.dat <- 
    sim_SP(com = 20, sp = 10, even_par = df.x[["even_par"]], even_mix = FALSE, com_mat = RA01, 
           mu = c(100, 0.5), sd = c(40, 0.2), r = df.x[["r"]],
           bio_m = NA, bio_sd = NA,
           dt = df.x[["dt"]], pheno = df.x[["pheno"]])
  
  # calculate SANPP
  prod.dat$ANPP<- with(prod.dat,
                       (M0*(exp(RGR*(tf-t0)) - 1)) )
  
  # summarise to the community-level
  prod.sum <- 
    prod.dat %>%
    group_by(com) %>%
    summarise(sum_M0 = sum(M0),
              cor_dom = cor(RGR, M0),
              CWM_pheno = mean( pi*(tf-t0) ),
              CWM_SLA = sum(pi*SLA),
              CWM_RGR = sum(pi*RGR),
              ANPP = (sum(ANPP)/df.x[["dt"]]) )
  
  # fit a linear model to get the RGR r2 value
  lm.x <- lm(ANPP ~ CWM_RGR, data = prod.sum)
  lm.x <- summary(lm.x)
  r2_RGR <- lm.x$r.squared
  
  # fit a linear model to get the SLA r2 value
  lm.y <- lm(ANPP ~ CWM_SLA, data = prod.sum)
  lm.y <- summary(lm.y)
  r2_SLA <- lm.y$r.squared
  
  # fit a linear model to get the SLA r2 value on log(ANPP)
  lm.z <- lm(log(ANPP) ~ CWM_SLA, data = prod.sum)
  lm.z <- summary(lm.z)
  r2_ln_SLA <- lm.z$r.squared
  
  # pull this into a data.frame
  df.out <- data.frame(M0_negative = any(prod.dat$M0 < 0),
                       cor_dom = median(prod.sum$cor_dom),
                       r2_SLA = round(r2_SLA, 3),
                       r2_ln_SLA = round(r2_ln_SLA, 3),
                       r2_RGR = round(r2_RGR, 3))
  
  # write to a list
  data.perfect[[i]] <- df.out
  
}

# convert full output into a data.frame
sim.out.df2 <- bind_cols(sim.df2, bind_rows(data.perfect))
head(sim.out.df2)

# remove any where the M0 was negative
sim.out.df <- 
  sim.out.df %>%
  filter(M0_negative != TRUE)

# calculate the average and standard deviation of the r2_SLA for the perfect simulation
m1 <- mean(sim.out.df2$r2_SLA)
print(m1)
sd1 <- sd(sim.out.df2$r2_SLA)
print(sd1)

# load the ggplot2 library
library(ggplot2)

i <- sample(1:length(data.sum), 1)
sim.df[i,]
par(mfrow = c(1, 2))
plot(data.sum[[i]]$CWM_SLA, (data.sum[[i]]$ANPP) )
plot(data.sum[[i]]$CWM_SLA, log(data.sum[[i]]$ANPP) )

mean(sim.out.df$r2_SLA)
mean(sim.out.df$r2_ln_SLA)

ggplot(data = sim.out.df %>% 
         filter(pheno == "equal") %>%
         mutate(`r - (RGR ~ SLA)` = as.character(r)),
       mapping = aes(x = bio_sd, y = r2_SLA, colour = `r - (RGR ~ SLA)`)) +
  geom_segment(mapping = aes(x = 0, xend = 0, y = m1-sd1, yend = m1+sd1),
              colour = "red") +
  geom_point(mapping = aes(x = 0, y = m1),
               colour = "red") +
  geom_jitter(alpha = 0.1, shape = 1, width = 2, height = 0.01) +
  ylab("r2 - (ANPP ~ CWM SLA)") +
  xlab("SD standing biomass (g m-2)") +
  geom_smooth(se = FALSE) +
  scale_colour_viridis_d(option = "C", end = 0.95) +
  scale_y_continuous(breaks = seq(0, 0.90, 0.10)) +
  ggtitle("Phenology: Equal") +
  theme_bw()

ggplot(data = sim.out.df %>% 
         filter(pheno == "random") %>%
         mutate(`r - (RGR ~ SLA)` = as.character(r)),
       mapping = aes(x = bio_sd, y = r2_SLA, colour = `r - (RGR ~ SLA)`)) +
  geom_segment(mapping = aes(x = 0, xend = 0, y = m1-sd1, yend = m1+sd1),
               colour = "red") +
  geom_point(mapping = aes(x = 0, y = m1),
             colour = "red") +
  geom_jitter(alpha = 0.1, shape = 1, width = 2, height = 0.01) +
  ylab("r2 - (ANPP ~ CWM SLA)") +
  xlab("SD standing biomass (g m-2)") +
  geom_smooth(se = FALSE) +
  scale_colour_viridis_d(option = "C", end = 0.95) +
  scale_y_continuous(breaks = seq(0, 0.90, 0.10)) +
  ggtitle("Phenology: Random") +
  theme_bw()

### END
