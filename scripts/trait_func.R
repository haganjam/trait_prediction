

# set-up the trait and parameter correlations

# how many species to simulate
sp <- 10

# how many communities to simulate
com <- 20

# set-up the parameters for the overall multivariate normal distribution for the traits
mu_SLA <- 100
sd_SLA <- 40

mu_RGR <- 0.5
sd_RGR <- 0.2

mu_H <- 200
sd_H <- 60

mu_M0 <- (1000/sp)
sd_M0 <- 20

r_SLA_RGR <- 0.5
r_H_M0 <- 0.3


# set-up error in the individual distributions of these traits
Eind_SLA <- 5
Eind_RGR <- 0.05
Eind_H <- 10
Eind_M0 <- 20

# write this into a function so we can simulate any trait and parameter that we want

# draw RGR and trait values
SLA_RGR_mu <- faux::rnorm_multi(n = sp, 
                                mu = c(mu_SLA, mu_RGR),
                                sd = c(sd_SLA, sd_RGR),
                                r = r_SLA_RGR, 
                                varnames = c("SLA", "RGR"),
                                empirical = TRUE)
                               
# round off the values
SLA_RGR_mu <- round(SLA_RGR_mu, 3)

# add a species column
SLA_RGR_mu <- cbind(data.frame(sp = paste0("sp_", 1:sp)), SLA_RGR_mu)

# draw an individual of each species for each community
SLA_RGR_ind <- vector("list", length = sp)
for(i in 1:sp) {
  
  # extract the first row
  df.x <- SLA_RGR_mu[i,]
  
  # sample from the multivariate normal distribution
  x <- faux::rnorm_multi(n = com,
                         mu = c( df.x[["SLA"]], df.x[["RGR"]]  ),
                         sd = c(Eind_SLA, Eind_RGR),
                         r = 0, 
                         varnames = c("SLA", "RGR"),
                         empirical = FALSE
  )
  
  # bind community and species information
  y <- cbind(data.frame(com = 1:com, sp = paste0("sp_", i)), x)
  
  # add to list
  SLA_RGR_ind[[i]] <- y
  
}

# bind into a large data.frame
SLA_RGR_ind <- dplyr::bind_rows(SLA_RGR_ind)










# draw RGR and trait values
H_M0_mu <- faux::rnorm_multi(n = sp, 
                             mu = c(mu_H, mu_M0),
                             sd = c(sd_H, sd_M0),
                             r = r_H_M0, 
                             varnames = c("H", "M0"),
                             empirical = TRUE
                             )

# round off the values
H_M0_mu <- round(H_M0_mu, 3)

# add a species column
H_M0_mu <- cbind(data.frame(sp = paste0("sp_", 1:sp)), H_M0_mu)

# draw an individual of each species for each community
H_M0_ind <- vector("list", length = sp)
for(i in 1:sp) {
  
  # extract the first row
  df.x <- H_M0_mu[i,]
  
  # sample from the multivariate normal distribution
  x <- faux::rnorm_multi(n = com,
                         mu = c( df.x[["SLA"]], df.x[["RGR"]]  ),
                         sd = c(Eind_SLA, Eind_RGR),
                         r = 0, 
                         varnames = c("SLA", "RGR"),
                         empirical = FALSE
  )
  
  # bind community and species information
  y <- cbind(data.frame(com = 1:com, sp = paste0("sp_", i)), x)
  
  # add to list
  SLA_RGR_ind[[i]] <- y
  
}

# bind into a large data.frame
SLA_RGR_ind <- dplyr::bind_rows(SLA_RGR_ind)




