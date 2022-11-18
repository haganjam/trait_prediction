#' 
#' @title sim_pheno()
#' 
#' @description Function to simulate species' mean phenological period from a Poisson
#' distribution and individual level phenological period from a normal distribution with
#' a defined standard deviation.
#' 
#' @param sp - number of species to simulate
#' @param com - number of different communities to simulate
#' @param pheno_mu - mean of the normal distribution for defining species' phenological period average
#' @param pheno_sd - sd of the normal distribution for defining species' phenological period
#' @param Eind_pheno - individual level variation in phenology
#' 
#' @return Returns a data.frame with community id, species name, phenological period
#' 

sim_pheno <- function(com = 20, sp = 10, 
                      pheno_mu = 60, pheno_sd = 5,
                      Eind_pheno) {
  
  # simulate a set of species' mean phenologies from a poisson distribution
  mu_pheno <- rnorm(n = sp, mean = pheno_mu, sd = pheno_sd)
  
  # make sure all are greater than zero
  mu_pheno[mu_pheno < 10] <- 10
  
  # simulate individual level phenologies
  pheno_ind <- vector("list", length = sp)
  for(i in 1:sp) {
    
    # draw individual phenology values around the average phenology
    x <- data.frame(t0f = round(rnorm(n = com, mean = mu_pheno[i], sd = Eind_pheno), 0)) 
    
    # make sure all the phenology values are greater than 1
    x[x < 1] <- 1
    
    # bind community and species information
    y <- cbind(data.frame(com = as.character(1:com), sp = paste0("sp_", i)), x)
    
    # add to list
    pheno_ind[[i]] <- y
    
  }
  
  # bind into large data.frame
  pheno_ind <- bind_rows(pheno_ind)
  
  # arrange by community
  pheno_ind <- arrange(pheno_ind, com)
  
  # return the data.frame
  return(pheno_ind)
  
}

### END
