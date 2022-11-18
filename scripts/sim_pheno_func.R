#' 
#' @title sim_pheno()
#' 
#' @description Function to simulate species' mean phenological period from a Poisson
#' distribution and individual level phenological period from a normal distribution with
#' a defined standard deviation.
#' 
#' @param sp - number of species to simulate
#' @param com - number of different communities to simulate
#' @param pheno_lambda - mean of the Poisson distribution for defining species' phenological period average
#' @param pheno_equal - if true, then all phenological periods for species are set at pheno_lambda 
#' @param pheno_sd - individual level variation in phenology
#' 
#' @return Returns a data.frame with community id, species name, phenological period
#' 

sim_pheno <- function(com = 20, sp = 10, 
                      pheno_lambda = 60, pheno_equal = FALSE, pheno_sd = 5) {
  
  # simulate a set of species' mean phenologies from a poisson distribution
  if(pheno_equal) {
    mu_pheno <- rep(sp, pheno_lambda)
  } else {
    mu_pheno <- rpois(n = sp, lambda = pheno_lambda)
  }
  
  # simulate individual level phenologies
  pheno_ind <- vector("list", length = sp)
  for(i in 1:sp) {
    
    x <- data.frame(pheno = round(rnorm(n = com, mean = mu_pheno[i], sd = pheno_sd), 0)) 
    
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
