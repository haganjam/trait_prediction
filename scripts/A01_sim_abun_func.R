#' 
#' @title sim_abun()
#' 
#' @description Function to simulate the abundance of sp species in com communities.
#' The function uses the Dirichlet distribution to simulate species relative abundances.
#' These relative abundances are then multiplied by N, the number of individuals which
#' we draw from a Poisson distribution.
#' 
#' @param sp - number of species to simulate
#' @param com - number of different communities to simulate
#' @param com_mat - custom community matrix with sites as rows and species as columns
#' @param even_par - evenness parameter corresponding to the alpha parameter in the Dirichlet distribution
#' @param even_mix - evenness parameter is the same for all communities (FALSE) or it differs (TRUE)
#' @param lambda - mean of the poisson distribution
#' @param lambda_equal - all communities have the same number of individuals
#' 
#' @return Returns a data.frame with community id, species name, number of individuals
#' and relative abundance.
#' 

sim_abund <- function(sp = 10, com = 20, com_mat = FALSE,
                      even_par = 0.5, even_mix = FALSE, 
                      lambda = 100, lambda_equal = TRUE) {
  
  # load the dplyr library
  library(dplyr)
  
  # simulate a relative abundance distribution
  if (is.matrix(com_mat)) { 
    
    RA <- com_mat
    
  } else {
    
    # get the alpha values randomly if even_mix == TRUE
    if (even_mix) {
      
      alpha <- runif(n = sp, min = 0.1, 20)
      
    } else {
      
      alpha <- rep(even_par, sp)
      
    }
    
    # simualate the relative abundance from the Dirichlet distribution
    RA <- gtools::rdirichlet(n = com, alpha = alpha)
    
    # if species is less than 0.01 then set to zero and round to two decimal places
    RA[RA < 0.01] <- 0
    RA <- round(RA, 2)
    
    # draw a total number of individuals from a Poisson distribution
    if (lambda_equal) {
      N <- rep(lambda, nrow(RA))
    } else {
      N <- rpois(n = nrow(RA), lambda = lambda)
    }
    
    # convert the relative abundance distribution into overall biomass
    for(i in 1:nrow(RA)) {
      RA[i,] <- round(RA[i,]*N[i], 0)
    }
    
  }
  
  # convert RA into a data.frame
  RA <- 
    apply(RA, 1, function(x) {
      
      df <- data.frame(sp = paste0("sp_", 1:length(x)))
      df$N <- x
      
      return(df)
      
    })
  RA <- dplyr::bind_rows(RA, .id = "com")
  
  # calculate the relative abundance as well
  RA <- 
    RA %>%
    group_by(com) %>%
    mutate(pi = N/sum(N)) %>%
    ungroup()
  
  # return the data.frame
  return(RA)
  
}

### END
