#' 
#' @title sim_traits()
#' 
#' @description Function to simulate pairs of correlated traits using multivariate
#' normal distributions. The function first simulates trait means for sp different species
#' from a multivariate normal distribution for two traits. Then, it simulates individual
#' level traits from a multivariate normal distribution of the two traits parameterised
#' with the species' trait means along with additional variance terms. The correlation of traits
#' within in an individual is set at zero.
#' 
#' @param sp - number of species to simulate
#' @param com - number of different communities to simulate
#' @param t1 - name of the first trait
#' @param t2 - name of the second trait
#' @param mu_t1 - mean of t1 in a multivariate normal distribution
#' @param mu_t2 - mean of t2 in a multivariate normal distribution
#' @param sd_t1 - sd of t1 in a multivariate normal distribution
#' @param sd_t2 - sd of t2 in a multivariate normal distribution
#' @param r - empirical correlation between t1 and t2
#' @param Eind_t1 - individual level sd in the t1 trait
#' @param Eind_t2 - individual level sd in the t2 trait
#' 
#' @return Functions returns a data.frame with the community id number, species name
#' and the values of t1 and t2.
#' 

sim_traits <- function(sp = 10, com = 20,
                       t1, t2, mu_t1, mu_t2, sd_t1, sd_t2, r,
                       Eind_t1, Eind_t2, 
                       t2_scale = 10) {
  
  # draw RGR and trait values
  t12_mu <- faux::rnorm_multi(n = sp, 
                              mu = c(mu_t1, mu_t2),
                              sd = c(sd_t1, sd_t2),
                              r = r, 
                              varnames = c(t1, t2),
                              empirical = TRUE
  )
  
  # round off the values
  t12_mu <- round(t12_mu, 3)
  
  # if less than 0, then set to zero
  t12_mu[t12_mu < 0] <- 0
  
  # add a species column
  t12_mu <- cbind(data.frame(sp = paste0("sp_", 1:sp)), t12_mu)
  
  # draw an individual of each species for each community
  t12_ind <- vector("list", length = sp)
  for(i in 1:sp) {
    
    # extract the first row
    df.x <- t12_mu[i,]
    
    # sample from the multivariate normal distribution
    x <- faux::rnorm_multi(n = com,
                           mu = c( df.x[[t1]], df.x[[t2]]  ),
                           sd = c(Eind_t1, Eind_t2),
                           r = 0, 
                           varnames = c(t1, t2),
                           empirical = FALSE
    )
    
    # bind community and species information
    y <- cbind(data.frame(com = (1:com), sp = paste0("sp_", i)), x)
    
    # add to list
    t12_ind[[i]] <- y
    
  }
  
  # bind into a large data.frame
  t12_ind <- dplyr::arrange( dplyr::bind_rows(t12_ind), com )
  
  # scale trait2
  t12_ind[[t2]] <- t12_ind[[t2]]/t2_scale
  
  # split by community
  t12_ind <- split(t12_ind, t12_ind$com)
  
  # remove the names of the list
  names(t12_ind) <- NULL
  
  # return this data.frame
  return(t12_ind)
  
}

### END
