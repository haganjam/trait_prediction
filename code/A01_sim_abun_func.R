#' 
#' @title sim_abun()
#' 
#' @description Function to simulate the abundance of sp species in com communities.
#' The function uses a lognormal distribution from the mobsim package to generate
#' a species abundance distribution.
#' 
#' @param sp - number of species to simulate
#' @param com - number of different communities to simulate
#' @param N - number of individuals per community
#' @param cv_abund - eveness parameter for the lognormal distribution
#' 
#' @return Returns a data.frame with community id, species name, number of individuals
#' and relative abundance.
#' 

# simulate a set of communities with some relative abundances
sim_abund <- function(sp, com, N, cv_abund = NA) {

  sad_list <- 
    
    lapply(1:com, function(x) {
      
      # sp_pool size determined by J and theta together (can we calculate this?)
      sad.x <- mobsim::sim_sad(s_pool = sp, n_sim = N, sad_type = c("lnorm"),
                               sad_coef = list("cv_abund" = cv_abund))
      
      # remove the class attributes
      attr(sad.x, "class") <- NULL
      names(sad.x) <- NULL
      
      # initialise an empty vector with 10 species
      sp_vec <- rep(0, sp)
      sp_vec[sample(1:sp, size = length(sad.x))] <- sad.x
      
      # wrap into a data.frame
      sp_df <- dplyr::tibble(com = x, sp = paste0("sp_", 1:sp), abund = sp_vec)
      
      # convert to relative abundance
      sp_df$abund <- with(sp_df, abund/sum(abund))
      
      return(sp_df)
      
    } )
  
  return(sad_list)

  }

### END
