
#' @title: Theoretical limits to the correlation between plant traits and productivity
#' 
#' @description: Here, we use the model of primary productivity based on relative growth
#' rates of individual species and their biomass proportions in the communities to test the
#' theoretical limits of the predictability of primary productivity from plant functional traits.
#' 
#' @details: Garnier et al. (2004, Ecology) proposed a formal, causal model linking species'
#' relative growth rates and phenology to primary productivity:
#' 
#' SANPP (g g-1 day-1) = log10( sum(pi)*exp(RGRi*(tf-t0)i) ) )/delta time
#' 
#' where: 
#' pi is standing, living biomass in a community of species i
#' RGR is the maximum relative growth rate (g g-1 day-1) determined experimentally
#' tf-t0 is the period of growth for species i
#' delta time is the period over which the measurements too place
#' 
#' Thus, we can use this equation to simulate SANPP based on known RGRs, pi's and
#' relative growth rates. 
#' 
#' We then assume as per Garnier et al. (2004) and Reich et al. (1997) that certain
#' easy-to-measure plant traits like SLA, leaf nitrogen etc. can predict RGR. 
#' 
#' We then can ask if there is uncertainty in that prediction along with phenological
#' uncertainty and uncertainty in the actual SANPP measurements, what the theoretical
#' limits to predicting SANPP from plant traits is.
#' 

# load relevant libraries
library(dplyr)

# communities
com <- 5

# species number
sp <- 5

# evenness parameter
even <- 3

# set the delta time period
dt <- 150

# each row is a different community
pi <- gtools::rdirichlet(n = com, alpha = rep(even, sp))
pi <- round(pi, 2)

# draw RGR values
RGR <- round(runif(n = sp, min = 0.05, 0.5), 2)
RGR <- data.frame(sp = paste0("sp_", 1:length(RGR)),
                  RGR = RGR)

# draw tf and t0 values
t0f <- lapply(1:sp, function(x) {
  y <- runif(n = 2, min = 0, max = dt)
  z <- sort(round(y, 0))
  
  df <- data.frame(sp = paste0("sp_", x))
  df$t0 <- z[1]
  df$tf <- z[2]
  
  return(df)
  
} )
t0f <- dplyr::bind_rows(t0f)

pi <- 
  apply(pi, 1, function(x) {
  df <- data.frame(sp = paste0("sp_", 1:length(x)))
  df$pi <- x
  
  return(df)
  
})
pi <- dplyr::bind_rows(pi, .id = "com")

# join these data.frames
prod.dat <- left_join(left_join(pi, t0f, by = "sp"), 
                      RGR, by = "sp"
                      )

# calculate SANPP
SANPP <- vector(length = nrow(prod.dat))
for(i in 1:nrow(prod.dat)) {
    
    SANPP[i] <- with(prod.dat[i,],
                     log10( sum( (pi*exp(RGR*(tf-t0))) ) )/dt )
    
  }

# add SANPP to the dataset
prod.dat$SANPP <- SANPP
      







