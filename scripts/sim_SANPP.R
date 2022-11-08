
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
#' notes:
#' 1. seems that one major part of the bad correlation is when the correlation
#' between species' pi and the relative growth rate is low
#' 
#' 2. but definitely not only that... the differences in growing season length
#' are also super important it seems
#' 
#' @function - sim_SANPP()
#' @param com - number of communities to simulate
#' @param sp - number of species to simulate
#' @param even - evenness parameter corresponding to the alpha parameter in the Dirichlet distribution
#' @param dt - growing season length (days)
#' @param mu - vector of mean values for SLA and RGR
#' @param sd - vector of sd values for SLA and RGR
#' @param r -  correlation coefficient (Pearson) between SLA and RGR
#' @param D - modifies the relative abundances according to the following rules:
#' D = "cor then species relative abundances correlate with their RGR
#' D = "dom" then one species completely dominates each community
#' D = "random" leaves the relative abundances unmodified
#' @param phen - either picks random phenology or keeps all phenologies equal for all species
#'  1. phen = "random"
#'  2. phen = "equal"
#' 

sim_SANPP <- function(com = 20, sp = 10, even = 0.5, dt = 150,
                      mu = c(100, 0.5), sd = c(40, 0.2), r = 0.5,
                      D = "random", phen = "random") {
  
  # load relevant libraries
  library(dplyr)
  library(betapart)
  library(faux)
  
  # each row is a different community
  pi <- gtools::rdirichlet(n = com, alpha = rep(even, sp))
  pi <- as.matrix(pi)
  pi[pi < 0.005] <- 0
  
  # draw RGR and trait values
  traits <- rnorm_multi(n = sp, 
                        mu = mu,
                        sd = sd,
                        r = r, 
                        varnames = c("SLA", "RGR"),
                        empirical = TRUE
  )
  
  # convert RGR to a data.frame
  traits <- data.frame(sp = paste0("sp_", 1:nrow(traits)),
                       RGR = traits$RGR,
                       SLA = traits$SLA
  )
  
  # modify abundances
  if (D == "cor") {
    
    pi <- 
      apply(pi, 1, function(x) {
        
        # sort out the ranks when there are duplicates
        if( length(x) != length(unique(x)) ) {
          
          # rank the values
          u <- sapply(x, function(vec) which(vec == sort(x)) )
          
          # extract any instances where many are matched
          uv <- u[lapply(u, length) > 1]
          
          # fill in the cases uniquely with multiple matches
          u[lapply(u, length)>1] <- unique(unlist(uv))
          
          # return an unlisted vector
          xs <- unlist(u)
          
        } else {
          
          xs <- rank(x)
          
        }
        
        # get the relative growth rates
        RGRs <- rank(traits$RGR)
        
        # match the rank of the growth rates with the relative abundances
        y <- x[order(match(xs, RGRs))]
        
        return(y)
        
      })
    
    # transpose the matrix
    pi <- t(pi)
    
  } else if (D == "dom") {
    
    pi <- 
      apply(pi, 1, function(x) {
        
        y <- x
        y[sample(1:length(y), 1)] <- 1
        y <- ifelse(y != 1, 0, y)
        
        return(y)
        
      })
    
    # transpose the matrix
    pi <- t(pi)
    
  }
  
  # draw tf and t0 values
  t0f <- lapply(1:sp, function(x) {
    
    if(phen == "random") {
      ymin <- runif(n = 1, min = 0, max = dt/5 )
      ymax <- runif(n = 1, min = (dt/5)*4, max = dt )
    } else if (phen == "equal") {
      ymin <- 0
      ymax <- 150
    }
    
    z <- sort(round(c(ymin, ymax), 0))
    
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
                        traits, by = "sp"
  )
  
  # calculate SANPP
  SANPP <- vector(length = nrow(prod.dat))
  for(i in 1:nrow(prod.dat)) {
    
    SANPP[i] <- with(prod.dat[i,],
                     (pi*exp(RGR*(tf-t0))) )
    
  }
  
  # add SANPP to the dataset
  prod.dat$SANPP <- SANPP
  
  # summarise to the community-level
  prod.sum <- 
    prod.dat %>%
    group_by(com) %>%
    summarise(CWM_SLA = sum(pi*SLA),
              CWM_RGR = sum(pi*RGR),
              SANPP = log10(sum(SANPP))/dt )
  
  # fit a linear model to get the RGR r2 value
  lm.x <- lm(SANPP ~ CWM_RGR, data = prod.sum)
  lm.x <- summary(lm.x)
  r2_RGR <- lm.x$r.squared
  
  # fit a linear model to get the SLA r2 value
  lm.y <- lm(SANPP ~ CWM_SLA, data = prod.sum)
  lm.y <- summary(lm.y)
  r2_SLA <- lm.y$r.squared
  
  # check if species with high RGR are dominating across all communities
  cor_dom1 <- cor(prod.dat$RGR, prod.dat$pi, method = "spearman")
  
  # check if species with high RGR are dominating in individual communities on average
  cor_dom2 <- 
    prod.dat %>%
    group_by(com) %>%
    summarise(cor1 = cor(pi, RGR, method = "spearman")) %>%
    pull(cor1) %>%
    mean()
  
  # phenological overlap
  phen_over <- beta.multi.abund(x = t0f[,-1], index.family = "bray")$beta.BRAY
  
  # pull this into a data.frame
  df.out <- data.frame(cor_dom1 = round(cor_dom1, 2),
                       cor_dom2 = round(cor_dom2, 2) ,
                       phen_over = 1 - round(phen_over, 2),
                       r2_SLA = round(r2_SLA, 3),
                       r2_RGR = round(r2_RGR, 3) )
  
  return(df.out)
  
}

# test the function
sim_SANPP(com = 20, sp = 10, even = 0.5, dt = 150,
          mu = c(100, 0.5), sd = c(40, 0.2), r = 0.5,
          D = "random", phen = "equal")

# set-up a data.frame for simulations
sim.df1 <- data.frame(D = c("dom", "random"),
                      phen = c("equal", "random"))
sim.df2 <- 
  expand.grid(rep = 1:10, 
              even = round(seq(0.5, 50, length.out = 20), 1),
              r = seq(0.05, 0.95, 0.05))
sim.df <- tibble(cbind(sim.df1, sim.df2))
sim.df <- 
  sim.df %>%
  arrange(rep)

# loop over these parameter values
sim.out <- vector("list", length = nrow(sim.df))
for(i in 1:nrow(sim.df)) {
  
  pars <- sim.df[i,]
  
  x <- 
    sim_SANPP(com = 20, sp = 10, even = pars[["even"]], dt = 150,
              mu = c(100, 0.5), sd = c(40, 0.2), r = pars[["r"]],
              D = pars[["D"]], phen = pars[["phen"]])
  
  sim.out[[i]] <- x
  
}

# convert the output into a data.frame
sim.out.df <- bind_cols(sim.df, bind_rows(sim.out))

# load the ggplot2 library
library(ggplot2)

ggplot(data = sim.out.df %>% mutate(r = as.character(r)),
       mapping = aes(x = phen_over, y = r2_SLA, colour = r)) +
  geom_point(alpha = 0.1) +
  geom_smooth(se = FALSE) +
  facet_wrap(vars(D, phen)) +
  theme_bw()

ggplot(data = sim.out.df %>% filter(D == "dom", phen == "equal"),
       mapping = aes(x = r, y = r2_SLA, colour = even)) +
  geom_jitter(width = 0.015, alpha = 0.5, shape = 1) +
  scale_colour_viridis_c(option = "C") +
  geom_smooth(colour = "black", size = 0.5) +
  ylab("r2: SANPP ~ CWM-SLA") +
  xlab("r: RGR ~ SLA") +
  theme_bw()

ggplot(data = sim.out.df %>% filter(D == "random", phen == "random"),
       mapping = aes(x = r, r2_SLA, colour = even )) +
  geom_jitter(width = 0.015, alpha = 0.5, shape = 1) +
  scale_colour_viridis_c(option = "C") +
  ylab("r2: SANPP ~ CWM-SLA") +
  xlab("r: RGR ~ SLA") +
  geom_smooth(colour = "black", size = 0.5) +
  theme_classic()

### END
