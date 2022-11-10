
# load relevant libraries
library(dplyr)
library(betapart)
library(faux)

# we can set even_mix = TRUE to see if this changes the results
# otherwise we stick with communities that vary in their dominance
sp = 5
com = 10
even_par = 0.5
even_mix = FALSE

com_mat = NA # allow the user to add their own community matrix (sites are rows, species are columns)

# parameters for the multivariate normal and the chosen correlation coefficient
mu = c(10, 0.5)
sd = c(4, 0.2)
r = 0.5

# phenology
pheno = "random" 
pheno = "equal"

# duration of measurements
dt <- 150

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
  
}

# if species is less than 0.01 then set to zero and round to two decimal places
RA[RA < 0.01] <- 0
RA <- round(RA, 2)

# convert RA into a data.frame
RA <- 
  apply(RA, 1, function(x) {
    
    df <- data.frame(sp = paste0("sp_", 1:length(x)))
    df$pi <- x
    
    return(df)
    
  })
RA <- dplyr::bind_rows(RA, .id = "com")

# draw RGR and trait values
traits <- faux::rnorm_multi(n = sp, 
                            mu = mu,
                            sd = sd,
                            r = r, 
                            varnames = c("SLA", "RGR"),
                            empirical = TRUE)

# convert RGR to a data.frame
traits <- data.frame(sp = paste0("sp_", 1:nrow(traits)),
                     RGR = traits$RGR,
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

# calculate phenological overlap
pheno_over <- beta.multi.abund(x = t0f[,-1], index.family = "bray")$beta.BRAY

# join these data.frames
prod.dat <- left_join(left_join(RA, t0f, by = "sp"), 
                      traits, 
                      by = "sp")

# calculate SANPP
prod.dat$SANPP<- with(prod.dat,
                      (pi*exp(RGR*(tf-t0))) )
  
# summarise to the community-level
prod.sum <- 
  prod.dat %>%
  group_by(com) %>%
  summarise(cor_dom = cor(RGR, pi),
            CWM_SLA = sum(pi*SLA),
            CWM_RGR = sum(pi*RGR),
            SANPP = (log(sum(SANPP))/dt)*1000 )

# fit a linear model to get the RGR r2 value
lm.x <- lm(SANPP ~ CWM_RGR, data = prod.sum)
lm.x <- summary(lm.x)
r2_RGR <- lm.x$r.squared

# fit a linear model to get the SLA r2 value
lm.y <- lm(SANPP ~ CWM_SLA, data = prod.sum)
lm.y <- summary(lm.y)
r2_SLA <- lm.y$r.squared

# pull this into a data.frame
df.out <- data.frame(cor_dom = median(prod.sum$cor_dom),
                     pheno_over = 1 - round(pheno_over, 2),
                     r2_SLA = round(r2_SLA, 3),
                     r2_RGR = round(r2_RGR, 3),
                     r_RGR = round(cor(prod.sum$CWM_RGR, prod.sum$SANPP), 3))

print(df.out)


plot(traits$RGR, )



return(df.out)


