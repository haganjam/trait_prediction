
# set the number of species and communities (sp = n communities)
sp <- 3

# two scenarios:

# scenario 1: Different species in different communities
abun_df <- 
  
  lapply(1:sp, function(x) {
    
    y <- rep(0, sp)
    y[x] <- 1
    
    df <- data.frame(scenario = 1,
                     com = x,
                     sp = paste0("sp_", 1:sp),
                     RA = y)
    
    return(df)
  
}) 
abun_df <- bind_rows(abun_df)

# draw RGR and SLA values and trait values
t12_mu <- faux::rnorm_multi(n = sp, 
                            mu = c(100, 0.05),
                            sd = c(20, 0.02),
                            r = 0.99, 
                            varnames = c("SLA", "RGR"),
                            empirical = TRUE
                            )
t12_mu <- bind_cols(data.frame(sp = paste0("sp_", 1:sp), t12_mu))

# plot the relationship between RGR and SLA 
plot(t12_mu[,c(2, 3)])

# get the phenology
pheno_df <- data.frame(sp = paste0("sp_", 1:sp),
                       t0f = 60)

# bind the phenology to the trait data
trait_df <- full_join(t12_mu, pheno_df, by = "sp")

# get the total standing biomass of each community
M0_df <- 
  data.frame(com = 1:sp, M0 = c(50, 100, 150))

# bind these data.frames into one large data.frame
scen.x <- 
  full_join( full_join(scen1_abun, M0_df, by = "com"),
           trait_df, 
           by = "sp") %>%
  mutate(M0 = M0*RA)
print(scen.x)

# calculate ANPP
scen.x$ANPP<- with(scen.x,
                  (M0*(exp(RGR*(t0f)) - 1)) )
print(scen.x)

# summarise to the community-level
scen.x <- 
  scen.x %>%
  group_by(com) %>%
  summarise(CWM_SLA = sum(RA*SLA),
            CWM_RGR = sum(RA*RGR),
            ANPP = (sum(ANPP)/90), .groups = "drop" )
print(scen.x)

# fit a linear model to get the RGR r2 value
lm.x <- lm(log(ANPP) ~ CWM_SLA, data = scen.x)
summary(lm.x)




