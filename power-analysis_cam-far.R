#### POWER ANALYSIS: cameras vs FAR ####


## Development version doesn't seem necessary:
#  library(devtools)
#  devtools::install_github("pitakakariki/simr")
#  with_libpaths(new = "C:/Users/sittersh/Documents/R/lib", install_github("pitakakariki/simr"))
#  library ("simr", lib.loc="~Documents/R/x86_64-pc-linux-gnu-library/3.4")
#  ~Documents/R/x86_64-pc-linux-gnu-library/3.4


#### The thing to do is ####

# Specify effect sizes for three of the treatment levels in relation to the reference level 
# Simulate data as desired, say 10 samples per treatment
# Perform a likelihood ratio test
# Repeat 1000 times to give a single power estimate with confidence intervals.
# The order of the levels must matter but at least the LR gives an overall number.

# s.tr.d = small effect size, plots per treatment, deer
# m.tr.d = medium effect size, plots per treatment vary, deer
# l.tr.d = large effect size, plots per treatment vary, deer


#### Power scenarios ####

# 2 x methods
# 2 x data types
# 5 x species
# 3 x effect sizes
# 2 x sample-size scenarios (treatment and swale)


#### Labels and effect sizes ####

#  FAR.y = raw count
#  FARp.y = count per day
#  FARm.y = count per day as proportion of the max. count per day
#  FARb.y = presence-absence data

## CB09 is the reference level.  Other treatments are 'open'.
## Effect sizes of interest: 
#    small (0.2) = qlogis -1.386294
#    medium (0.4) = qlogis -0.4054651
#    large (0.6) = qlogis 0.4054651
#    extra large (0.8) = qlogis 1.386294


#### GLMM error-solving ####

# https://rstudio-pubs-static.s3.amazonaws.com/33653_57fc7b8e5d484c909b615d8633c01d51.html


#### ANALYSIS ####


#### Packages and inputs ####

library(simr)
library(lme4)

setwd("~/Documents/git/Power-analysis_HSM")
cam <- read.table("Cam_T2.txt", header=TRUE) 
far <- read.table("FAR_T2.txt", header=TRUE) 


#### Presence-absence models ####

# GLM
gm1 <- glm (FARb.w ~ treatment, family = binomial, data=far)

# GLMM
gm2 <- glmer (FARb.w ~ treatment + (1|swale), family = binomial, data=far,
              control = glmerControl(optimizer="bobyqa"))
# Effect sizes
fixef(gm1)[2:4] <- c(-1.386294, -1.386294, -1.386294) # small
fixef(gm1)[2:4] <- c(-0.4054651, -0.4054651, -0.4054651) # medium
fixef(gm1)[2:4] <- c(0.4054651, 0.4054651, 0.4054651) # large


## BINOMIAL DEER
gm1 <- glmer (FARb.d ~ treatment + (1|swale), family = binomial, data=far,
              control = glmerControl(optimizer="bobyqa"))


power.it <- function (z) { # z = glmer
  require (simr)
  require (lme4)
  # SMALL
  fixef(z)[2:4] <- c(-1.386294, -1.386294, -1.386294) 
  # treatment
  ext.s.t <- extend (z, within="treatment", n = 160) 
  pc.s.t <- powerCurve (ext.s.t, within="treatment", breaks = c(5, 10, 20, 40, 80, 160))
  pc.s.t <- summary(pc.s.t)
  # swale
  ext.s.s <- extend (z, along="swale", n=64) 
  pc.s.s <- powerCurve(ext.s.s, along="swale", breaks = c(4, 8, 16, 32, 64))
  pc.s.s <- summary(pc.s.s)
  # MEDIUM
  fixef(z)[2:4] <- c(-0.4054651, -0.4054651, -0.4054651)
  # treatment
  ext.m.t <- extend (z, within="treatment", n = 160) 
  pc.m.t <- powerCurve (ext.m.t, within="treatment", breaks = c(5, 10, 20, 40, 80, 160)) 
  pc.m.t <- summary(pc.m.t)
  # swale
  ext.m.s <- extend (z, along="swale", n=64) 
  pc.m.s <- powerCurve(ext.m.s, along="swale", breaks = c(4, 8, 16, 32, 64))
  pc.m.s <- summary(pc.m.s)
  # LARGE
  fixef(z)[2:4] <- c(0.4054651, 0.4054651, 0.4054651)
  # treatment
  ext.l.t <- extend (z, within="treatment", n = 160) 
  pc.l.t <- powerCurve (ext.l.t, within="treatment", breaks = c(5, 10, 20, 40, 80, 160))
  pc.l.t <- summary(pc.l.t)
  # swale
  ext.l.s <- extend (z, along="swale", n=64) 
  pc.l.s <- powerCurve(ext.l.s, along="swale", breaks = c(4, 8, 16, 32, 64))
  pc.l.s <- summary(pc.l.s)
  # OUT
  out <- rbind(pc.s.t, pc.s.s, # small
                         pc.m.t, pc.m.s, # medium
                         pc.l.t, pc.l.s) # large
  return (out)
}

#scenario <- data.frame(scenario= factor(levels = c("pc.s.t_005", "pc.s.t_010", "pc.s.t_020", 
#                                                   "pc.s.t_040", "pc.s.t_080", "pc.s.t_160",
#                                                  "pc.s.s_04", "pc.s.s_08", "pc.s.s_s16",
#                                                 "pc.s.s_32", "pc.s.s_64",
#                                                "pc.m.t_005", "pc.m.t_010", "pc.m.t_020", 
#                                               "pc.m.t_040", "pc.m.t_080", "pc.m.t_160",
#                                              "pc.m.s_04", "pc.m.s_08", "pc.m.s_s16",
#                                             "pc.m.s_32", "pc.m.s_64",
#                                             "pc.l.t_005", "pc.l.t_010", "pc.l.t_020", 
#                                            "pc.l.t_040", "pc.l.t_080", "pc.l.t_160",
"pc.l.s_04", "pc.l.s_08", "pc.l.s_s16",
"pc.l.s_32", "pc.l.s_64")))

pit <- power.it (gm1)

# Tested on binomial deer

power.test <- function (z) {  ## power.test works; z = glmer. Changed last rbind to cbind
  require (simr)
  require (lme4)
  # LARGE
  fixef(z)[2:4] <- c(0.4054651, 0.4054651, 0.4054651)
  # treatment
  ext.l.t <- extend (z, within="treatment", n = 160) 
  pc.l.t <- powerCurve (ext.l.t, within="treatment", breaks = c(5, 160))
  pc.l.t <- summary(pc.l.t)
  # swale
  ext.l.s <- extend (z, along="swale", n=64) 
  pc.l.s <- powerCurve(ext.l.s, along="swale", breaks = c(4, 64))
  pc.l.s <- summary(pc.l.s)
  # OUT
  all.out <- rbind.data.frame(pc.l.t, pc.l.s) # large
  scenario <- data.frame(scenario= factor(levels = c("pc.l.t_005", "pc.l.t_160",
                                                     "pc.l.s_04", "pc.l.s_64")))
  out <- cbind (scenario, all.out)
  return (out)
}

pt <- power.test (gm1)

write.table(out, "deer_far_test.txt")
