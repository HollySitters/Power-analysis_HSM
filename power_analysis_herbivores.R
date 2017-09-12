library(simr)
library(lme4)

library(devtools)
devtools::install_github("pitakakariki/simr")
with_libpaths(new = "C:/Users/sittersh/Documents/R/lib", install_github("pitakakariki/simr"))
library ("simr", lib.loc="~Documents/R/x86_64-pc-linux-gnu-library/3.4")
 ~Documents/R/x86_64-pc-linux-gnu-library/3.4

## FAR input

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


setwd("~/Documents/git/Herbivore-survey-methods/5_Power analysis")
data <- read.table("Cam_T2.txt", header=TRUE) 
data <- read.table("Cam_T3.txt", header=TRUE) 
data <- read.table("FAR_T2-T3.txt", header=TRUE) 
data <- read.table("FAR_T2.txt", header=TRUE) 
data <- read.table("FAR_T3.txt", header=TRUE) 

library (lme4)

### Presence-absence models
gm2 <- glmer (FARb.w ~ treatment + (1|swale), family = binomial, data=data,
              control = glmerControl(optimizer="bobyqa"))

gm1 <- glm (FARb.w ~ treatment, family = binomial, data=data)

## Error-solving

https://rstudio-pubs-static.s3.amazonaws.com/33653_57fc7b8e5d484c909b615d8633c01d51.html

# Checking for singularity
tt <- getME(gm1,"theta")
ll <- getME(gm1,"lower")
min(tt[ll==0])


summary(gm1)$coef

## Camera frequency models
gm2 <- glmer (cbind (Camf.d, Days-Camf.d) ~ treatment + (1|swale)+ (1|obs), family = binomial, data=data,
              control = glmerControl(optimizer="bobyqa"))

gm3 <- glmer (cbind (Camf.y, Days-Camf.y) ~ treatment + (1|swale)+ (1|obs), 
              family = binomial (link = logit), data=data, control = glmerControl(optimizer="bobyqa"))

gm4 <- glmer (cbind (Camf.y, Days-Camf.y) ~ treatment + (1|obs), 
              family = binomial (link = logit), data=data, control = glmerControl(optimizer="bobyqa"))

## FAR frequency models
gm1 <- glmer (cbind (FAR.d, FAR.days-FAR.d) ~ treatment + (1|swale)+ (1|obs), family = binomial, data=data,
              control = glmerControl(optimizer="bobyqa"))
data$obs <- 1:nrow(data)
+ (1|obs)

#### Peter

fixef(gm1)[2:4] <- c(-1.386294, -1.386294, -1.386294) # small
fixef(gm1)[2:4] <- c(-0.4054651, -0.4054651, -0.4054651) # medium
fixef(gm1)[2:4] <- c(0.4054651, 0.4054651, 0.4054651) # large

powerSim(gm1, test=fixed("treatment", "lr"), nsim=50)



#### START HERE #####

# The thing to do is:
  
# Specify effect sizes for three of the treatment levels in relation to the reference level 
# Simulate data as desired, say 10 samples per treatment
# Perform a likelihood ratio test
# Repeat 1000 times to give a single power estimate with confidence intervals.
# The order of the levels must matter but at least the LR gives an overall number.

# s.tr.d = small effect size, plots per treatment, deer
# m.tr.d = medium effect size, plots per treatment vary, deer
# l.tr.d = large effect size, plots per treatment vary, deer


##### DEER #####

gm1 <- glmer (cbind (FAR.d, FAR.days-FAR.d) ~ treatment + (1|swale)+ (1|obs), family = binomial, data=data,
              control = glmerControl(optimizer="bobyqa"))

gm3 <- glmer (FARb.d ~ treatment + (1|swale), family = binomial, data=data,
              control = glmerControl(optimizer="bobyqa"))
gm4 <- glm (FARb.d ~ treatment, family = binomial, data=data)

#### SMALL

fixef(gm3)[2:4] <- c(-1.386294, -1.386294, -1.386294) # small
coef(gm4)[2:4] <- c(-1.386294, -1.386294, -1.386294) 

#### Plots per treatment
e1 <- extend (gm4, within="treatment", n = 160) 
pc1 <- powerCurve (e1, within="treatment", breaks = c(5, 10, 20, 40, 80, 160)) #1431
pc1 <- powerCurve (e1, within="treatment", breaks = c(5, 160)) 
s.tr.d4 <- summary(pc1)

#### Number of swales
e2 <- extend (gm1, along="swale", n=64) 
pc2 <- powerCurve(e2, along="swale", breaks = c(4, 8, 16, 32, 64))
s.sw.d <- summary(pc2)

#### MEDIUM

fixef(gm1)[2:4] <- c(-0.4054651, -0.4054651, -0.4054651) # medium

#### Plots per treatment
e1 <- extend (gm1, within="treatment", n = 160) 
pc1 <- powerCurve (e1, within="treatment", breaks = c(5, 10, 20, 40, 80, 160))
m.tr.d <- summary(pc1)

#### Number of swales
e2 <- extend (gm1, along="swale", n=64) 
pc2 <- powerCurve(e2, along="swale", breaks = c(4, 8, 16, 32, 64))
m.sw.d <- summary(pc2)

#### LARGE

fixef(gm1)[2:4] <- c(0.4054651, 0.4054651, 0.4054651) # large

#### Plots per treatment
e1 <- extend (gm1, within="treatment", n = 160) 
pc1 <- powerCurve (e1, within="treatment", breaks = c(5, 10, 20, 40, 80, 160))
l.tr.d <- summary(pc1)

#### Number of swales
e2 <- extend (gm1, along="swale", n=64) 
pc2 <- powerCurve(e2, along="swale", breaks = c(4, 8, 16, 32, 64))
l.sw.d <- summary(pc2)



##### KANGAROO #####

gm1 <- glmer (cbind (FAR.k, FAR.days-FAR.k) ~ treatment + (1|swale)+ (1|obs), family = binomial, data=data,
              control = glmerControl(optimizer="bobyqa"))

#### SMALL

fixef(gm1)[2:4] <- c(-1.386294, -1.386294, -1.386294) # small

#### Plots per treatment
e1 <- extend (gm1, within="treatment", n = 160) 
pc1 <- powerCurve (e1, within="treatment", breaks = c(5, 10, 20, 40, 80, 160))
s.tr.k <- summary(pc1)

#### Number of swales
e2 <- extend (gm1, along="swale", n=64) 
pc2 <- powerCurve(e2, along="swale", breaks = c(4, 8, 16, 32, 64))
s.sw.k <- summary(pc2)

#### MEDIUM

fixef(gm1)[2:4] <- c(-0.4054651, -0.4054651, -0.4054651) # medium

#### Plots per treatment
e1 <- extend (gm1, within="treatment", n = 160) 
pc1 <- powerCurve (e1, within="treatment", breaks = c(5, 10, 20, 40, 80, 160))
m.tr.k <- summary(pc1)

#### Number of swales
e2 <- extend (gm1, along="swale", n=64) 
pc2 <- powerCurve(e2, along="swale", breaks = c(4, 8, 16, 32, 64))
m.sw.k <- summary(pc2)

#### LARGE

fixef(gm1)[2:4] <- c(0.4054651, 0.4054651, 0.4054651) # large

#### Plots per treatment
e1 <- extend (gm1, within="treatment", n = 160) 
pc1 <- powerCurve (e1, within="treatment", breaks = c(5, 10, 20, 40, 80, 160))
l.tr.k <- summary(pc1)

#### Number of swales
e2 <- extend (gm1, along="swale", n=64) 
pc2 <- powerCurve(e2, along="swale", breaks = c(4, 8, 16, 32, 64))
l.sw.k <- summary(pc2)


##### RABBIT #####

gm1 <- glmer (cbind (FAR.r, FAR.days-FAR.r) ~ treatment + (1|swale)+ (1|obs), family = binomial, data=data,
              control = glmerControl(optimizer="bobyqa"))

#### SMALL

fixef(gm1)[2:4] <- c(-1.386294, -1.386294, -1.386294) # small

#### Plots per treatment
e1 <- extend (gm1, within="treatment", n = 160) 
pc1 <- powerCurve (e1, within="treatment", breaks = c(5, 10, 20, 40, 80, 160))
s.tr.r <- summary(pc1)

#### Number of swales
e2 <- extend (gm1, along="swale", n=64) 
pc2 <- powerCurve(e2, along="swale", breaks = c(4, 8, 16, 32, 64))
s.sw.r <- summary(pc2)

#### MEDIUM

fixef(gm1)[2:4] <- c(-0.4054651, -0.4054651, -0.4054651) # medium

#### Plots per treatment
e1 <- extend (gm1, within="treatment", n = 160) 
pc1 <- powerCurve (e1, within="treatment", breaks = c(5, 10, 20, 40, 80, 160))
m.tr.r <- summary(pc1)

#### Number of swales
e2 <- extend (gm1, along="swale", n=64) 
pc2 <- powerCurve(e2, along="swale", breaks = c(4, 8, 16, 32, 64))
m.sw.r <- summary(pc2)

#### LARGE

fixef(gm1)[2:4] <- c(0.4054651, 0.4054651, 0.4054651) # large

#### Plots per treatment
e1 <- extend (gm1, within="treatment", n = 160) 
pc1 <- powerCurve (e1, within="treatment", breaks = c(5, 10, 20, 40, 80, 160))
l.tr.r <- summary(pc1)

#### Number of swales
e2 <- extend (gm1, along="swale", n=64) 
pc2 <- powerCurve(e2, along="swale", breaks = c(4, 8, 16, 32, 64))
l.sw.r <- summary(pc2)


##### WOMBAT #####

# try running it with time included in the model
gm1 <- glm (FARb.w ~ treatment,# + (1|swale), 
              family = binomial, data=data)
              control = glmerControl(optimizer="bobyqa"))
gm2 <- glmer (FARb.w ~ treatment + (1|swale), 
            family = binomial, data=data, control = glmerControl(optimizer="bobyqa"))
tt <- getME(gm1,"theta")
ll <- getME(gm1,"lower")
min(tt[ll==0]) # 0.93 -indicates that random effect structure should be simplified
# proceed with GLM but replace fixef with coefficients

#### SMALL

fixef(gm2)[2:4] <- c(-1.386294, -1.386294, -1.386294) # small
coef(gm1)[2:4] <- c(-1.386294, -1.386294, -1.386294) # small

#### Plots per treatment
e1 <- extend (gm1, within="treatment", n = 160) 
# pc1 <- powerCurve (e1, within="treatment", breaks = c(5, 10, 20, 40, 80, 160))
pc1 <- powerCurve (e1, within="treatment", breaks = c(5, 20, 160))
s.tr.w5 <- summary(pc1)

#### Number of swales
e2 <- extend (gm1, along="swale", n=64) 
#pc2 <- powerCurve(e2, along="swale", breaks = c(4, 8, 16, 32, 64))
pc2 <- powerCurve(e2, along="swale", breaks = c(4, 16, 64))
s.sw.w4 <- summary(pc2)

#### MEDIUM

fixef(gm1)[2:4] <- c(-0.4054651, -0.4054651, -0.4054651) # medium

#### Plots per treatment
e1 <- extend (gm1, within="treatment", n = 160) 
pc1 <- powerCurve (e1, within="treatment", breaks = c(5, 10, 20, 40, 80, 160))
m.tr.w <- summary(pc1)

#### Number of swales
e2 <- extend (gm1, along="swale", n=64) 
pc2 <- powerCurve(e2, along="swale", breaks = c(4, 8, 16, 32, 64))
m.sw.w <- summary(pc2)

#### LARGE

fixef(gm1)[2:4] <- c(0.4054651, 0.4054651, 0.4054651) # large

#### Plots per treatment
e1 <- extend (gm1, within="treatment", n = 160) 
pc1 <- powerCurve (e1, within="treatment", breaks = c(5, 10, 20, 40, 80, 160))
l.tr.w <- summary(pc1)

#### Number of swales
e2 <- extend (gm1, along="swale", n=64) 
pc2 <- powerCurve(e2, along="swale", breaks = c(4, 8, 16, 32, 64))
l.sw.w <- summary(pc2)


##### WALLABY #####

gm1 <- glmer (cbind (FAR.y, FAR.days-FAR.y) ~ treatment + (1|swale)+ (1|obs), family = binomial, data=data,
              control = glmerControl(optimizer="bobyqa"))

#### SMALL

fixef(gm1)[2:4] <- c(-1.386294, -1.386294, -1.386294) # small

#### Plots per treatment
e1 <- extend (gm1, within="treatment", n = 160) 
pc1 <- powerCurve (e1, within="treatment", breaks = c(5, 10, 20, 40, 80, 160))
s.tr.y <- summary(pc1)

#### Number of swales
e2 <- extend (gm1, along="swale", n=64) 
pc2 <- powerCurve(e2, along="swale", breaks = c(4, 8, 16, 32, 64))
s.sw.y <- summary(pc2)

#### MEDIUM

fixef(gm1)[2:4] <- c(-0.4054651, -0.4054651, -0.4054651) # medium

#### Plots per treatment
e1 <- extend (gm1, within="treatment", n = 160) 
pc1 <- powerCurve (e1, within="treatment", breaks = c(5, 10, 20, 40, 80, 160))
m.tr.y <- summary(pc1)

#### Number of swales
e2 <- extend (gm1, along="swale", n=64) 
pc2 <- powerCurve(e2, along="swale", breaks = c(4, 8, 16, 32, 64))
m.sw.y <- summary(pc2)

#### LARGE

fixef(gm1)[2:4] <- c(0.4054651, 0.4054651, 0.4054651) # large

#### Plots per treatment
e1 <- extend (gm1, within="treatment", n = 160) 
pc1 <- powerCurve (e1, within="treatment", breaks = c(5, 10, 20, 40, 80, 160))
l.tr.y <- summary(pc1)

#### Number of swales
e2 <- extend (gm1, along="swale", n=64) 
pc2 <- powerCurve(e2, along="swale", breaks = c(4, 8, 16, 32, 64))
l.sw.y <- summary(pc2)



#### STICK TOGETHER ####

all <- rbind (s.tr.d, s.sw.d,
              s.tr.k, s.sw.k,
              s.tr.r, s.sw.r,
              s.tr.w, s.sw.w,
              s.tr.y, s.sw.y)
write.table (all, "farf_T3.txt", row.names = FALSE)

