# Structure of data (fixed and random effects)
tickdata <- expand.grid (chick = 1:3, brood = 1:2, location = 1:10)

# Work out the mean response (across all treatments)

# Work out the random effect variance associated with swale

# Then you can simulate response data

sim.tickdata <- function(...){
  sim.glmm(design.data = tickdata,
           fixed.eff = list(intercept = log(10)),
           rand.V = c(location = 1, brood = 0.7, chick = 0.3),
           distribution = "poisson")
}


## Calculate the margin of error using a Poisson GLM
## Calculate as half the 95% CI width as a percentage of the estimated mean tick count

sim.tick.err <- function(...){
  fit <-
    glmer(response ~ (1 | location) + (1 | brood) + (1 | chick), family = "poisson", data = sim.tickdata())
  intercept <- fixef(fit)
  estimate <- exp(intercept)
  se <- coef(summary(fit))[,"Std. Error"]
  ci.lo <- exp(intercept - qnorm(0.975) * se)
  ci.hi <- exp(intercept + qnorm(0.975) * se)
  as.vector(100 * 0.5 * (ci.hi - ci.lo) / estimate) # Output is labelled 'intercept' but it is the margin of error as a percentage
}

## Apply to 50 simulated datasets

sim.err <- sapply(1:50, sim.tick.err)
mean(sim.err) # gives margin of error as a percentage. 78% = an average 95% CI extends to +/-78% of the estimate. 
# This is (approx) the expected margin of error. 

# The margin of error has a wide distribution, so to be more cautious we could ask
#  what margin of error are we 80% confident of doing better than?
quantile(sim.err, 0.8)



