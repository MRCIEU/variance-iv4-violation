library("dplyr")
library("ggplot2")
library("grid")
library("gtable")
library("TwoSampleMR")
library("viridis")
library("genpwr")
source("funs.R")
set.seed(123)

n_obs <- 1000
n_sim <- 1000

# sim9
results <- data.frame()
for (i in 1:n_sim){
    r2_u <- 0.2 # U-Y main effect (small Cohen d)
    r2_x <- 0.2 # X-Y main effect (small Cohen d)
    r2_xu <- r2_x * 0.5 # X-Y interaction effect half the size of the main effect
    r2_z <- 0.01
    z_b <- sqrt(r2_z)
    z <- get_simulated_genotypes(0.25, n_obs); z <- scale(z)
    x <- z*z_b + rnorm(n_obs, sd=sqrt(1-r2_z))
    r <- summary(lm(x ~ z))$r.squared
    results <- rbind(results, data.frame(r))
}

# power sim
results <- data.frame()
for (i in 1:n_sim){
    af <- 0.4
    delta <- as.numeric(genpwr.calc(calc = "es", model = "linear", ge.interaction = NULL, N = n_obs, k = NULL, MAF = af, Power = 0.8, Alpha = 0.05, sd_y = 1, True.Model = "Additive", Test.Model = "Additive")$ES_at_Alpha_0.05)
    phi <- 0.5
    theta <- delta * phi
    # simulate covariates
    data <- data.frame(
        X = get_simulated_genotypes(af, n_obs)
    )

    # simulate outcome
    data$Y <- data$X * delta+ rnorm(n_obs)
    r <- summary(lm(Y ~ X, data=data))$r.squared
    results <- rbind(results, data.frame(r, varY=var(data$Y)))
}