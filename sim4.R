library("dplyr")
library("broom")
library("ggplot2")
library("TwoSampleMR")
source("funs.R")
set.seed(123)

n_sim <- 1000
n_obs <- 10000

results <- data.frame()
for (i in 1:n_sim){
    # simulate variables
    q <- 0.25
    p <- 1 - q
    z <- get_simulated_genotypes(q, n_obs)
    u <- rbinom(n_obs, 1, 0.5)
    b_z <- sqrt(0.05/(2*p*q))
    b_u <- sqrt(0.05/0.5^2)
    
    # simulate exposure
    x <- z*b_z + u*b_u + z*u + rnorm(n_obs)

    # simulate outcome
    y <- x + x*u + rnorm(n_obs)

    # check effect sizes
    results <- rbind(results, data.frame(
        var_x=var(x),
        var_y=var(y),
        zx_r2=summary(lm(x~z))$r.squared,
        ux_r2=summary(lm(x~u))$r.squared,
        zux_r2=summary(lm(x~z*u))$r.squared,
        cx_r2=summary(lm(x~c))$r.squared,
        xy_r2=summary(lm(y~x+c+u))$r.squared,
        uy_r2=summary(lm(y~u))$r.squared,
        uy_r2=summary(lm(y~u))$r.squared,
        xuy_r2=summary(lm(y~x*u))$r.squared,
        xuy_r2=summary(lm(y~x*u))$r.squared,
        cy_r2=summary(lm(y~c))$r.squared,
        r2_zu=r2_zu,
        r2_xu=r2_xu
    ))
}

# check values
results %>% dplyr::select(-r2_zu, -r2_xu) %>% apply(., 2, function(x) t.test(x) %>% tidy)