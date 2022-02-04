library("dplyr")
library("broom")
library("varGWASR")
library("ggplot2")
source("funs.R")
set.seed(123)

n_sim <- 200
n_obs <- 1000/5

# simulate MR estimate with IV-exp interaction effect
results <- data.frame()
for (i in 1:n_sim){
    data <- data.frame()
    for (b in seq(0.2, 1, 0.2)){
        c <- rnorm(n_obs)
        z <- get_simulated_genotypes(0.25, n_obs)
        u <- rnorm(n_obs)
        x <- z + c + rnorm(n_obs)
        y <- x*b + c + rnorm(n_obs)
        data <- rbind(data, data.frame(
            c, z, u, x, y, b
        ))
    }
    fit <- varGWASR::model(data.frame(z,x), "z", "x")
    results <- rbind(results, fit)
}