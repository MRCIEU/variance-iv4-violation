library("dplyr")
library("broom")
library("ggplot2")
library("TwoSampleMR")
source("funs.R")
set.seed(123)

n_sim <- 1000
n_obs <- 1000

phi <- 0

results <- data.frame()
for (i in 1:n_sim){
    z_b <- 0.1
    u_b <- 0.1
    z <- get_simulated_genotypes(0.25, n_obs)
    u <- rbinom(n_obs, 1, 0.5)
    x <- z*z_b + u*u_b + z*u*phi + rnorm(n_obs)
    results <- rbind(results, lm(x~z*u) %>% tidy)
}

# estimate power
binom.test(results %>% dplyr::filter(term == "z") %>% dplyr::filter(p.value < 0.05) %>% nrow, n_sim) %>% tidy
binom.test(results %>% dplyr::filter(term == "u") %>% dplyr::filter(p.value < 0.05) %>% nrow, n_sim) %>% tidy
binom.test(results %>% dplyr::filter(term == "z*u") %>% dplyr::filter(p.value < 0.05) %>% nrow, n_sim) %>% tidy