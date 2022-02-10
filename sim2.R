library("dplyr")
library("TwoSampleMR")
source("funs.R")
set.seed(123)

n_sim <- 1000
n_obs <- 1000
n_snps <- 14

results <- data.frame()
for (i in 1:n_sim){
    # confounder
    c <- rnorm(n_obs)
    # modifier
    u <- rbinom(n_obs, 1, 0.5)
    # MAF of Z1
    z1_q <- runif(1, min=0.05, max=0.5)
    # MAF of Zn
    zn_q <- runif(n_snps, min=0.05, max=0.5)
    # Z1
    z1 <- get_simulated_genotypes(z1_q, n_obs)
    # Zn
    zn <- apply(zn_q, 2, function(q) get_simulated_genotypes(q, n_obs))

    b <- c(runif(n_snps * 0.5, min=-0.05, max=-0.01), runif(n_snps * 0.5, min=0.01, max=0.05))
    x <- z1*0.2 + u*0.2 + z1*u*0.2 + rowSums(t(t(zn)*b)) + c + rnorm(n_obs)
    lm(x~zn + z1*u) %>% summary
    y <- x + x*u + c + rnorm(n_obs)
}