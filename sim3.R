library("dplyr")
library("TwoSampleMR")
source("funs.R")
set.seed(123)

#' Get the effect size for each SNP to have R^2 assuming standard Normal outcome
es <- function(maf){
    b <- sapply(maf, function(q) {p <- 1-q; v <- 2*p*q; return(sqrt(r2_z/v))})
    return(b)
}

#' Estimate R^2 and F-statistic
assoc <- function(term, x, y){
    fit <- summary(lm(y~x))
    return(data.frame(
        term,
        r.squared=fit$r.squared,
        fstatistic=fit$fstatistic[1] %>% as.numeric
    ))
}

n_sim <- 200
n_obs <- 10000
n_snps <- 10
r2_z <- 0.01 # per-SNP explained variance
r2_u <- 0.05 # explained variance of main effect of U on X
r2_c <- 0.05 # explained variance of main effect of C on X
x_b <- sqrt(0.05/1) # causal effect of X on Y

results <- data.frame()
for (r2_zu in seq(0.1, 0.1, 0)){
    for (r2_xu in seq(0.1, 0.1, 0)){
        for (i in 1:n_sim){
            # Z-X

            # Z1
            z1_q <- runif(1, min=0.05, max=0.5)
            z1 <- get_simulated_genotypes(z1_q, n_obs)
            # Zn
            zn_q <- runif(n_snps-1, min=0.05, max=0.5)
            zn <- sapply(zn_q, function(q) get_simulated_genotypes(q, n_obs))
            # modifier
            u <- rbinom(n_obs, 1, 0.5)
            # confounder
            c <- rnorm(n_obs)
            # betas
            z_b <- es(c(z1_q, zn_q))
            u_b <- sqrt(r2_u/0.5^2)
            c_b <- sqrt(r2_c/1)
            v_z1u <- var(z1*u) + 2*cov(z1, z1*u) + 2*cov(u, z1*u) # variance of interaction term # TODO expectation
            zu_b <- sqrt(r2_zu/v_z1u)
            # calculate explained variance
            ev <- r2_z*n_snps + r2_u + r2_c + r2_zu
            stopifnot(!ev > 1)
            # standard Normal exposure
            x <- 
                z1*z_b[1] + 
                rowSums(t(t(zn)*z_b[-1])) + 
                u*u_b +
                c*c_b +
                z1*u*zu_b +
                rnorm(n_obs, sd=sqrt(1 - ev))
            
            # X-Y

            # betas
            u_b <- sqrt(r2_u/0.5^2)
            c_b <- sqrt(r2_c/1)
            v_xu <- var(x*u) + 2*cov(x, x*u) + 2*cov(u, x*u) # variance of interaction term # TODO expectation
            xu_b <- sqrt(r2_xu/v_xu)
            # calculate explained variance
            ev <- r2_x + r2_u + r2_c + r2_xu
            stopifnot(!ev > 1)
            # standard Normal outcome
            y <- 
                x*x_b +
                u*u_b +
                c*c_b +
                x*u*xu_b +
                rnorm(n_obs, sd=sqrt(1 - ev))

            # MR
            

            # IV variance effects

        }
    }
}