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

n_sim <- 300
n_obs <- 50000
n_snps <- 10
r2_z <- 0.05 # per-SNP explained variance
r2_u <- 0.05 # explained variance of main effect of U on X
r2_c <- 0.05 # explained variance of main effect of C on X

results <- data.frame()
var_x <- rep(NA, n_sim)
for (r2_zu in seq(0.1, 0.1, 0.1)){
    for (i in 1:n_sim){
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
        v_z1u <- var(z1*u) + 2*cov(z1, z1*u) + 2*cov(u, z1*u) # variance of interaction term
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
        var_x[i] <- var(x)
        results <- rbind(results, assoc("z1", z1, x))
        results <- rbind(results, assoc("z2", zn[,1], x))
        results <- rbind(results, assoc("z3", zn[,2], x))
        results <- rbind(results, assoc("z4", zn[,3], x))
        results <- rbind(results, assoc("z5", zn[,4], x))
        results <- rbind(results, assoc("z6", zn[,5], x))
        results <- rbind(results, assoc("z7", zn[,6], x))
        results <- rbind(results, assoc("z8", zn[,7], x))
        results <- rbind(results, assoc("z9", zn[,8], x))
        results <- rbind(results, assoc("z10", zn[,9], x))
        results <- rbind(results, assoc("u", u, x))
        results <- rbind(results, assoc("c", c, x))
        fit <- lm(x ~ z1*u) %>% summary
        results <- rbind(results, data.frame(term="z1*u", r.squared=fit$r.squared, fstatistic=fit$fstatistic[1] %>% as.numeric))
    }
}

# check parameters
results %>% dplyr::group_by(term) %>% dplyr::summarize(r.squared=t.test(r.squared) %>% tidy)
results %>% dplyr::group_by(term) %>% dplyr::summarize(r.squared=t.test(fstatistic) %>% tidy)