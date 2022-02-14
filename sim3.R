library("dplyr")
library("TwoSampleMR")
library("varGWASR")
source("funs.R")
set.seed(123)

#' Get the effect size for each SNP to have R^2 assuming standard Normal outcome
es <- function(maf){
    b <- sapply(maf, function(q) {p <- 1-q; v <- 2*p*q; return(sqrt(r2_z/v))})
    return(b)
}

n_sim <- 50
n_obs <- 10000
n_snps <- 10
r2_z <- 0.01 # explained variance of main effect of Z on X
r2_u <- 0.05 # explained variance of main effect of U on X
r2_c <- 0.05 # explained variance of main effect of C on X
r2_x <- 0.05 # explained variance of main effect of X on Y

results <- data.frame()
for (r2_zu in seq(0, 0.1, 0.02)){
    for (r2_xu in seq(0, 0, 0)){
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
            x_b <- sqrt(r2_x/1)
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
            # IV-exp
            b_exp <- sapply(1:(n_snps-1), function(n) lm(x ~ zn[,n]) %>% tidy %>% dplyr::filter(term == "zn[, n]") %>% dplyr::pull(estimate))
            b_exp <- c(lm(x ~ z1) %>% tidy %>% dplyr::filter(term == "z1") %>% dplyr::pull(estimate), b_exp)
            se_exp <- sapply(1:(n_snps-1), function(n) lm(x ~ zn[,n]) %>% tidy %>% dplyr::filter(term == "zn[, n]") %>% dplyr::pull(std.error))
            se_exp <- c(lm(x ~ z1) %>% tidy %>% dplyr::filter(term == "z1") %>% dplyr::pull(std.error), se_exp)
            # IV-out
            b_out <- sapply(1:(n_snps-1), function(n) lm(y ~ zn[,n]) %>% tidy %>% dplyr::filter(term == "zn[, n]") %>% dplyr::pull(estimate))
            b_out <- c(lm(y ~ z1) %>% tidy %>% dplyr::filter(term == "z1") %>% dplyr::pull(estimate), b_out)
            se_out <- sapply(1:(n_snps-1), function(n) lm(y ~ zn[,n]) %>% tidy %>% dplyr::filter(term == "zn[, n]") %>% dplyr::pull(std.error))
            se_out <- c(lm(y ~ z1) %>% tidy %>% dplyr::filter(term == "z1") %>% dplyr::pull(std.error), se_out)
            # IV variance effects
            result1 <- varGWASR::model(data.frame(z1,x), "z1", "x")
            names(result1) <- paste0(names(result1), ".x")
            result2 <- varGWASR::model(data.frame(z1,y), "z1", "y")
            names(result2) <- paste0(names(result2), ".y")
            result <- cbind(result1, result2)
            # observational estimates
            b_obs <- lm(y ~ x) %>% tidy %>% dplyr::filter(term == "x") %>% dplyr::pull(estimate)
            b_obs_a <- lm(y ~ x + c + u) %>% tidy %>% dplyr::filter(term == "x") %>% dplyr::pull(estimate)
            # MR estimates
            ivw <- TwoSampleMR::mr_ivw(b_exp, b_out, se_exp, se_out)
            # store results
            result$b_obs <- b_obs
            result$b_obs_a <- b_obs_a
            result$b_mr <- ivw$b
            result$r2_zu <- r2_zu
            result$r2_xu <- r2_xu
            result$x_b <- x_b

            results <- rbind(results, result)
        }
    }
}