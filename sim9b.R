library("dplyr")
library("ggplot2")
library("grid")
library("gtable")
library("TwoSampleMR")
library("viridis")
source("funs.R")
set.seed(123)

n_sim <- 50

r2_u <- 0.05 # U-X main effect variance explained
or_u <- 0 # U-Y main effect OR
or_x <- 0.05 # X-Y main effect OR
or_xu <- or_x * 0 # X-Y interaction effect half the size of the main effect
lor_u <- log(or_u)
lor_x <- log(or_x)
lor_xu <- log(or_xu)

results <- data.frame()
for (n_obs in 10000){
    for (r2_z in 0.05){ # variance explained by main effect of Z-X
        for (phi in 2.5){ # size of Z-X interaction effect relative to main effect
            r2_zu <- r2_z * phi
            z_b <- sqrt(r2_z)
            u_b <- sqrt(r2_u)
            zu_b <- sqrt(r2_zu)
            for (i in 1:n_sim){
                # simulate variables
                z <- get_simulated_genotypes(0.25, n_obs); z <- scale(z)
                u <- rbinom(n_obs, 1, 0.5); u <- scale(u)
                x <- z*z_b + u*u_b + z*u*zu_b + rnorm(n_obs, sd=sqrt(1-(r2_z+r2_u+r2_zu))); x <- scale(x)
                a <- -2 + x*lor_x + u*lor_u + x*u*lor_xu
                pr <- 1/(1+exp(-a))
                y <- rbinom(n_obs, 1, pr)

                # MR effect
                iv_exp <- lm(x~z)
                iv_out <- glm(y~z, family="binomial")
                b_exp <- iv_exp %>% tidy %>% dplyr::filter(term=="z") %>% dplyr::pull(estimate)
                se_exp <- iv_exp %>% tidy %>% dplyr::filter(term=="z") %>% dplyr::pull(std.error)
                b_out <- iv_out %>% tidy %>% dplyr::filter(term=="z") %>% dplyr::pull(estimate)
                se_out <- iv_out %>% tidy %>% dplyr::filter(term=="z") %>% dplyr::pull(std.error)
                wald <- mr_wald_ratio(b_exp, b_out, se_exp, se_out, NULL)

                # IV-exp variance effect
                result <- model(data.frame(z, x), "z", "x")

                # store result
                result$n_obs <- n_obs
                result$r2_z <- r2_z
                result$r2_zu <- r2_zu
                result$lor_xu <- lor_xu
                result$b_mr <- wald$b
                result$b <- or_x
                result$phi <- phi
                result$bias <- result$b_mr / result$b

                results <- rbind(results, result)
            }
        }
    }
}