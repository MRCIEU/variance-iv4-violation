# sample size (X) vs power (Y) to detect increasing fold-bias (color)
# facet-grid with column and rows Z-X and X-Y variance explained
# two series - continuous vs binary outcomes

library("dplyr")
library("ggplot2")
library("TwoSampleMR")
library("viridis")
source("funs.R")
set.seed(123)

n_sim <- 50
r2_u <- 0.05
r2_x <- 0.05
r2_xu <- r2_x * 0.05

results <- data.frame()
for (n_obs in c(100000)){
    for (r2_z in seq(0.01, 0.05, 0.01)){
        for (phi in seq(0, 6, 0.5)){
            r2_zu <- r2_z * phi
            u_b <- sqrt(r2_u)
            z_b <- sqrt(r2_z)
            zu_b <- sqrt(r2_zu)
            xu_b <- sqrt(r2_xu)
            x_b <- sqrt(r2_x)
            for (i in 1:n_sim){
                # simulate variables
                z <- get_simulated_genotypes(0.25, n_obs); z <- scale(z)
                u <- rbinom(n_obs, 1, 0.5); u <- scale(u)
                x <- z*z_b + u*u_b + z*u*zu_b + rnorm(n_obs, sd=sqrt(1-(r2_z+r2_u+r2_zu))); x <- scale(x)
                y <- x*x_b + u*u_b + x*u*xu_b + rnorm(n_obs, sd=sqrt(1-(r2_x+r2_u+r2_xu))); y <- scale(y)

                # MR effect
                iv_exp <- lm(x~z)
                iv_out <- lm(y~z)
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
                result$r2_xu <- r2_xu
                result$b_mr <- wald$b
                result$b <- sqrt(r2_x)
                result$phi <- phi
                result$bias <- (result$b_mr - result$b) / result$b

                results <- rbind(results, result)
            }
        }
    }
}

# power
results %>% dplyr::group_by(n_obs, r2_z, phi) %>% dplyr::summarize(binom.test(sum(phi_p < 0.05), n()) %>% tidy)