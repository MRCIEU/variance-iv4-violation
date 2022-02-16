library("dplyr")
library("ggplot2")
library("TwoSampleMR")
library("viridis")
source("funs.R")
set.seed(123)

n_sim <- 50
n_obs <- 10000
r2_z <- 0.05
r2_x <- 0.05

results <- data.frame()
for (h in seq(-1, 1, .5)){
    for (s in seq(0, 1, 0.2)){
        for (i in 1:n_sim){
            # betas
            b0 <- 0.5
            b1 <- b0-h*s
            b2 <- b0-s

            # genotype
            z <- get_simulated_genotypes(0.25, n_obs)
            u0 <- z==0
            u1 <- z==1
            u2 <- z==2
            z_ev <- var(z*u0*b0 + z*u1*b1 + z*u2*b2)
            x_tv <- z_ev / r2_z

            # exposure           
            x <- z*u0*b0 + z*u1*b1 + z*u2*b2 + rnorm(n_obs, sd=sqrt(x_tv-z_ev))

            # outcome
            x_ev <- var(x)
            y_tv <- x_ev / r2_x
            y <- x + rnorm(n_obs, sd=sqrt(y_tv-x_ev))

            # estimate MR effect
            b_exp <- lm(x~z) %>% tidy %>% dplyr::filter(term=="z") %>% dplyr::pull(estimate)
            se_exp <- lm(x~z) %>% tidy %>% dplyr::filter(term=="z") %>% dplyr::pull(std.error)
            b_out <- lm(y~z) %>% tidy %>% dplyr::filter(term=="z") %>% dplyr::pull(estimate)
            se_out <- lm(y~z) %>% tidy %>% dplyr::filter(term=="z") %>% dplyr::pull(std.error)
            wald <- mr_wald_ratio(b_exp, b_out, se_exp, se_out, NULL)
            v_z <- model(data.frame(z, x), "z", "x")
            names(v_z) <- paste0(names(v_z), ".x")
            v_y <- model(data.frame(z, y), "z", "y")
            names(v_y) <- paste0(names(v_y), ".y")
            result <- cbind(v_z, v_y)
            result$b_mr <- wald$b
            result$se_mr <- wald$se_mr
            result$b <- sqrt(r2_x)
            result$h <- h
            result$s <- s

            results <- rbind(results, result)
        }
    }
}