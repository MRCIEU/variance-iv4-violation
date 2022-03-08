library("dplyr")
library("ggplot2")
library("grid")
library("gtable")
library("TwoSampleMR")
library("viridis")
source("funs.R")
set.seed(123)

n_sim <- 50
n_obs <- 10000
r2_z <- 0.05; z_b <- sqrt(r2_z)
r2_u <- 0.05; u_b <- sqrt(r2_u)
r2_zu <- 0; zu_b <- sqrt(r2_zu)
lor_x <- log(0.05)
lor_u <- log(0.05)
lor_xu <- log(0.025)

results <- data.frame()
for (i in 1:n_sim){
    # simulate variables
    z <- get_simulated_genotypes(0.25, n_obs); z <- scale(z)
    u <- rbinom(n_obs, 1, 0.5); u <- scale(u)
    x <- z*z_b + u*u_b + z*u*zu_b + rnorm(n_obs, sd=sqrt(1-(r2_z+r2_u+r2_zu))); x <- scale(x)
    a <- -2 + x*lor_x + u*0 + x*u*0
    pr <- 1/(1+exp(-a))
    y <- rbinom(n_obs, 1, pr)

    r2 <- summary(lm(x ~ z))$r.squared

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
    result$b_mr <- wald$b
    result$b <- lor_x
    result$phi <- phi
    result$bias <- result$b_mr / result$b
    result$r2 <- r2

    results <- rbind(results, result)
}