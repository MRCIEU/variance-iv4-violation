library("broom")
library("dplyr")
library("TwoSampleMR")
set.seed(1234)

n_obs <- 1000
n_sim <- 1000
b <- .5

bp <- function(x, y){
    fit1 <- lm(y ~ x)
    dsq <- resid(fit1)^2
    xsq <- x^2
    fit2 <- lm(dsq ~ x + xsq)
    fit3 <- lm(dsq ~ 1)
    f <- anova(fit2, fit3)
    p <- tidy(f) %>% pull(p.value) %>% nth(2)
    return(p)
}

results <- data.frame()
for (i in 1:n_sim){
    # simulate effects
    b_exp <- .5
    b_out <- b * b_exp
    z <- rbinom(n_obs, 2, .4)
    w <- rnorm(n_obs)
    u <- rnorm(n_obs)
    x <- z*b_exp + w*.2 + z*u*b_exp + rnorm(n_obs)
    y <- z*b_out + w*.4 + z*u*b_out + rnorm(n_obs)

    # estimate SNP effects
    exp_fit <- lm(x ~ z)
    out_fit <- lm(y ~ z)

    # check IV1
    exp_p <- exp_fit %>% tidy %>% pull(p.value) %>% nth(2)
    stopifnot(exp_p < 5e-8)

    # estimate MR effect
    mr <- mr_wald_ratio(
        exp_fit %>% tidy %>% pull(estimate) %>% nth(2), 
        out_fit %>% tidy %>% pull(estimate) %>% nth(2),
        exp_fit %>% tidy %>% pull(std.error) %>% nth(2), 
        out_fit %>% tidy %>% pull(std.error) %>% nth(2)
    )

    results <- rbind(results, as.data.frame(mr))
}