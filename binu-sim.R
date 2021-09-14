library("broom")
library("dplyr")
library("TwoSampleMR")
set.seed(1234)

n_obs <- 1000
n_sim <- 1000
b <- .5 # causal effect
phi <- 2

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
    b_exp <- .835 # IV-exp effect explains 5% variation in X
    b_out <- b * b_exp # IV-out effect
    z <- rbinom(n_obs, 2, .4) # IV
    w <- rnorm(n_obs) # confounder
    u <- rnorm(n_obs) # modifier
    x <- z*b_exp + w*b_exp*2 + z*u*b_exp*phi + rnorm(n_obs) # exposure
    y <- z*b_out + w*b_out*2 + z*u*b_out*phi + rnorm(n_obs) # outcome

    # estimate SNP effects
    exp_fit <- lm(x ~ z)
    out_fit <- lm(y ~ z)

    # check IV1
    exp_p <- exp_fit %>% tidy %>% pull(p.value) %>% nth(2)
    exp_r2 <- exp_fit %>% summary %>% .[["r.squared"]]

    if (exp_p > 5e-8){
        next
    }

    # estimate MR effect
    mr <- mr_wald_ratio(
        exp_fit %>% tidy %>% pull(estimate) %>% nth(2), 
        out_fit %>% tidy %>% pull(estimate) %>% nth(2),
        exp_fit %>% tidy %>% pull(std.error) %>% nth(2), 
        out_fit %>% tidy %>% pull(std.error) %>% nth(2)
    )

    res <- as.data.frame(mr)
    res$bp <- bp(z, x)
    res$exp_p <- exp_p
    res$exp_r2 <- exp_r2
    results <- rbind(results, res)
}

# check casual estimate is similar to effect
t.test(results$b, mu=b)
t.test(results$bp)
t.test(results$exp_r2)