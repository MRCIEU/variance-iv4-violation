library("broom")
library("dplyr")
set.seed(1234)

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

n_obs <- 1000
n_sim <- 1000
delta <- .2
phi <- .5
theta <- delta*phi

p_delta <- rep(NA, n_sim)
p_var <- rep(NA, n_sim)
for (i in 1:n_sim){
    x <- rnorm(n_obs)
    u <- rbinom(n_obs, 1, .5)
    y <- delta*x + theta*x*u + rnorm(n_obs)
    fit <- lm(y ~ x)
    p_delta[i] <- fit %>% tidy %>% pull(p.value) %>% nth(2)
    p_var[i] <- bp(x,y)
}

binom.test(sum(p_delta < 5e-8), n_sim)
binom.test(sum(p_var < 0.05), n_sim)

# test for subgroup effect
#lm(y[u==0] ~ x[u==0]) %>% summary
#lm(y[u==1] ~ x[u==1]) %>% summary