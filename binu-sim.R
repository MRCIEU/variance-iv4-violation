library("broom")
library("dplyr")
library("TwoSampleMR")
set.seed(1234)

n_obs <- 100000
n_sim <- 1000
b <- .5 # causal effect
phi <- .5 # size of interaction relative to main effect
delta <- .25 # IV-exp effect has 95% power alpha=5e-8
theta <- delta * phi # interaction effect size

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
p <- rep(NA, n_sim)
for (i in 1:n_sim){
    z <- rbinom(n_obs, 2, .4) # IV
    w <- rnorm(n_obs) # confounder
    u <- rbinom(n_obs, 1, .5) + 1 # modifier

    # subgroup treatment effects
    b1 <- theta * -0.5
    b2 <- theta

    # simulate IV-exp effect with interaction
    x <- b1 * z + b2 * z * u + delta * z + w*delta*.5 + rnorm(n_obs)

    # IV-out effect
    b_exp <- delta 
    b_out <- b * b_exp

    # simulate outcome
    y <- z*b_out + w*b_out*.5 + rnorm(n_obs)

    # estimate SNP effects
    exp_fit <- lm(x ~ z)
    out_fit <- lm(y ~ z)

    # check IV1
    exp_p <- exp_fit %>% tidy %>% pull(p.value) %>% nth(2)
    exp_r2 <- exp_fit %>% summary %>% .[["r.squared"]]
    p[i] <- exp_p

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
binom.test(sum(p<5e-8), n_sim)