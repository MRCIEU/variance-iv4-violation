library("dplyr")
library("broom")
set.seed(234)

n_sim <- 50
n_obs <- 10000
b <- log(5) # causal effect of X-Y in log OR

results <- data.frame()
for (i in 1:n_sim){
    z <- rbinom(n_obs, 2, .5) # instrument
    x <- z + rnorm(n_obs) # continuous exposure
    a <- -2 + x*b; pr <- 1/(1+exp(-a)); y <- rbinom(n_obs, 1, pr) # binary outcome
    
    # estimate IV-exposure effect
    iv_exp <- lm(x ~ z)
    b_exp <- iv_exp %>% tidy %>% dplyr::filter(term == "z") %>% dplyr::pull(estimate)
    f_exp <- summary(iv_exp)$fstatistic[1] %>% as.numeric
    # estimate IV-outcome effect
    b_out <- glm(y ~ z, family="binomial") %>% tidy %>% dplyr::filter(term == "z") %>% dplyr::pull(estimate)
    # estimate observational association
    b_obs <- glm(y ~ x, family="binomial") %>% tidy %>% dplyr::filter(term == "x") %>% dplyr::pull(estimate)

    # store results
    results <- rbind(results, data.frame(
        b_mr=b_out/b_exp, b_obs, f_exp
    ))
}

# check estmates
t.test(results$f_exp)
t.test(results$b_mr, mu=b)
t.test(results$b_obs, mu=b)