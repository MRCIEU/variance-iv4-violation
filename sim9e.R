library("dplyr")
library("broom")
source("funs.R")
set.seed(234)

n_sim <- 200
n_obs <- 10000

# log OR effects
b0 <- log(0.5)
b1 <- log(5) # causal effect of X-Y
b2 <- log(3) # causal effect of U-Y

results <- data.frame()
for (i in 1:n_sim){

    # simulate exposures 
    z <- get_simulated_genotypes(0.3, n_obs); z <- scale(z)
    x <- z*sqrt(0.05) + rnorm(n_obs, sd=sqrt(0.95)); x <- scale(x)

    # simulate outcome
    p <- b0 + x*b1
    p <- 1/(1+exp(-p))
    y <- rbinom(n_obs, 1, p)
    
    # IV estimates
    b_exp <- lm(x ~ z) %>% tidy %>% dplyr::filter(term == "z") %>% dplyr::pull("estimate")
    b_out <- glm(y ~ z, family="binomial") %>% tidy %>% dplyr::filter(term == "z") %>% dplyr::pull("estimate")

    # obs
    b_obs <- glm(y ~ x, family="binomial") %>% tidy %>% dplyr::filter(term == "x") %>% dplyr::pull("estimate")

    # calculate Wald estimate
    b_mr <- b_out / b_exp
    
    # store results
    results <- rbind(results, data.frame(
        b_mr, b_obs
    ))
}

# check estimates
t.test(results$b_mr, mu=b1)
t.test(results$b_obs, mu=b1)