library("broom")
set.seed(12)

# simulation to calculate variance explained

n_sim <- 200
n_obs <- 1000

results <- data.frame()
for (sd_x in seq(1, 10, 1)){
    r2 <- rep(NA, n_sim)
    for (i in 1:n_sim) {
        u <- rbinom(n_obs, 1, 0.5)
        x <- rnorm(n_obs, sd=sd_x)
        expl <- sd_x^2+0.5^2
        y <- x + u + rnorm(n_obs, sd=sqrt(expl/.2-expl))
        r2[i] <- summary(lm(y~x+u))$r.squared
    }
    res <- t.test(r2) %>% tidy
    res$sd <- sd_x
    results <- rbind(res, results)
}