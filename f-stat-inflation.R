library("broom")
library("lmtest")
library("sandwich")
source("funs.R")

n_obs <- 1000
n_sim <- 1000
bx <- .145
bu <- .09
maf <- .25

tv <- 2*(1-maf)*maf*bx^2 + 1*bu^2 + 1
ev <- 2*(1-maf)*maf*bx^2
r2 <- ev/tv
f <- get_f(r2, n_obs, 1)

results <- data.frame()
for (i in 1:n_sim){
    x <- rbinom(n_obs, 2, maf)
    u <- rnorm(n_obs)
    y <- x*bx + u*bu + rnorm(n_obs)
    fit <- lm(y ~ x)
    r2_est <- summary(fit)$r.squared
    f_est <- summary(fit)$fstatistic[1] %>% as.numeric
    results <- rbind(results, data.frame(r2_est, f_est, varY=var(y)))
}

results <- data.frame()
for (phi in seq(0, 0)){
    theta <- phi * bx
    for (i in 1:n_sim){
        x <- rbinom(n_obs, 2, .25)
        u <- rnorm(n_obs)
        y <- x*bx + u*bu + x*u*theta + rnorm(n_obs)
        fitx <- lm(y ~ x)
        r2 <- get_r2(
            fitx %>% tidy %>% dplyr::filter(term == "x") %>% dplyr::pull(estimate), 
            fitx %>% tidy %>% dplyr::filter(term == "x") %>% dplyr::pull(std.error), 
            .25, 
            n_obs
        )
        fx <- get_f(r2, n, n_snps)
        hc <- coeftest(fitx, vcov = vcovHC(fitx, type = "HC0"))
        fx_hc <- summary(hc)$fstatistic[1] %>% as.numeric
        results <- rbind(results, data.frame(phi, fx, fx_hc))
    }
}

# power
binom.test(sum(results$px<0.05), n_sim)
binom.test(sum(results$pu<0.05), n_sim)