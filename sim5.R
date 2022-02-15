library("dplyr")
library("broom")
library("ggplot2")
library("TwoSampleMR")
source("funs.R")
set.seed(123)

n_sim <- 1000
n_obs <- 10000

results <- data.frame()
for (i in 1:n_sim){
    u_p <- 0.5
    q <- 0.25
    p <- 1-q
    z_r2 <- 0.05
    z <- get_simulated_genotypes(q, n_obs)
    z_b <- sqrt(z_r2/(2*p*q))
    u_r2 <- 0.05
    u <- rbinom(n_obs, 1, u_p)
    u_b <- sqrt(u_r2/0.5^2)
    zu_v <- var(z*u) + 2*cov(z, z*u) + 2*cov(u, z*u)
    zu_r2 <- 0.1
    zu_b <- sqrt(zu_r2/zu_v)
    x <- z*z_b + u*u_b + z*u*zu_b + rnorm(n_obs, sd=sqrt(1-(z_r2+u_r2+zu_r2)))
    r2 <- summary(lm(x ~ z*u))$r.squared
    ev_zu <- q*u_p^2+u_p^2*(2*p*q)+(2*p*q)*u_p^2
    ecov_zzu <- u_p^2
    results <- rbind(results, data.frame(r2, var_x=var(x), var_all=var(z*z_b + u*u_b + z*u*zu_b), var_int=var(z*z_b + u*u_b + z*u*zu_b), ev_zu, v_zu=var(z*u), ecov_zzu, cov_zzu=cov(z, z*u)))
}


results <- data.frame()
for (q in seq(0.05, 0.5, 0.05)){
    for (i in 1:n_sim){
        u_p <- 0.5
        z <- get_simulated_genotypes(q, n_obs)
        u <- rbinom(n_obs, 1, u_p)
        zu <- z*u
        ecov <- mean(z*zu) - (q*2)*(q*2*u_p)
        results <- rbind(results, data.frame(cov=cov(z, zu), q, ecov))
    }
}

results <- data.frame()
r2_u <- 0.05
u_b <- sqrt(r2_u)
r2_z <- 0.05
z_b <- sqrt(r2_z)
r2_zu <- 0.1
zu_b <- sqrt(r2_zu)
for (i in 1:n_sim){
    z <- get_simulated_genotypes(0.25, n_obs)
    u <- rbinom(n_obs, 1, 0.5)
    z <- scale(z)
    u <- scale(u)
    x <- z*z_b + u*u_b + z*u*zu_b + rnorm(n_obs, sd=sqrt(1-(r2_u+r2_z+r2_zu)))
    r <- summary(lm(x ~ z*u))$r.squared
    results <- rbind(results, data.frame(var_x=var(x), r))
}