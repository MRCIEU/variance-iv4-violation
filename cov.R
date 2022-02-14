source("funs.R")
set.seed(13)

n_sim <- 1000
n_obs <- 1000
q <- 0.25
mu_u <- 5

results <- data.frame()
for (i in 1:n_sim){
    x <- get_simulated_genotypes(q, n_obs)
    u <- rnorm(n_obs, mean=mu_u)
    xu <- x*u
    y <- x + u + xu + rnorm(n_obs)
    results <- rbind(results, data.frame(
        evar_y=var(x) + var(u) + var(xu) + 2*cov(x, xu) + 2*cov(u, xu) + 1,
        ecov_xxu=mean(x*xu) - q*2 * (q*2 * mu_u),
        cov_xxu=cov(x, xu),
        var_y=var(y)
    ))
}