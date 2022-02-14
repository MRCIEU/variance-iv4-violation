library("dplyr")
library("TwoSampleMR")
source("funs.R")
set.seed(123)

n_sim <- 1000
n_obs <- 10000
n_snps <- 9
phi <- 0

for (i in 1:n_sim){
    # confounder
    c <- rnorm(n_obs)
    # modifier
    u <- rbinom(n_obs, 1, 0.5)
    # Z1
    z1_q <- runif(1, min=0.05, max=0.5)
    z1 <- get_simulated_genotypes(z1_q, n_obs)
    # Zn
    zn_q <- runif(n_snps, min=0.05, max=0.5)
    zn <- sapply(zn_q, function(q) get_simulated_genotypes(q, n_obs))
    # interaction effect size relative to main effect
    bx_zu <- b[1]*phi_zx
    # simulate exposure
    x <- z1*b[1] + rowSums(t(t(zn)*b[-1])) + z1*u*bx_zu + u*bx_u + c + rnorm(n_obs, sd=sqrt(evx(z1*b[1], u*bx_u, z1*u*bx_zu)))
    var_x[i] <- var(x)
    evar_x[i] <- as.numeric(pwr$sd_y[1])^2
    # power
    p <- lm(x ~ z1 + u + zn) %>% tidy %>% dplyr::pull(p.value)
    # store test P
    res_p <- rbind(res_p, as.data.frame(t(p)))
    # store test F
    f <- sapply(1:(n_snps+1), function(n) if (n==1) {summary(lm(x~z1))$fstatistic[1] %>% as.numeric} else {summary(lm(x~zn[,n-1]))$fstatistic[1] %>% as.numeric})
    res_f <- rbind(res_f, as.data.frame(t(f)))
    # causal effect of X on Y explaining 5% total variance
    #y <- x + u + x*u*phi_xy + rnorm(n_obs, sd=sqrt())
}

# Z-X
## calculate power
pwr <- do.call(rbind.data.frame, apply(res_p, 2, function(x) binom.test(sum(x < 5e-8), n_sim) %>% tidy)) %>% dplyr::select(estimate, conf.low, conf.high)
pwr$term <- lm(x ~ z1 + u + zn) %>% tidy %>% dplyr::pull(term)
## mean F-stat
f_stat <- do.call(rbind.data.frame, apply(res_f, 2, function(x) t.test(x) %>% tidy)) %>% dplyr::select(estimate, conf.low, conf.high)
## proportion of F < 10
f_stat_lt10 <- apply(res_f, 2, function(x) sum(x < 10))

# X-Y
