library("dplyr")
library("TwoSampleMR")
library("genpwr")
source("funs.R")
set.seed(123)

n_sim <- 1000
n_obs <- 10000
n_snps <- 9
phi_zx <- 0
phi_xy <- 0

# Main effect of U on X
bx_u <- 0.18

# var(UCE)
vuce <- 0.25*bx_u^2 + 1 + 1

# Func to estimate residual variance of X
evx <- function(x, u, xu){
    return(1 - (var(xu) + 2*cov(x, xu) + 2*cov(u, xu)))
}

res_p <- data.frame()
res_f <- data.frame()
var_x <- rep(NA, n_sim)
evar_x <- rep(NA, n_sim)
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
    # calculate SNP effect sizes assuming 80% power with P < 5e-8
    pwr <- genpwr.calc(
        calc = "es", 
        model = "linear", 
        N = n_obs, 
        Power = 0.8,
        MAF = c(z1_q, zn_q),
        Alpha = 5e-8,
        sd_y = sqrt(vuce)
    )
    b <- pwr %>% 
        dplyr::filter(Test.Model=="Additive", True.Model=="Additive") %>% 
        dplyr::pull("ES_at_Alpha_5e-08") %>% 
        as.numeric
    stopifnot(!any(is.na(b)))
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
