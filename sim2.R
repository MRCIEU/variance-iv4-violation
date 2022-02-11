library("dplyr")
library("TwoSampleMR")
library("genpwr")
source("funs.R")
set.seed(123)

n_sim <- 2000
n_obs <- 10000
n_snps <- 9
phi <- 0

# Main effect of U on X
bx_u <- 0.185

res_p <- data.frame()
res_f <- data.frame()
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
    # calculate SNP effect sizes
    pwr <- genpwr.calc(
        calc = "es", 
        model = "linear", 
        N = n_obs, 
        Power = 0.8,
        MAF = c(z1_q, zn_q), 
        Alpha = 5e-8, 
        sd_y = sqrt(
            (0.25*bx_u)^2 + 1 + 1 # main effect of U and C and E
        )
    )
    b <- pwr %>% dplyr::filter(Test.Model=="Additive", True.Model=="Additive") %>% dplyr::pull("ES_at_Alpha_5e-08") %>% as.numeric
    stopifnot(!any(is.na(b)))
    # interaction effect size relative to main effect
    bx_zu <- b[1]*phi
    # variance of Z1U
    vz1u <- z1_q^2*0.5^2 + 0.5^2*(2*(1-z1_q)*z1_q) + (2*(1-z1_q)*z1_q)*0.5^2
    # simulate exposure 
    x <- z1*b[1] + rowSums(t(t(zn)*b[-1])) + u*bx_u + z1*u*bx_zu + c + rnorm(n_obs, sd=sqrt(1-vz1u*bx_zu^2))
    # power
    p <- lm(x ~ z1 + u + zn) %>% tidy %>% dplyr::pull(p.value)
    # store test P
    res_p <- rbind(res_p, as.data.frame(t(p)))
    # store test F
    f <- sapply(1:(n_snps+1), function(n) if (n==1) {summary(lm(x~z1))$fstatistic[1] %>% as.numeric} else {summary(lm(x~zn[,n-1]))$fstatistic[1] %>% as.numeric})
    res_f <- rbind(res_f, as.data.frame(t(f)))
}

# calculate power
pwr <- do.call(rbind.data.frame, apply(res_p, 2, function(x) binom.test(sum(x < 5e-8), n_sim) %>% tidy)) %>% dplyr::select(estimate, conf.low, conf.high)
pwr$term <- lm(x ~ z1 + u + zn) %>% tidy %>% dplyr::pull(term)

# mean F-stat
f_stat <- do.call(rbind.data.frame, apply(res_f, 2, function(x) t.test(x) %>% tidy)) %>% dplyr::select(estimate, conf.low, conf.high)
# proportion of F < 10
f_stat_lt10 <- apply(res_f, 2, function(x) sum(x < 10))