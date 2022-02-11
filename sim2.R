library("dplyr")
library("TwoSampleMR")
library("genpwr")
source("funs.R")
set.seed(123)

n_sim <- 1000
n_obs <- 1000
n_snps <- 14

results <- data.frame()
for (i in 1:n_sim){
    # confounder
    c <- rnorm(n_obs)
    # modifier
    u <- rbinom(n_obs, 1, 0.5)
    # MAF of Z1
    z1_q <- runif(1, min=0.05, max=0.5)
    # MAF of Zn
    zn_q <- runif(n_snps, min=0.05, max=0.5)
    # Z1
    z1 <- get_simulated_genotypes(z1_q, n_obs)
    # Zn
    zn <- sapply(zn_q, function(q) get_simulated_genotypes(q, n_obs))   
    # calcualte SNP effect sizes
    b <- genpwr.calc(
        calc = "es", 
        model = "linear", 
        N = n_obs, 
        Power = 0.8, 
        MAF = c(z1_q, zn_q), 
        Alpha = 0.05, 
        sd_y = 1
    ) %>%
    dplyr::filter(Test.Model=="Additive", True.Model=="Additive") %>% dplyr::pull("ES_at_Alpha_0.05") %>% as.numeric
    # simulate exposure
    x <- z1*b[1] + rowSums(t(t(zn)*b[-1])) + c + rnorm(n_obs)
    # test power
    p <- lm(x~z1 + zn) %>% tidy %>% dplyr::pull(p.value)
    # store test statistics
    results <- rbind(results, as.data.frame(t(p)))
}

# calculate per-SNP power
apply(results, 2, function(x) binom.test(sum(x < 0.05), n_sim) %>% tidy)