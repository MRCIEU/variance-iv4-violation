library("dplyr")
library("TwoSampleMR")
library("genpwr")
source("funs.R")
set.seed(123)

# Model
# X=Z1+Z1U+U+Zn+C+E
# Y=X+U+XU+C+E
# Z-X has 80% power

n_sim <- 500
n_obs <- 1000
n_snps <- 14

# MAF of Z1
z1_q <- 0.25
# Main effect of Z1 on X
bx_z1 <- 0.3
# Main effect of U on X
bx_u <- 0.3
# interaction effect of Z1
bx_zu <- bx_z1*4

results <- data.frame()
for (i in 1:n_sim){
    # confounder
    c <- rnorm(n_obs)
    # modifier
    u <- rbinom(n_obs, 1, 0.5)
    # MAF of Zn
    zn_q <- runif(n_snps, min=0.05, max=0.5)
    # Z1
    z1 <- get_simulated_genotypes(z1_q, n_obs)
    # Zn
    zn <- sapply(zn_q, function(q) get_simulated_genotypes(q, n_obs))
    # calcualte SNP effect sizes
    pwr <- genpwr.calc(
        calc = "es", 
        model = "linear", 
        N = n_obs, 
        Power = 0.8,
        MAF = zn_q, 
        Alpha = 0.05, 
        sd_y = sqrt(
            (z1_q*2)^2*bx_z1^2 + (0.25*bx_u)^2 + (z1_q*2)^2*(0.5^2)*(bx_zu^2) + 1 + 1
        )
    )
    b <- pwr %>% dplyr::filter(Test.Model=="Additive", True.Model=="Additive") %>% dplyr::pull("ES_at_Alpha_0.05") %>% as.numeric
    # simulate exposure
    x <- z1*bx_z1 + u*bx_u + z1*u*bx_zu + rowSums(t(t(zn)*b)) + c + rnorm(n_obs)
    # test power
    p <- lm(x ~ z1*u + zn) %>% tidy %>% dplyr::pull(p.value)
    # store test statistics
    results <- rbind(results, as.data.frame(t(p)))
}

# calculate power
res <- do.call(rbind.data.frame, apply(results, 2, function(x) binom.test(sum(x < 0.05), n_sim) %>% tidy)) %>% dplyr::select(estimate, conf.low, conf.high)
res$term <- lm(x ~ z1*u + zn) %>% tidy %>% dplyr::pull(term)