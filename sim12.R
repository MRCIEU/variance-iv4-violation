library("dplyr")
library("ggplot2")
library("grid")
library("gtable")
library("TwoSampleMR")
library("viridis")
source("funs.R")
set.seed(123)

mr <- function(b_exp, b_out, se_exp, se_out, phi_p, n_obs, q){
    # define P-threshold using quantile
    p_var <- quantile(phi_p, q)

    # Drop top q SNPs with IV-exp variance effect
    b_exp <- b_exp[phi_p >= p_var]
    b_out <- b_out[phi_p >= p_var]
    se_exp <- se_exp[phi_p >= p_var]
    se_out <- se_out[phi_p >= p_var]

    # IVW
    ivw <- mr_ivw(b_exp, b_out, se_exp, se_out, NULL)
    ivw$q <- q
    ivw$n_obs <- n_obs

    return(ivw)
}

n_sim <- 10
r2_u <- 0.2 # U-Y main effect (small Cohen d)
r2_x <- 0.2 # X-Y main effect (small Cohen d)
r2_xu <- r2_x * 0.5 # X-Y interaction effect half the size of the main effect
n_snps <- 20 # number of SNPs for use as IVs in MR
n_isnps <- 5 # number of IVs with interaction effects on X
r2_z <- 0.05 # combined variance explained by the IV
phi <- 0.5 # interaction effect size relative to main effect
p_thresh <- 0.05 / n_snps
x_b <- sqrt(r2_x)
u_b <- sqrt(r2_u)
xu_b <- sqrt(r2_xu)

results <- data.frame()
f_stats <- data.frame()
for (n_obs in c(5000,10000,100000)){
    for (i in 1:n_sim){
        # select SNP betas
        if (n_obs == 5000){
            z_b <- rnorm(100000, sd=0.025)
            z_b <- sample(z_b[abs(z_b) > 0.085], n_snps)
        } else if (n_obs == 10000){
            z_b <- rnorm(100000, sd=0.025)
            z_b <- sample(z_b[abs(z_b) > 0.06], n_snps)
        } else if (n_obs == 100000){
            z_b <- rnorm(100000, sd=0.025)
            z_b <- sample(z_b[abs(z_b) > 0.005], n_snps)
        }
        
        # set size of interaction effect relative to main effect
        zu_b <- z_b * phi

        # simulate variables
        z <- sapply(1:n_snps, function(j) get_simulated_genotypes(0.25, n_obs)); z <- scale(z) # instruments
        u <- rnorm(n_obs) # modifier
        x <- sapply(1:n_obs, function(j) sum(z[j,]*z_b)) + u*u_b + sapply(1:n_obs, function(j) sum(z[j,1:n_isnps]*u[j]*zu_b[1:n_isnps])) + rnorm(n_obs, sd=sqrt(1-(sum(z_b^2) + u_b^2 + sum(zu_b[1:n_isnps]^2)))) # exposure
        y <- x*x_b + u*u_b + x*u*xu_b + rnorm(n_obs, sd=sqrt(1-(r2_x+r2_u+r2_xu))) # outcome
        varX <- var(x)
        varY <- var(y)

        # SNP estimates
        b_exp <- sapply(1:n_snps, function(i) lm(x~z[,i]) %>% tidy %>% dplyr::filter(term!="(Intercept)") %>% dplyr::pull(estimate))
        se_exp <- sapply(1:n_snps, function(i) lm(x~z[,i]) %>% tidy %>% dplyr::filter(term!="(Intercept)") %>% dplyr::pull(std.error))
        b_out <- sapply(1:n_snps, function(i) lm(y~z[,i]) %>% tidy %>% dplyr::filter(term!="(Intercept)") %>% dplyr::pull(estimate))
        se_out <- sapply(1:n_snps, function(i) lm(y~z[,i]) %>% tidy %>% dplyr::filter(term!="(Intercept)") %>% dplyr::pull(std.error))
        f_stat <- sapply(n_isnps+1:(n_snps-n_isnps), function(i) summary(lm(x~z[,i]))$fstatistic[1] %>% as.numeric)
        phi_p <- sapply(1:n_snps, function(i) model(data.frame(z, x), paste0("X", i), "x")$phi_p)

        # estimate IVW with subsets of instruemnts
        results <- rbind(results, mr(b_exp, b_out, se_exp, se_out, phi_p, n_obs, 0.75))
        results <- rbind(results, mr(b_exp, b_out, se_exp, se_out, phi_p, n_obs, 0.5))
        results <- rbind(results, mr(b_exp, b_out, se_exp, se_out, phi_p, n_obs, 0.25))
        results <- rbind(results, mr(b_exp, b_out, se_exp, se_out, phi_p, n_obs, 0.1))
        results <- rbind(results, mr(b_exp, b_out, se_exp, se_out, phi_p, n_obs, 0.05))
        results <- rbind(results, mr(b_exp, b_out, se_exp, se_out, phi_p, n_obs, 0))

        # store IV-exp f-stats
        f_stats <- rbind(f_stats, data.frame(f_stat, n_obs))
    }
}

# check f_stats
f_stats %>% dplyr::group_by(n_obs) %>% dplyr::summarize(t.test(f_stat) %>% tidy)

# plot bias
bias <- results %>% dplyr::group_by(n_obs, q) %>% dplyr::summarize(t.test(b) %>% tidy)
ggplot(bias, aes(x=q, y=estimate, ymin=conf.low, ymax=conf.high, shape=n_obs)) +
    geom_point() +
    geom_errorbar() +
    coord_flip()