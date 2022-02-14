library("dplyr")
library("broom")
library("TwoSampleMR")
source("funs.R")
set.seed(123)

n_sim <- 50
n_obs <- 1000

# Paper
# 1. Bias and coverage of MR estimates under varying IV4 violation
# 2. Power and T1E to detect #1 using variance test
# 3. Applied example; for positives compare with IV strength among covariates
# 4. Hypothesis-free application?

# Hartwig et al does bias, coverage, rejection rate
# compare with IV strength among strata of measured covariates

# NOSH violation
results1 <- data.frame()
for (phi_z in seq(0, 2, 0.5)){
    for (phi_x in seq(0, 2, 0.5)){
        for (i in 1:n_sim){
            c <- rnorm(n_obs)
            z <- get_simulated_genotypes(0.25, n_obs)
            u <- rbinom(n_obs, 1, 0.5)
            x <- z + u + z*u*phi_z + c + rnorm(n_obs)
            y <- x + u + x*u*phi_x + c + rnorm(n_obs)
            res <- get_est(z, x, y)
            res$phi_z <- phi_z
            res$phi_x <- phi_x
            results1 <- rbind(results1, res)
        }
    }
}

# Not NOSH violated
results2 <- data.frame()
for (phi_z in seq(0, 2, 0.5)){
    for (phi_x in seq(0, 2, 0.5)){
        for (i in 1:n_sim){
            c <- rnorm(n_obs)
            z <- get_simulated_genotypes(0.25, n_obs)
            u1 <- rbinom(n_obs, 1, 0.5)
            u2 <- rbinom(n_obs, 1, 0.5)
            x <- z + u1 + z*u1*phi_z + c + rnorm(n_obs)
            y <- x + u2 + x*u2*phi_x + c + rnorm(n_obs)
            res <- get_est(z, x, y)
            res$phi_z <- phi_z
            res$phi_x <- phi_x
            results2 <- rbind(results2, res)
        }
    }
}

# combine results
results1$NOSH <- T
results2$NOSH <- F
results <- rbind(
    results1, results2
)

# ACE bias
mr_mean <- results %>% dplyr::group_by(NOSH,phi_z,phi_x) %>% dplyr::summarize(t.test(b_mr) %>% tidy)
pdf("bias.pdf")
ggplot(mr_mean, aes(x=phi_z, y=estimate, ymin=conf.low, ymax=conf.high, group=NOSH)) +
    geom_point() + geom_errorbar(width=.05) + theme_classic() + geom_hline(yintercept=1, linetype="dashed", color="grey") +
    labs(y="Mean causal estimate (95% CI [Monte Carlo])",x="IV-X interaction effect size") + facet_grid(~phi_x)
dev.off()