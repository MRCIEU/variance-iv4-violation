library("dplyr")
library("broom")
library("varGWASR")
library("ggplot2")
library("TwoSampleMR")
library("GWASTools")
source("funs.R")
set.seed(123)

n_sim <- 1000
n_obs <- 10000

# Paper
# 1. Bias and coverage of MR estimates under varying IV4 violation
# 2. Power and T1E to detect #1 using variance test
# 3. Applied example; for positives compare with IV strength among covariates
# 4. Hypothesis-free application?

# Hartwig et al does bias, coverage, rejection rate
# compare with IV strength among strata of measured covariates

# Z-X homogeneity
# X-Y homogeneity
results1 <- data.frame()
for (pb in seq(0.5, 0.5)){
  for (i in 1:n_sim){
      c <- rnorm(n_obs)
      z <- get_simulated_genotypes(0.25, n_obs)
      u <- rbinom(n_obs, 1, pb)
      x <- z + c + rnorm(n_obs)
      y <- x + c + rnorm(n_obs)
      result <- get_est(z, x, y)
      result$pb <- pb
      results1 <- rbind(results1, result)
  }
}

# Z-X heterogeneity i.e. sex-specific effects of SNP-exposure relationship (eg SLC2A9 x sex on urate)
# X-Y homogeneity (which is considered implausible)
# Z has variance effect on X and Y but is underpowered to detect Z-Y heterogeneity
results2 <- data.frame()
for (pb in seq(0.5, 0.5)){
  for (i in 1:n_sim){
      c <- rnorm(n_obs)
      z <- get_simulated_genotypes(0.25, n_obs)
      u <- rbinom(n_obs, 1, pb)
      x <- z + z*u + c + rnorm(n_obs)
      y <- x + c + rnorm(n_obs)
      result <- get_est(z, x, y)
      result$pb <- pb
      results2 <- rbind(results2, result)
  }
}

# Z-X homogeneity (is plausible)
# X-Y heterogeneity (is plausible)
# Z has no variance effect on X but does on Y but is underpowered to detect Z-Y heterogeneity
results3 <- data.frame()
for (pb in seq(0.5, 0.5)){
  for (i in 1:n_sim){
      c <- rnorm(n_obs)
      z <- get_simulated_genotypes(0.25, n_obs)
      u <- rbinom(n_obs, 1, pb)
      x <- z + c + rnorm(n_obs)
      y <- x*0.5 + x*u + c + rnorm(n_obs)
      result <- get_est(z, x, y)
      result$pb <- pb
      results3 <- rbind(results3, result)
  }
}

# Z-X heterogeneity (is plausible)
# X-Y heterogeneity (is plausible)
# Modifier is the same
# Z-X and Z-Y are strong variance effects
results4 <- data.frame()
for (pb in seq(0.5, 0.5)){
  for (i in 1:n_sim){
      c <- rnorm(n_obs)
      z <- get_simulated_genotypes(0.25, n_obs)
      u <- rbinom(n_obs, 1, pb)
      x <- z + z*u + c + rnorm(n_obs)
      y <- x*0.5 + x*u + c + rnorm(n_obs)
      result <- get_est(z, x, y)
      result$pb <- pb
      results4 <- rbind(results4, result)
  }
}

# Z-X heterogeneity (is plausible)
# X-Y heterogeneity (is plausible)
# Modifier is not the same
results5 <- data.frame()
for (pb in seq(0.5, 0.5)){
  for (i in 1:n_sim){
      c <- rnorm(n_obs)
      z <- get_simulated_genotypes(0.25, n_obs)
      u1 <- rbinom(n_obs, 1, pb)
      u2 <- rbinom(n_obs, 1, pb)
      x <- z + z*u1 + c + rnorm(n_obs)
      y <- x*.5 + x*u2 + c + rnorm(n_obs)
      result <- get_est(z, x, y)
      result$pb <- pb
      results5 <- rbind(results5, result)
  }
}

# combine results
results1$sim <- 1
results2$sim <- 2
results3$sim <- 3
results4$sim <- 4
results5$sim <- 5
results <- rbind(
    results1,results2,results3,results4,results5
)

# ACE bias
mr_mean <- results %>% dplyr::group_by(sim) %>% dplyr::summarize(t.test(b_mr) %>% tidy)
pdf("bias.pdf")
ggplot(mr_mean, aes(x=sim, y=estimate, ymin=conf.low, ymax=conf.high)) +
    geom_point() + geom_errorbar(width=.05) + theme_classic() + geom_hline(yintercept=1, linetype="dashed", color="grey") +
    labs(y="Average causal estimate (95% CI)",x="Scenario")
dev.off()

# ACE coverage
results$lci <- results$b_mr - (results$se_mr * 1.96)
results$uci <- results$b_mr + (results$se_mr * 1.96)
results$h1 <- (results$lci <= 1 & results$uci >= 1)
ci <- results %>% dplyr::group_by(sim) %>% dplyr::summarize(binom.test(sum(h1), n()) %>% tidy)

pdf("coverage.pdf")
ggplot(ci, aes(x=sim, y=estimate, ymin=conf.low, ymax=conf.high)) +
    geom_point() + geom_errorbar(width=.05) + theme_classic() + geom_hline(yintercept=0.95, linetype="dashed", color="grey") +
    labs(y="Causal estimate coverage (95% CI)",x="Scenario")
dev.off()

# T1E of variance test
save_qq <- function(p, file){
  pdf(file)
  qqPlot(p)
  dev.off()
}
save_qq(results %>% dplyr::filter(sim==1) %>% dplyr::pull(phi_p.zx), "zx_sim1.pdf")
save_qq(results %>% dplyr::filter(sim==1) %>% dplyr::pull(phi_p.zy), "zy_sim1.pdf")
save_qq(results %>% dplyr::filter(sim==2) %>% dplyr::pull(phi_p.zx), "zx_sim2.pdf")
save_qq(results %>% dplyr::filter(sim==2) %>% dplyr::pull(phi_p.zy), "zy_sim2.pdf")
save_qq(results %>% dplyr::filter(sim==3) %>% dplyr::pull(phi_p.zx), "zx_sim3.pdf")
save_qq(results %>% dplyr::filter(sim==3) %>% dplyr::pull(phi_p.zy), "zy_sim3.pdf")
save_qq(results %>% dplyr::filter(sim==4) %>% dplyr::pull(phi_p.zx), "zx_sim4.pdf")
save_qq(results %>% dplyr::filter(sim==4) %>% dplyr::pull(phi_p.zy), "zy_sim4.pdf")
save_qq(results %>% dplyr::filter(sim==5) %>% dplyr::pull(phi_p.zx), "zx_sim5.pdf")
save_qq(results %>% dplyr::filter(sim==5) %>% dplyr::pull(phi_p.zy), "zy_sim5.pdf")

# -log10 MR P value mean and SD
se <- results %>% dplyr::group_by(sim) %>% dplyr::summarize(t.test(se_mr) %>% tidy)
pdf("mr_se.pdf")
ggplot(se, aes(x=sim, y=estimate, ymin=conf.low, ymax=conf.high)) +
    geom_point() + geom_errorbar(width=.05) + theme_classic() + 
    labs(y="Causal estimate SE (95% CI)",x="Scenario")
dev.off()