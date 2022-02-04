library("dplyr")
library("broom")
library("varGWASR")
library("ggplot2")
source("funs.R")
set.seed(123)

n_sim <- 200
n_obs <- 1000

# simulate linear IV-exp variance effect and estimate MR point SD
results <- data.frame()
for (phi in seq(0, 5, 0.5)){
  for (i in 1:n_sim){
      c <- rnorm(n_obs)
      z <- get_simulated_genotypes(0.25, n_obs)
      u <- rbinom(n_obs, 2, 0.5)
      x <- z + c + rnorm(n_obs, sd=1+(z*phi))
      y <- x + c + rnorm(n_obs)
      exp_fit <- lm(x ~ z)
      b_exp <- exp_fit %>% tidy %>% dplyr::filter(term == "z") %>% dplyr::select("estimate") %>% as.numeric
      f_exp <- summary(exp_fit)$fstatistic[['value']]
      b_out <- lm(y ~ z) %>% tidy %>% dplyr::filter(term == "z") %>% dplyr::select("estimate") %>% as.numeric
      b_mr <- b_out/b_exp
      b_obs <- lm(y ~ x) %>% tidy %>% dplyr::filter(term == "x") %>% dplyr::select("estimate") %>% as.numeric
      fit <- varGWASR::model(data.frame(z,x), "z", "x")
      fit$b_mr <- b_mr
      fit$b_obs <- b_obs
      fit$lin_var <- var(x[z==2]) - var(x[z==0])
      fit$f_exp <- f_exp
      fit$phi <- phi
      results <- rbind(results, fit)
  }
}
results$iv1 <- results$f_exp > 10

# estimate MR SD per phi
tbl <- results %>% dplyr::filter(iv1) %>% dplyr::group_by(phi) %>% dplyr::summarize(phi_x1=mean(phi_x1), phi_x2=mean(phi_x2), sd=sd(b_mr))

# plot coef 1 vs 2 and then use MR point estimate SD
pdf("plot.pdf")
ggplot(tbl, aes(x=phi_x1, y=phi_x2, size=sd)) +
  geom_point() +
  theme_classic()
dev.off()