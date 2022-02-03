library("dplyr")
library("broom")
library("varGWASR")
library("ggplot2")
set.seed(123)

n_sim <- 1000
n_obs <- 1000

#' Function to simulate genotypes in HWE
#' @param q Recessive/alternative allele frequency
#' @param n_obs Number of observations to return
get_simulated_genotypes <- function(q, n_obs){
  p <- 1 - q
  x <- sample(c(0, 1, 2), n_obs, prob=c(p^2, 2 * p * q, q^2), replace=T)
  return(x)
}

# simulate MR estimate with IV-exp interaction effect
results <- data.frame()
for (pb in seq(0, 1, 0.2)){
  for (i in 1:n_sim){
      c <- rnorm(n_obs)
      z <- get_simulated_genotypes(0.25, n_obs)
      u <- rbinom(n_obs, 2, pb)
      x <- z*u + c + rnorm(n_obs)
      y <- x + c + rnorm(n_obs)
      b_exp <- lm(x ~ z) %>% tidy %>% dplyr::filter(term == "z") %>% dplyr::select("estimate") %>% as.numeric
      b_out <- lm(y ~ z) %>% tidy %>% dplyr::filter(term == "z") %>% dplyr::select("estimate") %>% as.numeric
      b_mr <- b_exp / b_out
      b_obs <- lm(y ~ x) %>% tidy %>% dplyr::filter(term == "x") %>% dplyr::select("estimate") %>% as.numeric
      fit <- varGWASR::model(data.frame(z,x), "z", "x")
      fit$b_mr <- b_mr
      fit$b_obs <- b_obs
      fit$pb <- pb
      results <- rbind(results, fit)
  }
}

# plot
pdf("plot.pdf")
ggplot(results, aes(x=phi_p, y=b_mr)) +
  geom_point() +
  facet_wrap(~ pb, scale="free") +
  labs(x="Variance test P", y="MR point estimate") +
  theme_classic()
dev.off()